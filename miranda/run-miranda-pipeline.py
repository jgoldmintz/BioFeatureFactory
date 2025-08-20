import os
import subprocess
import argparse
import re
import pandas as pd
import concurrent.futures
import multiprocessing
#from dependencies.tools import get_mutation_data

# Global flag for cleanup
cleanup_raw_files = False

def read_fasta(inf, aformat="FIRST", duplicate="replace"):
    data = {}
    with open(inf, "r") as fa:
        name = ""
        for line in fa.readlines():
            if "#" in line:
                continue
            if ">" in line:
                if aformat.upper() == "NCBI":
                    name = re.search(">[a-zA-Z]+_?\d+(\.\d+)*", line).group(0)
                elif aformat.upper() in ["FIRST", "WORD"]:
                    name = line.split()[0]
                else:
                    name = line.strip()
                name = name[1:].strip()
                if name in data.keys():
                    if duplicate.lower() in ["append", "a"]:  # simply add to existing sequence
                        pass
                    elif duplicate.lower() in ["replace", "r"]:  # reset sequence to empty
                        data[name] = ""
                    elif duplicate.lower() in ["separate", "s"]:  # add underscore+number to end of sequence name
                        matches = re.findall("/_\d+$/", name)
                        if matches != None and len(matches) > 0:
                            num = int(max(matches)[1:])
                            name = name[:-len(str(num))] + str(num + 1)
                            data[name] = ""
                        else:
                            name = name + "_2"
                            data[name] = ""
                else:
                    data[name] = ""
            else:
                data[name] = data[name] + line.strip()
    return data

def get_mutation_data(ntposnt):
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1]) - 1  # Convert to 0-based index
    return position, (original_nt, mutant_nt)

def run_single_miranda(args):
    """Run miranda for a single sequence - defined at module level for pickling"""
    k, v, miranda_dir, outdir, mirna_db = args

    # Use unique temp file names to avoid conflicts
    tmp_file = os.path.join(miranda_dir, f'tmp_{os.getpid()}_{k}.fasta')
    out_file = os.path.join(outdir, f'{k}-miranda.out')

    # Skip if output already exists (for resuming interrupted runs)
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return k, "Already processed"

    try:
        with open(tmp_file, 'w') as f:
            f.write('>' + k + '\n' + v)

        # Use subprocess.run for better control
        cmd = ['./miranda', mirna_db, os.path.basename(tmp_file)]
        #print(f"DEBUG: Running command: {' '.join(cmd)} in directory {miranda_dir}")
        #print(f"DEBUG: Temp file exists: {os.path.exists(tmp_file)}")
        result = subprocess.run(
            cmd,
            cwd=miranda_dir,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per miranda run
        )

        #print(f"DEBUG: Return code: {result.returncode}")
        #print(f"DEBUG: STDOUT length: {len(result.stdout)}")
        #print(f"DEBUG: STDERR: {result.stderr}")
        
        # Write output directly
        with open(out_file, 'w') as f:
            f.write(result.stdout)

        return k, "Success"
    except subprocess.TimeoutExpired:
        return k, "Error: Timeout"
    except Exception as e:
        return k, f"Error: {str(e)}"
    finally:
        # Clean up temp file
        if os.path.exists(tmp_file):
            try:
                os.remove(tmp_file)
            except:
                pass


def miranda_parallel(seq, miranda_dir, outdir, mirna_db, max_workers=None):
    """Run miranda in parallel for multiple sequences within a gene"""
    if not seq:  # Handle empty sequence dict
        return []

    if max_workers is None:
        # Use more workers for CPU-bound tasks
        max_workers = min(multiprocessing.cpu_count() - 1, len(seq), 16)

    print(f"    Using {max_workers} workers for {len(seq)} sequences")

    # Prepare arguments for each task
    task_args = [(k, v, miranda_dir, outdir, mirna_db) for k, v in seq.items()]

    # Use ProcessPoolExecutor for CPU-bound tasks
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Process in batches to show progress
        batch_size = 100
        all_results = []

        for i in range(0, len(task_args), batch_size):
            batch = task_args[i:i + batch_size]
            results = list(executor.map(run_single_miranda, batch))
            all_results.extend(results)

            # Show progress
            completed = min(i + batch_size, len(task_args))
            print(f"      Processed {completed}/{len(task_args)} sequences")

            # Log any errors
            for k, status in results:
                if status.startswith("Error"):
                    print(f"        {k}: {status}")

    return all_results


def miranda_sequential(seq, miranda_dir, outdir, mirna_db):
    """Sequential version as backup"""
    for k, v in seq.items():
        out_file = os.path.join(outdir, f'{k}-miranda.out')

        # Skip if already processed
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            print(f"    Skipping {k} - already processed")
            continue

        tmp_file = os.path.join(miranda_dir, 'tmp.fasta')
        with open(tmp_file, 'w') as f:
            f.write('>' + k + '\n' + v)
        
        # Use same approach as parallel version
        cmd = ['./miranda', mirna_db, 'tmp.fasta']
        result = subprocess.run(
            cmd,
            cwd=miranda_dir,
            capture_output=True,
            text=True,
            timeout=300
        )
        
        # Write output directly
        with open(out_file, 'w') as f:
            f.write(result.stdout)


def run_miranda(directory, outdir, miranda_dir, is_mutant, mirna_db, mutation_dir=None, use_parallel=True,
                max_workers=None):
    """Process genes sequentially, but mutations within each gene in parallel"""
    processed_genes = 0

    # Create output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # Set default mutation directory if not provided
    if mutation_dir is None and not is_mutant:
        possible_paths = [
            os.path.join(os.path.dirname(directory), 'mut', 'transcript'),
            os.path.join(os.path.dirname(os.path.dirname(directory)), 'FASTA_files', 'mut', 'transcript'),
            '/Users/jgoldmintz/disease_associated_synonymous_mutations/FASTA_files/mut/transcript/'
        ]

        for path in possible_paths:
            if os.path.exists(path):
                mutation_dir = path
                print(f"Found mutation directory: {mutation_dir}")
                break

        if mutation_dir is None:
            raise FileNotFoundError(
                "Could not find mutation directory. Please specify --mutation_dir explicitly. "
                f"Searched: {possible_paths}"
            )

    for file in os.scandir(directory):
        if file.is_file() and file.name.endswith('.fasta'):
            processed_genes += 1
            print(f"\nProcessing gene file {processed_genes}: {file.name}")

            try:
                if not is_mutant:
                    # For wildtype mode - run miranda ONCE and duplicate results
                    mut_filename = file.name.rsplit('.')[0] + '_all_muts_nt.fasta'
                    mut_filepath = os.path.join(mutation_dir, mut_filename)

                    if not os.path.exists(mut_filepath):
                        print(f"  Warning: Mutation file not found: {mut_filepath}")
                        continue

                    # Load mutation file ONLY to get mutation IDs (not using the sequences!)
                    mutation_data = read_fasta(mut_filepath)
                    mutation_ids = list(mutation_data.keys())  # We only need the IDs

                    # Load reference sequence - this is what we'll actually run miranda on
                    refseq = read_fasta(os.path.join(directory, file.name))

                    if not mutation_ids or not refseq or 'transcript' not in refseq:
                        print(f"  No sequences found in files")
                        continue

                    gene_name = file.name.rsplit('.')[0]
                    print(f"  Found {len(mutation_ids)} mutation IDs for gene {gene_name}")

                    # Check if we need to run miranda (skip if all output files exist)
                    files_to_create = []
                    for mutation_id in mutation_ids:
                        out_file = os.path.join(outdir, f'{mutation_id}-miranda.out')
                        if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                            files_to_create.append((mutation_id, out_file))

                    if not files_to_create:
                        print(f"  All output files already exist, skipping...")
                        continue

                    print(f"  Need to create {len(files_to_create)} output files")
                    print(f"  Running miranda ONCE on wildtype sequence...")

                    # Run miranda once on the REFERENCE sequence
                    tmp_file = os.path.join(miranda_dir, f'tmp_{gene_name}_wt.fasta')

                    # Write REFERENCE sequence to temp file
                    with open(tmp_file, 'w') as f:
                        f.write('>transcript\n' + refseq['transcript'])

                    # Run miranda once on REFERENCE sequence
                    cmd = ['./miranda', mirna_db, tmp_file]
                    result = subprocess.run(
                        cmd,
                        cwd=miranda_dir,
                        capture_output=True,
                        text=True,
                        timeout=300
                    )
                    miranda_output = result.stdout

                    # Clean up temp file
                    if os.path.exists(tmp_file):
                        os.remove(tmp_file)

                    # Now create output files for each mutation ID
                    print(f"  Creating {len(files_to_create)} wildtype output files...")
                    for mutation_id, out_file in files_to_create:
                        with open(out_file, 'w') as f:
                            f.write(miranda_output)

                    print(f"    Done! Created {len(files_to_create)} output files")
                else:
                    # For mutant mode - process sequences as-is
                    seq = read_fasta(os.path.join(directory, file.name))

                    if not seq:
                        print(f"  No sequences found in file")
                        continue

                    print(f"  Found {len(seq)} sequences to process")

                    # Warning for large sequence counts
                    if len(seq) > 1000:
                        print(f"  WARNING: Processing {len(seq)} sequences will create {len(seq)} miranda output files.")
                        response = input("  Do you want to automatically delete raw miranda files after parsing? (y/n): ").lower().strip()
                        global cleanup_raw_files
                        cleanup_raw_files = response in ['y', 'yes']
                        if cleanup_raw_files:
                            print("  Raw miranda files will be deleted after parsing to save disk space.")
                        else:
                            print("  Raw miranda files will be kept.")

                    if use_parallel:
                        print(f"  Processing in parallel...")
                        miranda_parallel(seq, miranda_dir, outdir, mirna_db, max_workers)
                    else:
                        print(f"  Processing sequentially...")
                        miranda_sequential(seq, miranda_dir, outdir, mirna_db)

            except Exception as e:
                print(f"  Error processing {file.name}: {str(e)}")
                with open('gene_processing_errors.log', 'a') as f:
                    f.write(f'Error processing gene file {file.name}: {str(e)}\n')

    print(f"\nCompleted processing {processed_genes} gene files")

def parse_miranda_file(miranda_file_path, transcript_pos, pkey):
    """Parse a single miranda output file directly in Python - follows awk script logic"""
    matched_rows = []
    
    try:
        with open(miranda_file_path, 'r') as f:
            lines = f.readlines()
        
        # State variables (like awk script)
        current_mirna = ""
        query_seq = ""
        ref_seq = ""
        found_hit = False
        
        for line in lines:
            line = line.strip()
            
            # Start of a new miRNA block
            if line.startswith('Read Sequence:') and 'hsa-' in line:
                # Reset for new miRNA
                current_mirna = line.split('Read Sequence:')[1].split()[0]
                query_seq = ""
                ref_seq = ""
                found_hit = False
            
            # Capture query sequence (entire line after "Query:")
            elif 'Query:' in line:
                query_seq = line[line.index('Query:') + 6:].strip()  # Everything after "Query:"
                # Clean up: extract sequence between quotes and remove 3'/5' markers
                if "'" in query_seq:
                    parts = query_seq.split("'")
                    if len(parts) >= 2:
                        seq = parts[1].strip()
                        query_seq = seq.replace("3", "").replace("5", "").replace(" ", "")
            
            # Capture reference sequence (entire line after "Ref:")
            elif 'Ref:' in line:
                ref_seq = line[line.index('Ref:') + 4:].strip()  # Everything after "Ref:"
                # Clean up: extract sequence between quotes and remove 3'/5' markers
                if "'" in ref_seq:
                    parts = ref_seq.split("'")
                    if len(parts) >= 2:
                        seq = parts[1].strip()
                        ref_seq = seq.replace("3", "").replace("5", "").replace(" ", "")
            
            # Capture hit data from ">>" lines
            elif line.startswith('>>'):
                parts = line.split('\t')
                if len(parts) >= 10:
                    try:
                        # Extract data from the hit line
                        tot_score = parts[2].strip()
                        tot_energy = parts[3].strip()
                        max_score = parts[4].strip()
                        max_energy = parts[5].strip()
                        strand = parts[6].strip()
                        len1 = parts[7].strip()
                        len2 = parts[8].strip()
                        positions_str = parts[9].strip()
                        
                        # Parse positions
                        positions = []
                        for pos_part in positions_str.split():
                            try:
                                positions.append(int(pos_part.strip()))
                            except ValueError:
                                continue
                        
                        # Check if our target transcript position is in the positions
                        if transcript_pos in positions:
                            matched_row = {
                                'Pkey': pkey,
                                'miRNA': current_mirna,
                                'Transcript-Position': transcript_pos,
                                'Tot Score': tot_score,
                                'Tot Energy': tot_energy,
                                'Max Score': max_score,
                                'Max Energy': max_energy,
                                'Strand': strand,
                                'Len1': len1,
                                'Len2': len2,
                                'query-sequence': query_seq,
                                'ref-sequence': ref_seq
                            }
                            matched_rows.append(matched_row)
                            found_hit = True
                            
                    except (IndexError, ValueError):
                        continue  # Skip malformed lines
    
    except Exception as e:
        print(f"  Error parsing file {miranda_file_path}: {e}")
    
    return matched_rows


def parse_miranda(outdir, mutation_data_dir=None):
    """Parse miranda output files and return filtered results"""
    all_results = []
    csv_cache = {}  # Cache CSV files to avoid repeated reads
    processed_files = 0

    # Set default mutation data directory if not provided
    if mutation_data_dir is None:
        possible_paths = [
            os.path.join(os.getcwd(), 'mutations', 'combined', 'transcript'),
            os.path.join(os.path.dirname(os.getcwd()), 'mutations', 'combined', 'transcript'),
            './mutations/combined/transcript',
            '../mutations/combined/transcript'
        ]

        for path in possible_paths:
            if os.path.exists(path):
                mutation_data_dir = path
                print(f"Found mutation data directory: {mutation_data_dir}")
                break

        if mutation_data_dir is None:
            raise FileNotFoundError(
                "Could not find mutation data directory. Please specify --mutation_data_dir explicitly. "
                f"Searched: {possible_paths}"
            )

    print("Starting to parse miranda output files...")

    for file in os.scandir(outdir):
        if file.is_file() and file.name.endswith('-miranda.out'):
            processed_files += 1
            if processed_files % 1000 == 0:  # Progress indicator
                print(f"Parsed {processed_files} files...")

            try:
                # Extract components from filename
                filename_base = file.name.replace('-miranda.out', '')
                parts = filename_base.split('-')
                genename = parts[0]
                file_ntposnt = parts[1]

                # Load transcript mutation mapping
                csv_file = os.path.join(mutation_data_dir, f'combined_{genename}.csv')
                
                if not os.path.exists(csv_file):
                    print(f"  Warning: CSV file not found: {csv_file}")
                    continue

                if csv_file not in csv_cache:
                    csv_cache[csv_file] = pd.read_csv(csv_file, sep='\t')
                transcript_mapping = csv_cache[csv_file]

                # File contains transcript mutation, get corresponding mutant for pkey
                mutant_match = transcript_mapping[transcript_mapping['transcript'] == file_ntposnt]['mutant']
                if not mutant_match.empty:
                    mutant_ntposnt = mutant_match.reset_index(drop=True)[0]
                    transcript_ntposnt = file_ntposnt
                    pkey = f"{genename}-{mutant_ntposnt}"
                else:
                    print(f"  No mapping found for transcript {file_ntposnt}")
                    continue
                
                # Extract transcript position using get_mutation_data
                transcript_pos, transcript_mut = get_mutation_data(transcript_ntposnt)

                # Parse miranda output directly in Python (no subprocess)
                matched_rows = parse_miranda_file(file.path, transcript_pos, pkey)
                
                # Convert to DataFrame and add to results
                if matched_rows:
                    matched_df = pd.DataFrame(matched_rows)
                    all_results.append(matched_df)

            except Exception as e:
                with open('parse_errors.log', 'a') as f:
                    f.write(f'Error parsing file {file.name}: {str(e)}\n')

    print(f"Completed parsing {processed_files} files")
    return all_results


def save_results(all_results, is_mutant, outdir=None):
    """Save parsed results to file"""
    if is_mutant:
        output_file = 'filtered_mut_miranda.tsv'
    else:
        output_file = 'filtered_wt_miranda.tsv'

    if all_results:
        print(f"Processing {len(all_results)} files with matches...")
        final_data = pd.concat(all_results, ignore_index=True)
        final_data.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(final_data)} total records to {output_file}")

        # Cleanup raw miranda files if requested
        global cleanup_raw_files
        if cleanup_raw_files and outdir:
            print("Cleaning up raw miranda output files...")
            miranda_files = [f for f in os.listdir(outdir) if f.endswith('-miranda.out')]
            deleted_count = 0
            for file in miranda_files:
                try:
                    os.remove(os.path.join(outdir, file))
                    deleted_count += 1
                except:
                    pass
            print(f"Deleted {deleted_count} raw miranda files to save disk space.")

        return final_data
    else:
        print("No matches found across all files")
        return pd.DataFrame()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="miranda processing pipeline with parallelization")
    
    # Pipeline mode selection
    parser.add_argument("--mode", choices=['full-pipeline', 'miranda-only', 'parse-only'], 
                        default='full-pipeline',
                        help="Pipeline mode: full-pipeline (run miranda + parse), miranda-only (just run miranda), parse-only (just parse existing output)")
    
    # Input/Output arguments
    parser.add_argument("-i", "--input", help="input directory of the fasta files (use full path)")
    parser.add_argument("-o", "--output", help="output directory of miranda")
    
    # Miranda-specific arguments
    parser.add_argument("-m", "--miranda_dir", help="path to directory containing the miranda executable")
    parser.add_argument("-d", "--mirna_db", help="path to mirna database")
    parser.add_argument('-M', '--is_mutant', action='store_true',
                        help='process mutant sequences (default: wildtype)')
    
    # Parsing-specific arguments
    parser.add_argument("--mutation_data_dir", help="path to transcript mutation data directory (default: auto-detect)")
    
    # Performance arguments
    parser.add_argument("--no-parallel", action='store_true', help="disable parallel processing")
    parser.add_argument("--max-workers", type=int, help="maximum number of parallel workers")

    args = parser.parse_args()
    
    # Validate required arguments based on mode
    if args.mode in ['full-pipeline', 'miranda-only']:
        if not all([args.input, args.output, args.miranda_dir, args.mirna_db]):
            parser.error("For miranda execution modes, -i, -o, -m, and -d are required")
    
    if args.mode in ['full-pipeline', 'parse-only']:
        if not args.output:
            parser.error("For parsing modes, -o (output directory) is required")
    
    use_parallel = not args.no_parallel
    max_workers = args.max_workers if args.max_workers else None

    print(f"Starting miranda pipeline in '{args.mode}' mode")
    print(f"Processing {'mutant' if args.is_mutant else 'wildtype'} sequences")
    
    if args.mode in ['full-pipeline', 'miranda-only']:
        print(f"Parallel processing: {'enabled' if use_parallel else 'disabled'}")
        if max_workers and use_parallel:
            print(f"Max workers set to: {max_workers}")

        # Convert miRNA database to absolute path
        mirna_db_abs = os.path.abspath(args.mirna_db)
        
        # Run miranda
        run_miranda(
            args.input,
            args.output,
            args.miranda_dir,
            args.is_mutant,
            mirna_db_abs,
            use_parallel=use_parallel,
            max_workers=max_workers
        )
        print("Miranda execution completed!")

    if args.mode in ['full-pipeline', 'parse-only']:
        print("\nStarting miranda output parsing...")
        
        # Parse miranda results
        results = parse_miranda(args.output, args.mutation_data_dir)
        
        if results:
            # Save results
            final_data = save_results(results, args.is_mutant, args.output)
            print(f"Parsing completed! Saved {len(final_data)} total records")
        else:
            print("No matching results found during parsing")

    print(f"\n{args.mode} pipeline completed successfully!")