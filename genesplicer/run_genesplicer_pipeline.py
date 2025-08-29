import os
import subprocess
import argparse
import sys
import pandas as pd
import concurrent.futures
import multiprocessing
from pathlib import Path

script_dir = Path(__file__).parent.absolute()
dependencies_dir = script_dir.parents[2] / 'dependencies'
sys.path.insert(0, str(dependencies_dir))

from utility import read_fasta, get_mutation_data

def genesplicer_parallel(seq, genesplicer_dir, outdir, max_workers=None):
    """Run genesplicer in parallel for multiple sequences within a gene"""
    if not seq:  # Handle empty sequence dict
        return []

    if max_workers is None:
        # Limit workers to avoid overwhelming the system
        max_workers = min(multiprocessing.cpu_count() // 2, len(seq), 8)

    def run_single_genesplicer(item):
        k, v = item
        # Use unique temp file names to avoid conflicts
        tmp_file = os.path.join(genesplicer_dir, f'tmp_{os.getpid()}_{k}.fasta')
        try:
            with open(tmp_file, 'w') as f:
                f.write('>' + k + '\n' + v)
            cmd = f'cd {genesplicer_dir} && ./genesplicer {os.path.basename(tmp_file)} ../human -f {outdir}{k}-genesplicer.out'
            result = subprocess.getoutput(cmd)
            return k, result
        except Exception as e:
            return k, f"Error: {str(e)}"
        finally:
            # Clean up temp file
            if os.path.exists(tmp_file):
                try:
                    os.remove(tmp_file)
                except:
                    pass  # Ignore cleanup errors

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(run_single_genesplicer, seq.items()))

    return results


def genesplicer_sequential(seq, genesplicer_dir, outdir):
    """Sequential version as backup"""
    for k, v in seq.items():
        with open(os.path.join(genesplicer_dir, 'tmp.fasta'), 'w') as f:
            f.write('>' + k + '\n' + v)
        cmd = 'cd ' + genesplicer_dir + '&& ./genesplicer tmp.fasta ../human -f ' + outdir + k + '-genesplicer.out'
        subprocess.getoutput(cmd)


def run_genesplicer(directory, outdir, genesplicer_dir, is_mutant, mutation_dir=None, use_parallel=True):
    """Process genes sequentially, but mutations within each gene in parallel"""
    processed_genes = 0

    # Set default mutation directory if not provided
    if mutation_dir is None and not is_mutant:
        possible_paths = [
            os.path.join(os.path.dirname(directory), 'mut', 'genomic'),
            os.path.join(os.path.dirname(os.path.dirname(directory)), 'FASTA_files', 'mut', 'genomic')
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
        if file.is_file():
            processed_genes += 1
            print(f"Processing gene file {processed_genes}: {file.name}")

            try:
                if not is_mutant:
                    # Construct mutation file path
                    mut_filename = file.name.rsplit('.')[0] + '_all_muts_nt.fasta'
                    mut_filepath = os.path.join(mutation_dir, mut_filename)

                    if not os.path.exists(mut_filepath):
                        print(f"  Warning: Mutation file not found: {mut_filepath}")
                        continue

                    mutseq = read_fasta(mut_filepath)
                    refseq = read_fasta(os.path.join(directory, file.name))

                    # Build complete sequence dict for this gene
                    seq = {}
                    for k in mutseq.keys():
                        seq[k] = refseq['genomic']

                    # Use parallel processing for mutations within this gene
                    if use_parallel and len(seq) > 1:
                        print(f"  Processing {len(seq)} mutations in parallel...")
                        genesplicer_parallel(seq, genesplicer_dir, outdir)
                    else:
                        print(f"  Processing {len(seq)} mutations sequentially...")
                        genesplicer_sequential(seq, genesplicer_dir, outdir)
                else:
                    seq = read_fasta(os.path.join(directory, file.name))
                    if use_parallel and len(seq) > 1:
                        print(f"  Processing {len(seq)} sequences in parallel...")
                        genesplicer_parallel(seq, genesplicer_dir, outdir)
                    else:
                        print(f"  Processing {len(seq)} sequences sequentially...")
                        genesplicer_sequential(seq, genesplicer_dir, outdir)

            except Exception as e:
                print(f"Error processing {file.name}: {str(e)}")
                with open('gene_processing_errors.log', 'a') as f:
                    f.write(f'Error processing gene file {file.name}: {str(e)}\n')

    print(f"Completed processing {processed_genes} gene files")


def parse_genesplicer(outdir, mutation_data_dir=None):
    """Parse genesplicer output files and return filtered results"""
    all_results = []
    csv_cache = {}  # Cache CSV files to avoid repeated reads
    processed_files = 0

    # Set default mutation data directory if not provided
    if mutation_data_dir is None:
        possible_paths = [
            os.path.join(os.getcwd(), 'mutations', 'combined'),
            os.path.join(os.path.dirname(os.getcwd()), 'mutations', 'combined'),
            '../mutations/combined',
            './mutations/combined'
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

    print("Starting to parse GeneSplicer output files...")

    for file in os.scandir(outdir):
        if file.is_file():
            processed_files += 1
            if processed_files % 1000 == 0:  # Progress indicator
                print(f"Parsed {processed_files} files...")

            try:
                genename = file.name.rsplit('-')[0]
                ntposnt = file.name.rsplit('-')[1]
                genomic_pos, genomic_mut = get_mutation_data(ntposnt)

                # Cache CSV reads to avoid repeated file I/O
                csv_file = os.path.join(mutation_data_dir, f'combined_{genename}.csv')

                if not os.path.exists(csv_file):
                    print(f"  Warning: CSV file not found: {csv_file}")
                    continue

                if csv_file not in csv_cache:
                    csv_cache[csv_file] = pd.read_csv(csv_file, sep='\t')
                ref = csv_cache[csv_file]

                match_ref_ntposnt = ref[(ref['genomic'] == ntposnt)]['mutant']
                if match_ref_ntposnt.empty:
                    continue  # Skip if no match found

                match_ref_ntposnt = match_ref_ntposnt.reset_index(drop=True)[0]

                # Read and process data
                d = pd.read_csv(os.path.join(outdir, file.name), sep=' ', header=None)
                d = d.rename(columns={0: 'End5', 1: 'End3', 2: 'Score', 3: 'confidence', 4: 'splice_site_type'})

                # Filter for matches
                matched_donor_rows = d[(d['End5'] == genomic_pos) & (d['splice_site_type'] == 'donor')]
                matched_acceptor_rows = d[(d['End3'] == genomic_pos) & (d['splice_site_type'] == 'acceptor')]
                matched_data = pd.concat([matched_donor_rows, matched_acceptor_rows], ignore_index=True)

                # Only process if we have matches
                if not matched_data.empty:
                    matched_data = matched_data.assign(
                        **{'Gene Name': genename + '-' + match_ref_ntposnt, 'Wt nt': genomic_mut[0],
                           'Genomic Position': genomic_pos, 'Mut nt': genomic_mut[1]}
                    )
                    matched_data = matched_data[['Gene Name', 'Wt nt', 'Position', 'Mut nt'] +
                                                [col for col in matched_data.columns if
                                                 col not in ['Gene Name', 'Wt nt', 'Position', 'Mut nt']]]

                    all_results.append(matched_data)

            except Exception as e:
                with open('parse_errors.log', 'a') as f:
                    f.write(f'Error parsing file {file.name}: {str(e)}\n')

    print(f"Completed parsing {processed_files} files")
    return all_results


def save_results(all_results, is_mutant):
    """Save parsed results to file"""
    if is_mutant:
        output_file = 'filtered_mut_genesplicer.tsv'
    else:
        output_file = 'filtered_wt_genesplicer.tsv'

    if all_results:
        print(f"Processing {len(all_results)} files with matches...")
        final_data = pd.concat(all_results, ignore_index=True)
        final_data.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(final_data)} total records to {output_file}")
        return final_data
    else:
        print("No matches found across all files")
        return pd.DataFrame()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="GeneSplicer processing pipeline with parallelization")
    parser.add_argument("-i", "--input", help="input directory of the fasta files (use full path)")
    parser.add_argument("-o", "--output", help="output directory of GeneSplicer", required=True)
    parser.add_argument("-g", "--genesplicer_dir", help="path to directory containing the genesplicer executable")
    parser.add_argument('-m', '--is_mutant', action='store_true', help='process mutant sequences (default: wildtype)')
    parser.add_argument("-p", "--pipeline", help="which pipeline to run: full/GS/parse_GS",
                        choices=['full', 'GS', 'parse_GS'], required=True)
    parser.add_argument("--no-parallel", action='store_true', help="disable parallel processing")
    parser.add_argument("--mutation_dir", help="directory containing mutation FASTA files (for wildtype processing)")
    parser.add_argument("--mutation_data_dir", help="directory containing mutation CSV files (for parsing)")

    args = parser.parse_args()

    # Validate required arguments based on pipeline
    if args.pipeline in ['full', 'GS']:
        if not args.input:
            parser.error("--input is required for 'full' and 'GS' pipelines")
        if not args.genesplicer_dir:
            parser.error("--genesplicer_dir is required for 'full' and 'GS' pipelines")

    # Determine if using parallel processing
    use_parallel = not args.no_parallel

    print(f"Starting pipeline: {args.pipeline}")
    if args.pipeline in ['full', 'GS']:
        print(f"Processing {'mutant' if args.is_mutant else 'wildtype'} sequences")
        print(f"Parallel processing: {'enabled' if use_parallel else 'disabled'}")

    if args.pipeline == 'full':
        print("Running full pipeline: GeneSplicer + Parsing")
        run_genesplicer(args.input, args.output, args.genesplicer_dir, args.is_mutant,
                        args.mutation_dir, use_parallel)
        results = parse_genesplicer(args.output, args.mutation_data_dir)
        save_results(results, args.is_mutant)
    elif args.pipeline == 'GS':
        print("Running GeneSplicer only")
        run_genesplicer(args.input, args.output, args.genesplicer_dir, args.is_mutant,
                        args.mutation_dir, use_parallel)
    else:  # parse_GS
        print("Running parsing only")
        results = parse_genesplicer(args.output, args.mutation_data_dir)
        save_results(results, args.is_mutant)

    print("Pipeline completed successfully!")