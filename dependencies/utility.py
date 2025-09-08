import inspect
import warnings
import re
import os
import math
import tempfile
import subprocess
import shutil

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

def trim_muts(ntPosnt):
    mut_list = []
    with open(ntPosnt, 'r') as inf:
        for i in inf:
            mut_list.append(i.replace('*',''))

        return [i.replace('\n','') for i in mut_list[1:len(mut_list)]]

def get_mutation_data(ntposnt):
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1]) - 1  # Convert to 0-based index
    return position, (original_nt, mutant_nt)

def get_mutation_data_bioAccurate(ntposnt):
    original_nt = ntposnt[0]
    mutant_nt = ntposnt[-1]
    if int(ntposnt[1:-1]) == 1:
        position = int(ntposnt[1:-1])
    else:
        position = int(ntposnt[1:-1])
    return position, (original_nt, mutant_nt)

def convert_position(seq1, seq2, position1, space="-"):
    error = None

    if position1 == 0:
        warnings.warn("\033[93m" + "Position given is not 1-indexed in function " + str(inspect.stack()[1].function) + "\033[0m")
        error = "Position given is not 1-indexed\n" + str(inspect.stack()[1].function)

    i1 = 0
    i2 = 0
    increment = 0
    while i1 < int(position1) and increment < len(seq1) and increment < len(seq2):
        if seq1[increment] != space:
            i1 += 1
        if seq2[increment] != space:
            i2 += 1
        increment += 1
    if not seq1[increment - 1] == seq1.replace(space, "")[i1 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq1[increment - 1] + " " + seq1.replace(space, "")[i1 - 1] + "\033[0m")
    elif seq2[increment - 1] != space and not seq2[increment - 1] == seq2.replace(space, "")[i2 - 1]:
        warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq2[increment - 1] + " " + seq2.replace(space, "")[i2 - 1] + "\033[0m")
    if seq2[increment - 1] == space:
        error = "Sequence 1 position aligns with a gap in sequence 2\n" + str(inspect.stack()[1].function)
    return (i2, error)

def split_fasta_into_batches(fasta_file, batch_size=100):
    """Split a FASTA file into smaller batches for processing
    
    Args:
        fasta_file: Path to input FASTA file
        batch_size: Number of sequences per batch (default: 100)
        
    Returns:
        List of batch file paths created
    """
    try:
        # Read all sequences from FASTA file
        sequences = read_fasta(fasta_file)
        total_sequences = len(sequences)
        
        if total_sequences == 0:
            return []
        
        print(f"Splitting {total_sequences} sequences into batches of {batch_size}")
        
        # Calculate number of batches needed
        num_batches = math.ceil(total_sequences / batch_size)
        
        batch_files = []
        sequence_items = list(sequences.items())
        
        for i in range(num_batches):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, total_sequences)
            
            # Create batch sequences
            batch_sequences = dict(sequence_items[start_idx:end_idx])
            batch_count = len(batch_sequences)
            
            # Create temporary batch file
            batch_filename = fasta_file.replace('.fasta', f'_batch{i+1}.fasta')
            
            with open(batch_filename, 'w') as f:
                for seq_name, sequence in batch_sequences.items():
                    f.write(f">{seq_name}\n")
                    # Write sequence in lines of 80 characters
                    for j in range(0, len(sequence), 80):
                        f.write(sequence[j:j+80] + "\n")
            
            batch_files.append(batch_filename)
            print(f"Created batch {i+1}/{num_batches}: {batch_filename} ({batch_count} sequences)")
        
        return batch_files
        
    except Exception as e:
        print(f"Error splitting FASTA file {fasta_file}: {e}")
        return []

def combine_batch_outputs(batch_output_files, final_output_file, format_type='netnglyc'):
    """
    Combine multiple batch outputs into a single file for backwards compatibility
    
    Args:
        batch_output_files: List of individual batch output files
        final_output_file: Path for combined output file
        format_type: Output format type ('netnglyc' or 'netphos')
        
    Returns:
        bool: True if combination successful
    """
    try:
        if not batch_output_files:
            return False
        
        print(f"Combining {len(batch_output_files)} batch outputs (format: {format_type})...")
        
        if format_type == 'netnglyc':
            return _combine_glycosylation_outputs(batch_output_files, final_output_file)
        elif format_type == 'netphos':
            return _combine_phosphorylation_outputs(batch_output_files, final_output_file)
        else:
            raise ValueError(f"Unsupported format_type: {format_type}")
            
    except Exception as e:
        print(f"Error combining batch outputs: {e}")
        return False

def _combine_glycosylation_outputs(batch_output_files, final_output_file):
    """Combine glycosylation prediction batch outputs (NetNGlyc format)"""
    try:
        # Count total sequences for header
        total_sequences = 0
        all_sequence_sections = []
        all_prediction_lines = []
        
        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()
                
                # Extract sequences count from header
                lines = content.split('\n')
                for line in lines:
                    if 'amino acids' in line and line.startswith('>'):
                        try:
                            seq_count = int(line.split()[1])
                            total_sequences += seq_count
                        except:
                            total_sequences += 1  # Fallback
                        break
                
                # Collect sequence display sections and prediction lines
                in_sequence_section = False
                in_prediction_section = False
                
                for line in lines:
                    if line.strip().startswith('>') and 'amino acids' in line:
                        continue  # Skip individual headers
                    
                    if 'SeqName' in line and 'Position' in line:
                        in_prediction_section = True
                        continue
                    
                    if line.startswith('    ') and len(line.strip()) > 10:
                        in_sequence_section = True
                        all_sequence_sections.append(line)
                    elif in_prediction_section and line.strip() and not line.startswith('#'):
                        all_prediction_lines.append(line)
                        
            except Exception as e:
                print(f"Error reading batch file {batch_file}: {e}")
                continue
        
        # Write combined output
        with open(final_output_file, 'w') as f:
            # Write header with total sequence count
            f.write(f">{os.path.basename(final_output_file).replace('.out', '')}\t{total_sequences} amino acids\n\n")
            
            # Write prediction header
            f.write("SeqName                 Position  Potential  N-Glyc result  Comment\n")
            f.write("=" * 70 + "\n")
            
            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")
            
            # Write sequence sections
            if all_sequence_sections:
                f.write("\n")
                for line in all_sequence_sections:
                    f.write(line + "\n")
        
        print(f"Combined {len(batch_output_files)} batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")
        
        return True
        
    except Exception as e:
        print(f"Error combining glycosylation outputs: {e}")
        return False

def _combine_phosphorylation_outputs(batch_output_files, final_output_file):
    """Combine phosphorylation prediction batch outputs (NetPhos format)"""
    try:
        total_sequences = 0
        all_prediction_lines = []
        
        for i, batch_file in enumerate(batch_output_files):
            try:
                with open(batch_file, 'r') as f:
                    content = f.read()
                
                # Extract sequences count from header  
                lines = content.split('\n')
                for line in lines:
                    if 'amino acids' in line and line.startswith('>'):
                        try:
                            seq_count = int(line.split()[1])
                            total_sequences += seq_count
                        except:
                            total_sequences += 1  # Fallback
                        break
                
                # Collect prediction lines
                for line in lines:
                    if line.startswith('# ') and len(line.split()) >= 7:
                        # This is a prediction line
                        all_prediction_lines.append(line)
                        
            except Exception as e:
                print(f"Error reading phosphorylation batch file {batch_file}: {e}")
                continue
        
        # Write combined phosphorylation output
        with open(final_output_file, 'w') as f:
            # Write header
            gene_name = os.path.basename(final_output_file).replace('.out', '').replace('-netphos', '')
            f.write(f">{gene_name}\t{total_sequences} amino acids\n")
            f.write("#\n")
            f.write("#  prediction results\n")
            f.write("#\n")
            f.write("# Sequence\t\t   # x   Context     Score   Kinase    Answer\n")
            f.write("# " + "-" * 67 + "\n")
            
            # Write all predictions
            for line in all_prediction_lines:
                f.write(line + "\n")
        
        print(f"Combined {len(batch_output_files)} phosphorylation batch files into {final_output_file}")
        print(f"Total sequences: {total_sequences}")
        
        return True
        
    except Exception as e:
        print(f"Error combining phosphorylation outputs: {e}")
        return False

def run_docker_command(docker_image, fasta_file, command_template, output_file, timeout=300):
    """Generic Docker execution wrapper
    
    Args:
        docker_image: Docker image name
        fasta_file: Input FASTA file path
        command_template: Command template (e.g., "./ape -m netphos {input}")
        output_file: Output file path
        timeout: Command timeout in seconds
        
    Returns:
        tuple: (success, output_content, error_message)
    """
    work_dir = tempfile.mkdtemp()
    
    try:
        # Copy the FASTA file to the work directory
        docker_input = os.path.join(work_dir, "input.fasta")
        shutil.copy2(fasta_file, docker_input)
        
        # Build Docker command
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{work_dir}:/data",
            docker_image
        ]
        
        # Add command from template
        command = command_template.format(input="/data/input.fasta")
        docker_cmd.extend(command.split())
        
        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        
        if result.returncode == 0:
            # Save output to file
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            return True, result.stdout, None
        else:
            return False, result.stdout, result.stderr
            
    except subprocess.TimeoutExpired:
        return False, "", f"Docker command timed out after {timeout} seconds"
    except Exception as e:
        return False, "", str(e)
    finally:
        # Clean up work directory
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir, ignore_errors=True)

def extract_mutation_from_sequence_name(seq_name):
    """Extract mutation ID from sequence name (e.g., 'ZFP36-C330T' → 'C330T')
    
    Args:
        seq_name: Sequence name in format 'GENE-MUTATION' or just 'GENE'
        
    Returns:
        tuple: (gene, mutation_id) or (gene, None) if no mutation found
    """
    if '-' in seq_name:
        parts = seq_name.rsplit('-', 1)
        if len(parts) == 2:
            return parts[0], parts[1]  # gene, mutation_id
    
    return seq_name, None

def get_mutation_data_bioAccurate_unified(aaposaa):
    """Extract amino acid position and amino acids from mutation notation (e.g., 'K541E' -> 541, ('K', 'E'))
    
    This is a unified version that works for both NetPhos and NetNGlyc pipelines.
    
    Args:
        aaposaa: Amino acid mutation string (e.g., 'K541E', 'Y110F')
        
    Returns:
        tuple: (position, (original_aa, mutant_aa)) or (None, None) if invalid
    """
    import re
    
    # Skip stop codons
    if 'Stop' in aaposaa or 'Sto' in aaposaa:
        return None, None
    
    # Match pattern: [Letter][Numbers][Letter]
    match = re.match(r'[A-Z](\d+)[A-Z]', aaposaa)
    if match:
        position = int(match.group(1))
        original_aa = aaposaa[0]
        mutant_aa = aaposaa[-1]
        return position, (original_aa, mutant_aa)
    return None, None

def process_single_mutation_for_sequence(seq_name, predictions, mapping_dict, is_mutant=True, tool_type='netphos'):
    """Process predictions for one sequence against its specific mutation
    
    Args:
        seq_name: Sequence name (e.g., 'ZFP36-C330T')
        predictions: List of predictions for this sequence
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant sequences (should be True for single-mutation processing)
        tool_type: 'netphos' or 'netnglyc' for tool-specific field handling
        
    Returns:
        list: Filtered predictions with pkeys for the specific mutation
    """
    if not is_mutant:
        raise ValueError("process_single_mutation_for_sequence should only be used for mutant sequences")
    
    # Extract mutation ID from sequence name
    gene, mutation_id = extract_mutation_from_sequence_name(seq_name)
    if mutation_id is None:
        return []
    
    # Look up this mutation in the mapping
    if mutation_id not in mapping_dict:
        return []
    
    aaposaa = mapping_dict[mutation_id]  # e.g., "Y110F"
    
    # Parse amino acid position and mutation info
    position_data = get_mutation_data_bioAccurate_unified(aaposaa)
    if position_data[0] is None:
        return []
    
    aa_pos = position_data[0]  # e.g., 110
    aa_tuple = position_data[1]  # e.g., ('Y', 'F')
    target_aa = aa_tuple[1]  # F for mutant
    
    # Filter predictions for this specific mutation
    results = []
    for pred in predictions:
        # Get position from prediction
        if tool_type == 'netphos':
            pred_pos = pred['pos']
            pred_aa = pred['amino_acid']
        elif tool_type == 'netnglyc':
            pred_pos = pred['position']
            pred_aa = pred['sequon'][0] if pred['sequon'] else None
        else:
            raise ValueError(f"Unsupported tool_type: {tool_type}")
        
        # Check if this prediction matches the mutation position and amino acid
        if pred_pos == aa_pos and pred_aa == target_aa:
            # Create pkey for this match
            pkey = f"{gene}-{mutation_id}"
            
            # Add pkey and fix Gene field to prediction
            result_pred = pred.copy()
            result_pred['pkey'] = pkey
            # Fix Gene field to just gene name (not gene-mutation)
            if tool_type == 'netphos':
                result_pred['Gene'] = gene
            elif tool_type == 'netnglyc':
                result_pred['Gene'] = gene
            results.append(result_pred)
    
    return results

def parse_predictions_with_mutation_filtering(predictions, mapping_dict, is_mutant, threshold=0.0, yes_only=False, tool_type='netphos'):
    """Universal prediction filtering logic for both NetPhos and NetNGlyc
    
    Args:
        predictions: List of all predictions
        mapping_dict: Dictionary mapping mutation_id -> aaposaa
        is_mutant: Whether processing mutant (True) or wildtype (False) sequences
        threshold: Score threshold for filtering
        yes_only: Only include predictions with 'YES' answer (NetPhos only)
        tool_type: 'netphos' or 'netnglyc' for tool-specific handling
        
    Returns:
        list: Filtered predictions with pkeys
    """
    results = []
    
    if is_mutant:
        # Group predictions by sequence name for single-mutation processing
        seq_predictions = {}
        for pred in predictions:
            seq_name = pred.get('Gene', '') if tool_type == 'netphos' else pred.get('seq_name', '')
            if seq_name not in seq_predictions:
                seq_predictions[seq_name] = []
            seq_predictions[seq_name].append(pred)
        
        # Process each sequence separately with its specific mutation
        for seq_name, seq_preds in seq_predictions.items():
            seq_results = process_single_mutation_for_sequence(
                seq_name, seq_preds, mapping_dict, is_mutant=True, tool_type=tool_type
            )
            
            # Apply additional filters
            for result in seq_results:
                # Apply threshold filter
                score_field = 'score' if tool_type == 'netphos' else 'potential'
                if result[score_field] < threshold:
                    continue
                    
                # Apply yes_only filter (NetPhos only)
                if yes_only and tool_type == 'netphos' and result.get('answer') != 'YES':
                    continue
                
                results.append(result)
    
    else:
        # Wildtype processing - use existing bulk logic (not implemented here)
        # This should use the existing wildtype processing logic from each pipeline
        raise NotImplementedError("Wildtype processing should use existing pipeline-specific logic")
    
    return results
