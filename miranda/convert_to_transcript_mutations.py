import os
import sys
import argparse
import subprocess
from pathlib import Path

script_dir = Path(__file__).parent.absolute()
dependencies_dir = script_dir.parents[2] / 'dependencies'
sys.path.insert(0, str(dependencies_dir))

from utility import trim_muts, get_mutation_data, read_fasta, convert_position

def align_seqs(file):
    seq = read_fasta(file)
    with open('tmp.tmp', 'w') as tmp:
        tmp.write('>ORF\n'+seq['ORF']+'\n>transcript\n'+seq['transcript'])
    cmd = 'mafft tmp.tmp > tmpMSA.fasta'
    subprocess.getoutput(cmd)
    seq = read_fasta('tmpMSA.fasta')
    #subprocess.getoutput('rm tmpMSA.fasta')
    return seq

def get_transcript_ntposnt(indir, outdir, mutdir):
    """
    Convert ORF mutation positions to transcript positions using pre-aligned sequences.

    Args:
        indir: Directory containing FASTA files with aligned sequences
        outdir: Directory to write transcript mutation CSV files
        mutdir: Directory containing ORF mutation CSV files
    """

    # Ensure output directory exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for file in os.scandir(indir):
        if not file.name.endswith('.fasta'):
            continue

        gene = file.name.split(".")[0]
        print(f'Processing {gene}...')

        print('     -- aligning sequences')
        allseq = align_seqs(file.path)

        orf_aligned = allseq['ORF']
        transcript_aligned = allseq['transcript']

        print(f'  - ORF aligned length: {len(orf_aligned)}')
        print(f'  - Transcript aligned length: {len(transcript_aligned)}')

        # Get mutation list
        mut_file = os.path.join(mutdir, f'{gene}_mutations.csv')
        if not os.path.exists(mut_file):
            print(f'Warning: No mutation file found for {gene}, skipping...')
            continue

        print(f'Converting ORF positions to transcript positions for {gene}...')
        mut_list = trim_muts(mut_file)

        converted_count = 0
        error_count = 0

        # Write converted mutations
        with open(os.path.join(outdir, f'{gene}_transcript_mutations.csv'), 'w') as f:
            f.write('mutant\n')
            for ntposnt in mut_list:
                pos, mut = get_mutation_data(ntposnt)
                transcriptPos, error = convert_position(orf_aligned, transcript_aligned, pos)

                if error:
                    print(f'    - Error converting {ntposnt}: {error}')
                    error_count += 1
                    subprocess.getoutput(f'mv tmpMSA.fasta error_{gene}.fasta')
                    # You might want to skip writing this mutation or write a placeholder
                    continue

                f.write(f'{mut[0]}{transcriptPos}{mut[1]}\n')
                converted_count += 1

        print(f'  - Successfully converted {converted_count} positions for {gene}')
        if error_count > 0:
            print(f'  - {error_count} positions could not be converted (gaps in transcript)')


def main():
    parser = argparse.ArgumentParser(
        description='Convert ORF mutations to transcript positions using pre-aligned sequences.'
    )
    parser.add_argument(
        'indir',
        help='Input directory of FASTA files containing aligned sequences'
    )
    parser.add_argument(
        'outdir',
        help='Output directory for transcript mutation CSV files'
    )
    parser.add_argument(
        'mutdir',
        help='Directory containing ORF mutation CSV files (format: GENE_mutations.csv)'
    )

    args = parser.parse_args()

    # Validate input directories
    if not os.path.exists(args.indir):
        print(f"Error: Input directory '{args.indir}' does not exist")
        return

    if not os.path.exists(args.mutdir):
        print(f"Error: Mutation directory '{args.mutdir}' does not exist")
        return

    get_transcript_ntposnt(args.indir, args.outdir, args.mutdir)
    print("\nConversion complete!")


if __name__ == '__main__':
    main()