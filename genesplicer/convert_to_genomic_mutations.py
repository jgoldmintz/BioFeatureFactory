import os
import sys
from pathlib import Path
import argparse
import subprocess

script_dir = Path(__file__).parent.absolute()
dependencies_dir = script_dir.parents[2] / 'dependencies'
sys.path.insert(0, str(dependencies_dir))

from utility import convert_position, trim_muts, get_mutation_data, read_fasta


def align_seqs(file):
    seq = read_fasta(file)
    with open('tmp.tmp', 'w') as tmp:
        tmp.write('>ORF\n'+seq['ORF']+'\n>genomic\n'+seq['genomic'])
    cmd = 'mafft tmp.tmp > tmpMSA.fasta && rm tmp.tmp'
    subprocess.getoutput(cmd)
    seq = read_fasta('tmpMSA.fasta')
    subprocess.getoutput('rm tmpMSA.fasta')
    return seq

def get_genomic_ntposnt(indir, outdir, mutdir):
    for file in os.scandir(indir):
        gene = file.name.split(".")[0]
        #print('Retrieving genomic sequences for {}'.format(gene))
        print('Aligning sequences')
        allseq = align_seqs(file)
        orf_aligned = allseq['ORF']
        genomic_aligned = allseq['genomic']
        print('Converting genomic positions')
        mut_list = trim_muts(mutdir+'/'+gene+'_mutations.csv')

        with open(outdir+'/'+gene+'_genomic_mutations.csv','w') as f:
            f.write('mutant\n')
            for ntposnt in mut_list:
                pos, mut =get_mutation_data(ntposnt)
                genomicPos, error  = convert_position(orf_aligned,genomic_aligned,pos)
                f.write(mut[0]+str(genomicPos)+mut[1]+'\n')

        print('All positions converted for {}'.format(gene))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert ORF mutations to genomic.')
    parser.add_argument('indir', help='Input directory of FASTA files assumed -- used to retrieve genenames')
    parser.add_argument('outdir', help='Output directory of new genomic mutations -- assumes that the mutations'
                                       'in this directory are in csv files where the first row is "mutant" followed by '
                                       'ntposnt')
    parser.add_argument('mutdir', help='Mutation directory of ORF mutations')

    args = parser.parse_args()
    get_genomic_ntposnt(args.indir, args.outdir, args.mutdir)