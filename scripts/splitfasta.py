#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Splits a protein FASTA file into a given number of roughly equal-sized chunks
# and writes each chunk as a separate FASTA file (chunk_NNN.faa) to the output
# directory. Arguments:
# --chunks : number of chunks to split the input sequences into
# --input  : input protein FASTA file
# --outdir : output directory for the chunk FASTA files
# USAGE:
# ./splitfasta.py --chunks N_CHUNKS --input INPUT_FASTA --outdir OUTPUT_DIR

#-------------------------------------------------------------------------------
import argparse
from lib.fasta import read_fasta

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --chunks : number of chunks to split the input sequences into
       --input  : input protein FASTA file
       --outdir : output directory for the chunk FASTA files
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--chunks', type=int, required=True,
        help='')
    parser.add_argument('--input', type=str, required=True,
        help='')
    parser.add_argument('--outdir', type=str, required=True,
        help='')

    args = parser.parse_args()
    return args

def main():
    '''Reads all sequences from the input FASTA file, divides them into the
       requested number of roughly equal-sized chunks, and writes each chunk
       to a separate file (chunk_NNN.faa) in the output directory.
    '''
    
    args = parse_args()
    seqs = read_fasta(args.input)
    seqs = list( seqs.values() )
    chunk_len = int( round( len(seqs) / args.chunks ) )
    for i in range(0, len(seqs), chunk_len):
        chunk = seqs[i:i+chunk_len]
        buf = ''.join([ seq.fasta() for seq in chunk ])
        with open(f'{args.outdir}/chunk_{i:03d}.faa', 'w') as f:
            f.write(buf)

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

