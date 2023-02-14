#!/usr/bin/env python3

#-------------------------------------------------------------------------------
import argparse, sys
from os import linesep as eol
from itertools import count
from lib.fasta import fasta, fasta_meta

#-------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', type=str,
        help='fasta file with sequences for clustering')
    parser.add_argument('--output', type=str, required=True,
        help='fasta file with representative sequences')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    ngen = count()
    seqs = {}
    meta_keys = 'seqid asmacc srcid start end'.split(' ')
    
    fin    = sys.stdin if args.input is None else open(args.input)
    frepr  = open(args.output, 'w')
    fclust = open(args.output + '.tsv', 'w')
    fclust.write('clust\t' + '\t'.join(meta_keys) + eol)
    
    header, seq = None, ''
    for line in fin:
        if line[0] == '>':
            line = line[1:].rstrip(eol)
            if header is not None and seq != '' and 'X' not in seq:
                key = seq[1:]
                if key not in seqs:
                    seqs[key] = next(ngen)
                    frepr.write(fasta(header, '', seq))
                clust = seqs[key]
                meta = fasta_meta(header, meta_keys)
                meta = '\t'.join( meta[key] for key in meta_keys )
                fclust.write(f'{clust}\t{meta}{eol}')
            header = line
            seq = ''
        elif header is not None:
            seq += line.rstrip(eol)
    
    fclust.close()
    frepr.close()
    if fin is not sys.stdin:
       fin.close()

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

