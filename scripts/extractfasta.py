#!/usr/bin/env python3

#-------------------------------------------------------------------------------
import argparse
import pandas as pd
from lib import fasta

#-------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--seqdir', type=str, required=True,
        help='directory with per assembly seqence data')
    parser.add_argument('--hmmres', type=str, required=True,
        help='final hmmsearch results')
    parser.add_argument('--output', type=str, required=True,
        help='output file for extracted sequences')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    hmm_df = pd.read_csv(args.hmmres, sep='\t')
    hmm_df.drop_duplicates('tname', inplace=True)
    fout = open(args.output, 'w')
    for asmacc, sub_df in hmm_df[['asmacc', 'tname']].groupby('asmacc'):
        seqs = fasta.read_fasta(f'{args.seqdir}/{asmacc}.faa', sub_df['tname'].to_numpy())
        for seqid in sub_df['tname'].sort_values():
            fout.write(seqs[seqid].fasta())
    fout.close()

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

