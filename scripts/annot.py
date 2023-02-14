#!/usr/bin/env python3

#-------------------------------------------------------------------------------
import argparse
import pandas as pd
from os import linesep as eol
from lib import fasta

#-------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--sigres', type=str, required=True,
        help='final signalp results for Gram+ and Gram- bacteria')
    parser.add_argument('--hmmres', type=str, required=True,
        help='final hmmsearch results')
    parser.add_argument('--seqs', type=str, required=True,
        help='fasta file with seqs')
    parser.add_argument('--output', type=str, required=True,
        help='GFF3 file with annotations')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    seqs = fasta.read_fasta(args.seqs)
    
    hmm_df = pd.read_csv(args.hmmres, sep='\t')
    hmm_df = hmm_df[ hmm_df['tname'].isin(seqs) ]
    
    with open(args.sigres) as f:
        f.readline()
        names = f.readline()[2:-1].split('\t')
    sig_df = pd.read_csv(args.sigres, names=names, comment='#', sep='\t')
    sig_df.dropna(inplace=True)
    if sig_df.shape[0] > 0:
        cols = sig_df['CS Position'].str.split('(?:\-)|(?:\ )', expand=True)
        sig_df['pos'] = cols[2].astype(int)
        sig_df['prob'] = cols[7].astype(float)
        sig_df.sort_values(['ID', 'prob'], ascending=False, inplace=True)
        sig_df.drop_duplicates('ID', inplace=True)
        sig_df.set_index('ID', drop=True, inplace=True)
    
    fout = open(args.output, 'w')
    fout.write(f'##gff-version 3{eol}')
    for tname, sub_df in hmm_df.groupby('tname'):
        fout.write(f'##sequence-region {tname} 1 {len(seqs[tname])}{eol}')
        if tname in sig_df.index:
            srcid = sub_df['srcid'].iloc[0]
            pos   = sig_df.loc[tname, 'pos']
            prob  = (sig_df.loc[tname, 'prob']*100).round(0).astype(int)
            name  = f'Name=SP ({prob:d}%)'
            line  = tname, srcid, 'SP', '1', str(pos), '.', '+', \
                    '.', name
            fout.write('\t'.join(line) + eol)
        for _, (qname, qacc, group, srcid, env_from, env_to) in \
            sub_df['qname qacc group srcid env_from env_to'.split()].iterrows():
            name = f'Name={group} ({qname} [{qacc}])'
            line = tname, srcid, group, str(env_from), str(env_to), '.', '+', \
                   '.', name
            fout.write('\t'.join(line) + eol)
    fout.close()

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
