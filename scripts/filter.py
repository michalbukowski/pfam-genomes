#!/usr/bin/env python3

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--qcovt', type=float, required=True,
        help='query (profile) covarage by target sequence')
    parser.add_argument('--iEvalue', type=float, required=True,
        help='independent E-value threshold for each domain hit')
    parser.add_argument('--cols', type=str, required=True,
        help='names of columns reported by hmmsearch')
    parser.add_argument('--leave', type=str, required=True,
        help='names of columns that are to be preserved')
    parser.add_argument('--arch', type=str, required=True,
        help='expected domain arrangement as regex, starts with \'GRP|\' or ' +
             '\'ACC|\' depending whether arrangement describes domains by '   +
             'their PFAM accession numbers or their group names')
    parser.add_argument('--domdata', type=str, required=True,
        help='TSV file with data on searched domains')
    parser.add_argument('--domtbl', type=str, required=True,
        help='per domain (domtblout) hmmsearch tabular output data')
    parser.add_argument('--output', type=str, required=True,
        help='TSV output file with filtered data')

    args = parser.parse_args()
    return args

def basic_filter(args, hmm_df):
    hmm_df = pd.concat([hmm_df.loc[:, :21], hmm_df.loc[:, 22:].apply(
        lambda row: ' '.join(row), axis=1)], axis=1)
    hmm_df.columns = args.cols.split(' ')
    for col in 'srcid start end asmacc'.split():
        hmm_df[col] = hmm_df['desc'].str.extract(f'{col}=([^\ {eol}]+)')
    hmm_df['qcovt'] = ((hmm_df['hmm_to'] - hmm_df['hmm_from']).abs() + 1) \
                      / hmm_df['qlen']
    hmm_df = hmm_df[
        (hmm_df['qcovt'] >= args.qcovt) &
        (hmm_df['i-Evalue'] <= args.iEvalue)
    ]
    return hmm_df

def arch_filter(args, hmm_df, dom_df):
    col = 'group' if args.arch.startswith('GRP|') else 'qacc'
    args.arch = args.arch[4:]
    hmm_df = hmm_df.merge(dom_df[['pfm_acc', 'group']], left_on='qacc',
                          right_on='pfm_acc', how='left')
    hmm_df.sort_values(['asmacc', 'tname', 'env_from'], inplace=True)
    agg_df = hmm_df[['tname', 'group', 'qacc']].groupby('tname').agg(
        lambda values: '-'.join(values))
    agg_df = agg_df.loc[ agg_df[col].str.fullmatch(args.arch) ]
    hmm_df = hmm_df[ hmm_df['tname'].isin(agg_df.index) ]
    return hmm_df

def main():
    args = parse_args()
    if not any(args.arch.startswith(prefix) for prefix in ['GRP|', 'ACC|']):
        raise Exception('The value of --arch must start either with \'GRP|\' ' +
                        'or \'ACC|\'')
    try:
        hmm_df = pd.read_csv(args.domtbl, header=None, comment='#', sep='\s+')
    except pd.errors.EmptyDataError:
        hmm_df = pd.DataFrame(columns=args.leave.split())
        hmm_df.to_csv(args.output, index=False, sep='\t')
        sys.exit()
    dom_df = pd.read_csv(args.domdata, sep='\t')
    print(f'Found {hmm_df.shape[0]:>3} results, ', end='')
    hmm_df = basic_filter(args, hmm_df)
    hmm_df = arch_filter(args, hmm_df, dom_df)
    print(f'of which {hmm_df.shape[0]:>3} remained')
    hmm_df = hmm_df[args.leave.split()]
    hmm_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
