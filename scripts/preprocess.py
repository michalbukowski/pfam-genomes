#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Perform pre-processing and E-value filtering on raw HMMsearch results,
# and labels hits with domain group names provided in input data. Saves only
# selected columns in a TSV file. Arguments:
# --qcovt   : query domain covarage by target sequence threshold
# --iEvalue : independent E-value threshold for each domain hit
# --cols    : names for columns reported by hmmsearch
# --leave   : names of columns that are to be preserved
# --domdata : a TSV file with data describing searched domains, it must
#             link Pfam accession versions (pfam_acc column) with names of
#             groups those are assigned to (group column)
# --domtbl  : per domain (domtblout) hmmsearch tabular output data file
# --output  : final TSV output file with filtered data
#USAGE:
# ./preprocess.py --qcovt   QCOVT_THRESH    --iEvalue IEVALUE_THRESH \
#                 --cols    HMMOUT_COLNAMES --leave   LEAVE_COLNAMES \
#                 --domdata INPUT_DOMDESC   --domtbl  DOM_HMMOUT     \
#                 --output  PREPROC_HMMOUT

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --qcovt   : query domain covarage by target sequence threshold
       --iEvalue : independent E-value threshold for each domain hit
       --cols    : names for columns reported by hmmsearch
       --leave   : names of columns that are to be preserved
       --domdata : a TSV file with data describing searched domains, it must
                   link Pfam accession versions (pfam_acc column) with names of
                   groups those are assigned to (group column)
       --domtbl  : per domain (domtblout) hmmsearch tabular output data file
       --output  : final TSV output file with filtered data
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--qcovt', type=float, required=True,
        help='Query domain covarage by target sequence threshold')
    parser.add_argument('--iEvalue', type=float, required=True,
        help='Independent E-value threshold for each domain hit')
    parser.add_argument('--cols', type=str, required=True,
        help='Names of columns reported by hmmsearch')
    parser.add_argument('--leave', type=str, required=True,
        help='Names of columns that are to be preserved')
    parser.add_argument('--domdata', type=str, required=True,
        help='TSV file with data on searched domains it must link Pfam '  +
             'accession versions (pfam_acc column) with names of groups ' +
             'those are assigned to (group column)')
    parser.add_argument('--domtbl', type=str, required=True,
        help='Per domain (domtblout) hmmsearch tabular output data')
    parser.add_argument('--output', type=str, required=True,
        help='TSV output file with preprocessed data')

    args = parser.parse_args()
    return args

def main():
    '''The entry point function that perform pre-processing and E-value filtering
       on raw HMMsearch results, and labels hits with domain group names provided
       in input data. Saves only selected columns in a TSV file.
    '''
    # Parse args.
    args = parse_args()
    # Try to read raw HMMsearch per domain output. Split the content in respect
    # to white-space substrings and give columns the requested names (args.cols).
    # If file is empty (no relevant hits found), save and emtpy DataFrame with
    # requested output columns (args.leave).
    try:
        hmm_df = pd.read_csv(args.domtbl, header=None, comment='#', sep='\s+')
    except pd.errors.EmptyDataError:
        hmm_df = pd.DataFrame(columns=args.leave.split())
        hmm_df.to_csv(args.output, index=False, sep='\t')
        sys.exit()
    # Colums were created using '\s+' separator and that means the last column
    # with metadata on sequecne, which values contains spaces, is split too.
    # Take 21 first columns and combine the rest into one. Then extract sequence
    # metadata to new, propely named, columns.
    hmm_df = pd.concat([hmm_df.loc[:, :21], hmm_df.loc[:, 22:].apply(
        lambda row: ' '.join(row), axis=1)], axis=1)
    hmm_df.columns = args.cols.split(' ')
    for col in 'srcid start end asmacc clustid'.split():
        hmm_df[col] = hmm_df['desc'].str.extract(f'{col}=([^\ {eol}]+)')
    # Calculate query domain coverage by target protein sequence (qcovt) and
    # filter against it and independent E-value.
    hmm_df['qcovt'] = ((hmm_df['hmm_to'] - hmm_df['hmm_from']).abs() + 1) \
                      / hmm_df['qlen']
    hmm_df = hmm_df[
        (hmm_df['qcovt'] >= args.qcovt) &
        (hmm_df['i-Evalue'] <= args.iEvalue)
    ]
    # Load to a DataFrame searched domain data, merge with HMMsearch data.
    dom_df = pd.read_csv(args.domdata, sep='\t')
    hmm_df = hmm_df.merge(dom_df[['pfm_name', 'group']], left_on='qname',
                          right_on='pfm_name', how='left')
    hmm_df['group'].fillna('Unassigned', inplace=True)
    # Select the requested columns (args.leave).
    hmm_df = hmm_df[args.leave.split()]
    # Sort hmm_df rows by assembly accession, target protein sequecne id/name and
    # the star position of matched domain. It is necessary not to assign to
    # a protein wrong architecture (domain order) by other scripts.
    hmm_df.sort_values('asmacc tname env_from'.split(), inplace=True)
    #Save the preprocessed data.
    hmm_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

