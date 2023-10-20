#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Merges preprocessed HMMsearch results in respect to target query id,
# that gives one row - one protein with an architecture assigned, i.e. sequence
# of domains. Next filters the results in respect to the requested domain
# architecture. Arguments:
# --arch    : expected domain arrangement as regex, starts with
#             'GRP|' or 'ACC|' depending whether arrangement describes
#             domains by their Pfam accession versions or their group names
# --hmmres  : preprocessed HMMsearch results
# --output  : final TSV output file with filtered data
#USAGE:
# ./filter.py --arch    DOMARCH_REGEX --hmmres INPUT_HMMRES    \
#             --output  FINAL_HMMOUT

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --arch    : expected domain arrangement as regex, starts with
                    'GRP|' or 'ACC|' depending whether arrangement describes
                    domains by their Pfam accession versions or their group names
       --hmmres  : preprocessed HMMsearch results
       --output  : final TSV output file with filtered data
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--arch', type=str, required=True,
        help='Expected domain arrangement as regex, starts with \'GRP|\', ' +
             '\'NAM|\' or \'ACC\' depending whether arrangement describes ' +
             'domains by their Pfam accession versions, Pfam domain names or' +
             'their group names')
    parser.add_argument('--hmmres', type=str, required=True,
        help='Preprocessed HMMsearch results')
    parser.add_argument('--output', type=str, required=True,
        help='TSV output file with filtered data')

    args = parser.parse_args()
    return args

def main():
    '''The entry point function that merges preprocessed HMMsearch results
       in respect to target query id, that gives one row - one protein with
       an architecture assigned, i.e. sequence of domains. Next filters
       the results in respect to the requested domain architecture.
    '''
    # Parse args and make sure that architecture regex contains a proper prefix.
    args = parse_args()
    if not any(args.arch.startswith(prefix) for prefix in ['GRP|', 'NAM|', 'ACC|']):
        raise Exception('The value of --arch must start either with \'GRP|\', ' +
                        '\'NAM|\' or \'ACC\'')
    # Read preprocessed HMMserch results.
    hmm_df = pd.read_csv(args.hmmres, sep='\t')
    print(f'Found {hmm_df.shape[0]:>3} results, ', end='')
    # Use either group column, qname or qacc (Pfam accession version) depending
    # on input regex prefix. Remove the prefix from regex.
    if args.arch.startswith('GRP|'):
        col = 'group'
    elif args.arch.startswith('NAM|'):
        col = 'qname'
    else:
        col = 'qacc'
    args.arch = args.arch[4:]
    # In the preprocessed data rows have been already sorted by assembly accession,
    # target protein sequecne id/name and the star position of matched domain to
    # avoid assigning to a protein wrong architecture (domain order).
    # Having that done, group rows in respect to protein sequecne id/name,
    # join values of remaining colums (domain group and Pfam accession version)
    # with single dash. That creates for every target protein sequence a string
    # describing its domain architecture in regard to domain groups as well as
    # domain Pfam accession versions.
    agg_df = hmm_df['tname group qacc qname'.split()].groupby('tname').agg(
        lambda values: '-'.join(values))
    # Filter the aggregated rows in respect to either group or qacc column
    # (whichever is selected in command line arguments) using architecure
    # regex. Use agg_df index to select relevant rows from the original hmm_df.
    agg_df = agg_df.loc[ agg_df[col].str.fullmatch(args.arch) ]
    hmm_df = hmm_df[ hmm_df['tname'].isin(agg_df.index) ]
    print(f'of which {hmm_df.shape[0]:>3} remained')
    # Save the filtered data.
    hmm_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

