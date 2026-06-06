#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Filters integrated HMMsearch and SignalP results by domain architecture.
# For each protein, aggregates its domain hits into a dash-joined architecture
# string ordered by position, then retains only those proteins whose architecture
# matches the provided regex. Outputs a list of matching protein IDs. Arguments:
# --arch    : architecture regex prefixed with 'GRP|', 'NAM|' or 'ACC|' to
#             indicate whether domains are identified by group names, Pfam domain
#             names or Pfam accession numbers, respectively
# --allres  : integrated HMMsearch and SignalP results TSV
# --output  : single-column TSV with protein IDs matching the architecture
#USAGE:
# ./archfilter.py --arch DOMARCH_REGEX --allres INTEGRATED_RES --output TNAMES_TSV

#-------------------------------------------------------------------------------
import argparse
import pandas as pd
from os import linesep as eol

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --arch    : expected domain arrangements as a set of regex, each of which
                   starts with 'GRP|', 'NAM|' or 'ACC|' depending whether domain
                   arrangement describes domains by their group names, Pfam
                   domain names or accession numbers
       --allres  : final filtered HMMsearch and SignalP results
       --output  : final TSV output file with filtered data
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--arch', type=str, required=True,
        help='expected domain arrangements as a set of regex, each of which ' + \
             'starts with \'GRP|\', \'NAM|\' or \'ACC|\' depending whether '  + \
             'domain arrangement describes domains by their group names, '    + \
             'Pfam domain names or accession numbers')
    parser.add_argument('--allres', type=str, required=True,
        help='Final filtered HMMsearch and SignalP results')
    parser.add_argument('--output', type=str, required=True,
        help='TSV output file with filtered data')

    args = parser.parse_args()
    return args

def main():
    '''Reads integrated results, aggregates domain hits per protein into
       a dash-joined architecture string ordered by position, filters proteins
       whose architecture matches the --arch regex, and writes the matching
       protein IDs to the output TSV.
    '''
    # Parse the command-line args and make sure that the architecture regex
    # contains a proper prefix.
    args = parse_args()
    if not any(args.arch.startswith(prefix) for prefix in ['GRP|', 'NAM|', 'ACC|']):
        raise Exception('The value of --arch must start either with \'GRP|\', ' +
                        f'\'NAM|\' or \'ACC|\': {args.arch}')
    # Read preprocessed HMMserch results.
    all_df = pd.read_csv(args.allres, sep='\t')
    print(f'Found {all_df.shape[0]:>3} results, ', end='')
    # Use either group column, qname or qacc (Pfam accession version) depending
    # on input regex prefix. Remove the prefix from regex.
    if args.arch.startswith('GRP|'):
        col = 'group'
    elif args.arch.startswith('NAM|'):
        col = 'qname'
    else:
        col = 'qacc'
    args.arch = args.arch[4:]

    # Sort the final DataFrame with respect to domains locations (order domains).
    all_df.loc[all_df['group'] == 'SP', 'qacc'] = 'SP'
    all_df.sort_values('tname env_from env_to'.split(), inplace=True)
    
    # Having the DataFrame sorted, group rows in respect to protein sequecne id/name,
    # join values of remaining colums (domain group and Pfam accession version)
    # with single dash. That creates for every target protein sequence a string
    # describing its domain architecture in regard to domain groups as well as
    # domain Pfam accession versions.
    agg_df = all_df['tname group qacc qname'.split()].groupby('tname').agg(
        lambda values: '-'.join(values))
    print(agg_df)

    # Filter the aggregated rows in respect to either group or qacc column
    # (whichever is selected in command line arguments) using architecure
    # regex. Use agg_df index to select relevant rows from the original hmm_df.
    agg_df = agg_df.loc[ agg_df[col].str.fullmatch(args.arch) ]
    all_df = all_df[ all_df['tname'].isin(agg_df.index) ]
    print(f'of which {all_df.shape[0]:>3} remained')
    # Save the filtered data.
    all_df['tname'].to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

