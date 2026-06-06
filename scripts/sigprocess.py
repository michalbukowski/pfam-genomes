#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Parses raw SignalP output, extracts the cleavage site position and signal
# peptide probability for each sequence, deduplicates per sequence ID keeping
# the highest-probability prediction, and writes the results to a TSV file.
# Arguments:
# --pthresh : signal peptide probability threshold (parsed but reserved for
#             future filtering; currently all detected signal peptides are kept)
# --sigres  : raw SignalP output file
# --output  : processed SignalP results in TSV format
# USAGE:
# ./sigprocess.py --pthresh PROB_THRESH --sigres SIGNALP_RES --output PROCESSED_RES

#-------------------------------------------------------------------------------
import argparse
import pandas as pd

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --pthresh : signal peptide probability threshold
       --sigres  : raw SignalP output file
       --output  : processed SignalP results in TSV format
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--pthresh', type=float, required=True,
        help='Probability threshold for signal peptide detection')
    parser.add_argument('--sigres', type=str, required=True,
        help='Final SignalP results for Gram+ and Gram- bacteria')
    parser.add_argument('--output', type=str, required=True,
        help='Processed SignalP results in TSV format')

    args = parser.parse_args()
    return args

def main():
    '''Parses raw SignalP output, extracts cleavage site position and probability,
       deduplicates per sequence ID keeping the highest-probability hit, and
       writes the results to a TSV file.
    '''
    
    # Parse the command-line args
    args = parse_args()
    
    # Open SingalP results and parse them into a DataFrame. Drop duplicates in
    # regard to analysed sequence id (ID), leve those with highest propabilities
    # of possesing an N-terminal. Save the filtered DataFrame.
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
    sig_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

