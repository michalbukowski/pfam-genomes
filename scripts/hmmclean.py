#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Removes from HMMsearch results domain hits that fall entirely within the
# signal peptide cleavage region identified by SignalP (env_to <= cleavage
# position), preventing signal peptides from being double-counted as domains.
# Arguments:
# --sigres  : processed SignalP results TSV
# --hmmres  : all-domains HMMsearch results TSV
# --output  : HMMsearch results with signal-peptide-overlapping hits removed
# USAGE:
# ./hmmclean.py --sigres SIGNALP_RES --hmmres HMMSEARCH_RES --output CLEANED_RES

#-------------------------------------------------------------------------------
import argparse
from sys import exit
import pandas as pd
#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --sigres  : processed SignalP results TSV
       --hmmres  : all-domains HMMsearch results TSV
       --output  : HMMsearch results with signal-peptide-overlapping hits removed
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--sigres', type=str, required=True,
        help='Final SignalP results for Gram+ and Gram- bacteria')
    parser.add_argument('--hmmres', type=str, required=True,
        help='Final HMMsearch results')
    parser.add_argument('--output', type=str, required=True,
        help='HMMsearch results cleaned from SP duplicates')

    args = parser.parse_args()
    return args

def main():
    '''Loads HMMsearch and SignalP results and, for each protein with a SignalP
       hit, drops HMM domain hits whose coordinates fall entirely within the
       signal peptide cleavage region. Writes the cleaned results to the output.
    '''
    
    # Parse the command-line args
    args = parse_args()
    
    # Load HMMsearch results.
    hmm_df = pd.read_csv(args.hmmres, sep='\t')
    
    # Load SignalP results.
    sig_df = pd.read_csv(args.sigres, sep='\t')
    
    # If HMMsearch or SignalP results are empty, save the former as merged and exit.
    if hmm_df.shape[0] == 0 or sig_df.shape[0] == 0:
        hmm_df.to_csv(args.output, index=False, sep='\t')
        exit(0)
    
    # Remove duplicated SP results from HMMsearch results.
    for tname, subhmm_df in hmm_df.groupby('tname'):
        subhmm_df = subhmm_df['env_from env_to'.split()].sort_values('env_from')
        subsig_df = sig_df.loc[sig_df['ID'] == tname]
        if subsig_df.shape[0] > 0:
            assert subsig_df.shape[0] == 1, f'Unexpected subsig_df length: {subsig_df.shape[0]}'
            subsig_df = subsig_df.squeeze()
            pos = subsig_df.at['pos']
            for ind, (env_from, env_to) in subhmm_df.iterrows():
                if env_to <= pos:
                    hmm_df.drop(ind, inplace=True)
                else:
                    break
    hmm_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

