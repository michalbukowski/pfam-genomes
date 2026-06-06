#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Combines HMMsearch results and SignalP results into a single unified table.
# For each protein that has both HMM domain hits and a SignalP hit, the signal
# peptide is appended as an extra row with qname='SP', qacc=confidence%, and
# group='SP'. Arguments:
# --sigres  : preprocessed SignalP results TSV
# --hmmres  : filtered and cleaned HMMsearch results TSV
# --output  : combined TSV with domain hits and signal peptide rows
# USAGE:
# ./join.py --sigres SIGNALP_RES --hmmres HMMSEARCH_RES --output JOINED_RES

#-------------------------------------------------------------------------------
import argparse
from sys import exit
import pandas as pd
#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --sigres  : preprocessed SignalP results TSV
       --hmmres  : filtered and cleaned HMMsearch results TSV
       --output  : combined TSV with domain hits and signal peptide rows
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
    '''Loads and deduplicates HMMsearch and SignalP results, appends a signal
       peptide row for each protein present in both, and writes the combined
       table to the output TSV file.
    '''
    args = parse_args()
    
    # Load the final preprocessed HMMsearch results, and final SignalP results.
    hmm_df    = pd.read_csv(args.hmmres, sep='\t')
    sig_df    = pd.read_csv(args.sigres, sep='\t')
    
    # Remove duplicated rows since the same result might appear more than
    # once anytime a chosen architecture is a subset of another one.
    # TODO: is it necessary???
    hmm_df.drop_duplicates('tname qacc env_from env_to'.split(), inplace=True)
    sig_df.drop_duplicates('ID', inplace=True)
    
    # Add data from sig_df to hmm_df.
    hmm_gb = hmm_df.groupby('tname')
    print('HMM results groupped')
    sig_gb = sig_df.groupby('ID')
    print('SignalP results groupped')
    sp_data = []
    for tname, subhmm_df in hmm_gb:
        if tname in sig_gb.groups:
            srcid, start, end, asmacc, clustid = subhmm_df[
                'srcid start end asmacc clustid'.split()
            ].iloc[0]
            subsig_df = sig_gb.get_group(tname)
            assert subsig_df.shape[0] == 1,  'Unexpected subsig_df length: ' + \
                                            f'{subsig_df.shape[0]}'
            subsig_df = subsig_df.squeeze()
            qname     = 'SP'
            qacc      = f'{subsig_df.at["prob"]*100:0.0f}%'
            group     = 'SP'
            env_from  = 1
            env_to    = subsig_df.at['pos']
            vals = tname, srcid, start, end, asmacc, qname, qacc, clustid, \
                   group, env_from, env_to
            sp_data.append(vals)

    cols = 'tname srcid start end asmacc qname qacc clustid group'.split() + \
           'env_from env_to'.split()
    sp_df = pd.DataFrame(sp_data, columns=cols)
    
    all_df = pd.concat([hmm_df, sp_df])
    all_df.to_csv(args.output, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

