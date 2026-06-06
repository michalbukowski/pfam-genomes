#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Fetches HMM profiles by name from a Pfam database HMM file and writes the
# selected profiles to a new HMM file. Domain names to fetch are read from
# a TSV file listing domains of interest. Arguments:
# --domains : TSV file listing domains of interest
# --namecol : column in the domains TSV file that holds domain names
# --pfamdb  : Pfam database HMM file
# --output  : output HMM file with profiles for the selected domains
#USAGE:
# ./hmmfetch.py --domains DOMAINS_TSV  --namecol DOM_NAME_COL \
#               --pfamdb PFAM_HMM_FILE --output  DOM_HMM_FILE

#-------------------------------------------------------------------------------
import argparse
import pyhmmer
import pandas as pd

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --domains : summary TSV file for domains of interest
       --namecol : column in the domains TSV file that holds domain names
       --pfamdb  : PFAM database HMM domains file
       --output  : HMM output file with profiles describing selected domains
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--domains', type=str, required=True,
        help='Summary TSV file for domains of interest')
    parser.add_argument('--namecol', type=str, required=True,
        help='Column in the domains TSV file that holds domain names')
    parser.add_argument('--pfamdb', type=str, required=True,
        help='PFAM database HMM domains file')
    parser.add_argument('--output', type=str, required=True,
        help='HMM output file with profiles describing selected domains')

    args = parser.parse_args()
    return args

def main():
    '''The entry point function that fetches by their names HMM profiles
       indicated in the input domains table from a PFAM database HMM file,
       and saves them to the otuput HMM file.
    '''
    args = parse_args()
    dom_df = pd.read_csv(args.domains, sep='\t')
    hmm_fin  = pyhmmer.plan7.HMMFile(args.pfamdb)
    hmm_fout = open(args.output, 'wb')
    for hmm in hmm_fin:
        if hmm.name in dom_df[args.namecol].values:
            hmm.write(hmm_fout)
    hmm_fout.close()
    hmm_fin.close()
#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

