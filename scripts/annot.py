#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Creates a GFF3 file with annotations for the final set of protein sequences
# based on integrated HMMsearch and SignalP results. For each sequence, generates
# GFF3 feature records describing the positions and identities of domain and
# signal peptide hits. Arguments:
# --allres  : integrated HMMsearch and SignalP results TSV
# --seqs    : FASTA file with the final protein sequences to annotate
# --output  : output GFF3 file with domain and signal peptide annotations
# USAGE:
# ./annot.py --allres INTEGRATED_RES --seqs ANALYSED_SEQS --output GFF3_ANNOTS

#-------------------------------------------------------------------------------
import argparse
import pandas as pd
from os import linesep as eol
from lib.fasta import read_fasta

GFF3_HEAD = f'##gff-version 3{eol}'

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --allres  : integrated HMMsearch and SignalP results TSV
       --seqs    : FASTA file with the final protein sequences to annotate
       --output  : output GFF3 file with domain and signal peptide annotations
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--allres', type=str, required=True,
        help='Final integrated HMMsearch and SignalP results')
    parser.add_argument('--seqs', type=str, required=True,
        help='FASTA file with analysed sequences')
    parser.add_argument('--output', type=str, required=True,
        help='GFF3 file with annotations')

    args = parser.parse_args()
    return args

def main():
    '''Loads integrated HMMsearch and SignalP results, filters to sequences
       present in the input FASTA, and writes a GFF3 file with one feature
       record per domain or signal peptide hit, sorted by position.
    '''
    # Parse command line arguments, load analysed sequences from FASTA file.
    args = parse_args()
    seqs = read_fasta(args.seqs)
    
    # Load HMMsearch results
    all_df = pd.read_csv(args.allres, sep='\t')
    all_df = all_df[ all_df['tname'].isin(seqs) ]
    
    # Sort the final DataFrame with respect to domains locations (order domains).
    all_df.sort_values('tname env_from env_to'.split(), inplace=True)
    
    # Process SignalP hits that are and filtered HMM search hits to generate
    # a GFF3 file desribing positions and kids of relevant domains in
    # the analysed sequenes.
    fout = open(args.output, 'w')
    fout.write(GFF3_HEAD)
    for tname, sub_df in all_df.groupby('tname'):
        for _, (qname, qacc, group, srcid, env_from, env_to) in \
            sub_df['qname qacc group srcid env_from env_to'.split()].iterrows():
            name = f'Name={group} ({qname} [{qacc}])'
            line = tname, srcid, group, str(env_from), str(env_to), '.', '+', \
                   '.', name
            fout.write('\t'.join(line) + eol)
    fout.close()

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

