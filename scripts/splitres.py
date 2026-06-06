#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Extracts a subset of sequences and their annotations, selected by a list of
# protein IDs, from the final integrated sequence and GFF3 annotation files.
# Used to produce per-architecture output files in the last pipeline step.
# Arguments:
# --allseqs   : FASTA file with all final protein sequences
# --allannots : GFF3 file with all annotations
# --seqids    : single-column TSV with protein IDs to extract
# --outseqs   : output FASTA file with selected sequences
# --outannots : output GFF3 file with annotations for selected sequences
#USAGE:
# ./splitres.py --allseqs ALL_SEQS --allannots ALL_ANNOTS --seqids SEQ_IDS \
#               --outseqs OUT_SEQS --outannots OUT_ANNOTS

#-------------------------------------------------------------------------------
import argparse, sys
import pandas as pd
from os import linesep as eol
from lib.fasta import read_fasta

GFF3_HEAD = f'##gff-version 3{eol}'
GFF3_COLS = 'seqid source type start end score strand phase attributes'.split()

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --allseqs   : FASTA file with all final protein sequences
       --allannots : GFF3 file with all annotations
       --seqids    : single-column TSV with protein IDs to extract
       --outseqs   : output FASTA file with selected sequences
       --outannots : output GFF3 file with annotations for selected sequences
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--allseqs', type=str, required=True,
        help='[Lorem ipsum]')
    parser.add_argument('--allannots', type=str, required=True,
        help='[Lorem ipsum]')
    parser.add_argument('--seqids', type=str, required=True,
        help='[Lorem ipsum]')
    parser.add_argument('--outseqs', type=str, required=True,
        help='[Lorem ipsum]')
    parser.add_argument('--outannots', type=str, required=True,
        help='[Lorem ipsum]')

    args = parser.parse_args()
    return args

def main():
    '''Reads the list of protein IDs, extracts matching sequences from the
       input FASTA and matching annotations from the input GFF3, and writes
       them to the respective output files.
    '''
    args = parse_args()
    
    seqids = pd.read_csv(args.seqids).squeeze()
    seqs = read_fasta(args.allseqs, seqids=seqids)
    with open(args.outseqs, 'w') as f:
        for seq in seqs.values():
            f.write(seq.fasta())
            
    df_gff3 = pd.read_csv(args.allannots, comment='#', names=GFF3_COLS, sep='\t')
    df_gff3_sel = df_gff3[df_gff3['seqid'].isin(seqids.values)]
    with open(args.outannots, 'w') as f:
        f.write(GFF3_HEAD)
        df_gff3_sel.to_csv(f, header=False, index=False, sep='\t')

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

