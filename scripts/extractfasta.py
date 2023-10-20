#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Based on final, filtered HMMsearch results extracts relevant sequences from
# files with ORF protein sequences (per each assembly accession) located in given
# directory and saves them to the output file. Arguments:
# --seqdir : directory with per assembly ORF seqence data in FASTA format
# --hmmres : final, filtered HMMsearch results in TSV format
# --output : output FASTA file for extracted sequences
# USAGE:
# ./extractfasta.py --seqdir SEQ_DIR --hmmres DOM_HMMOUT --output FASTA_OUT

#-------------------------------------------------------------------------------
import argparse
import pandas as pd
from lib.fasta import read_fasta

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --seqdir : directory with per assembly ORF seqence data in FASTA format
       --hmmres : final, filtered HMMsearch results in TSV format
       --output : output FASTA file for extracted sequences
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--seqdir', type=str, required=True,
        help='Directory with per assembly ORF seqence data in FASTA format')
    parser.add_argument('--hmmres', type=str, required=True,
        help='Final, filtered HMMsearch results in TSV format')
    parser.add_argument('--output', type=str, required=True,
        help='Output FASTA file for extracted sequences')

    args = parser.parse_args()
    return args

def main():
    '''The entry point function that based on final, filtered HMMsearch results
       extracts relevant sequences from files with ORF protein sequences
       (per each assembly accession) located in given directory and saves them
       to the output file.
    '''
    # Load final, filtered HMMsearch results to a DataFrame.
    args = parse_args()
    hmm_df = pd.read_csv(args.hmmres, sep='\t')
    fout = open(args.output, 'w')
    # Group hmm_df by asmacc (assembly accession) to read a given input sequence
    # data only once. iterate over tname, select seuences from seqs dictionary
    # { seqid : sequence }, save them to the output file in FASTA format.
    for asmacc, sub_df in hmm_df[['asmacc', 'tname']].groupby('asmacc'):
        seqs = read_fasta(f'{args.seqdir}/{asmacc}.faa', sub_df['tname'].to_numpy())
        for seqid in sub_df['tname'].sort_values():
            fout.write(seqs[seqid].fasta())
    fout.close()

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

