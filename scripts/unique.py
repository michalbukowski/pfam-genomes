#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Performs fast dictionary-based clustering of sequences in order to obtain
# a non-redundant set. Writes clustering information to an extra TSV file
# by changing the suffix of the FASTA output file name to '_clusts.tsv', and
# for cluster item counts to '_lenghts.tsv'. Arguments:
# --input  : protein FASTA file with sequences to be clustered
# --output : protein FASTA file with representative sequences
# USAGE:
# ./uniqueseqs.py [--input INPUT_FASTA] --output OUTPUT_FASTA

#-------------------------------------------------------------------------------
import argparse, sys
from os import linesep as eol
from os.path import sep, extsep
from glob import glob
from itertools import count
from lib.fasta import fasta, fasta_meta

#-------------------------------------------------------------------------------
# A helper class for cluster data storage: cluster id and cluster sequence count.
class Clust:
    def __init__(self, clustid, length):
        self.clustid = clustid
        self.length  = length

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --input  : protein FASTA file with sequences to be clustered,
                  wildcard * is allowed for integrating data from multiple files
       --output : protein FASTA file with representative sequences
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input', type=str,
        help='Path to a FASTA file with sequences for clustering, ' +
             'wildcard * is allowed for integrating data from multiple files.')
    parser.add_argument('--output', type=str, required=True,
        help='FASTA file with representative sequences')
    
    args = parser.parse_args()
    sep_pos      = args.output.rfind(sep)
    extsep_pos   = args.output.rfind(extsep)
    pos = extsep_pos if extsep_pos > sep_pos else None
    args.outbase = args.output[:pos]
    
    return args

def main():
    '''The entry point function that performs fast dictionary-based clustering
       of sequences in order to obtain a non-redundant set. Writes clustering
       information to an extra TSV file by changing the suffix of the FASTA
       output file name to '_clusts.tsv', and for cluster item counts to
       '_lenghts.tsv'
    '''
    # Parse command line arguments, create an open range integer generator for
    # naming subsequent clusters and a dictionary for representative sequecnes
    # {sequence : cluster id}
    args = parse_args()
    ngen = count()
    clusts = {}
    
    # Open input and output files, write column names to TSV files that
    # hold clustering information (cluster id and metadata of subsequent
    # sequences) and sequence counts (cluster id and sequecne count).
    # If input file path and path with wildcards are not provided,
    # data are read from stdin.
    fpaths = glob(args.input)
    frepr  = open(args.output, 'w')
    fclust = open(args.outbase + '_clusts.tsv', 'w')
    fcount = open(args.outbase + '_lengths.tsv', 'w')
    
    meta_keys = 'seqid asmacc srcid start end'.split(' ')
    fclust.write('clustid\t' + '\t'.join(meta_keys) + eol)
    col_names = 'clustid length'.split(' ')
    fcount.write('\t'.join(col_names) + eol)
    
    # A helper function for the next block.
    def process_seq():
        seqid  = header.split(' ')[0]
        seqkey = seq[1:]
        if seqkey not in clusts:
            clusts[seqkey] = Clust(seqid, 1)
            frepr.write(fasta(header, '', seq))
        else:
            clusts[seqkey].length += 1
        clustid = clusts[seqkey].clustid
        meta = fasta_meta(header, meta_keys)
        meta = '\t'.join( meta[key] for key in meta_keys )
        fclust.write(f'{clustid}\t{meta}{eol}')
    
    # Read input FASTA file(s) line by line to extract subsequent sequences and their
    # metadata. Use the helper function process_seq() when a sequence is completly
    # extracted and is not a key in the seqs dict (except the first letter), add
    # it to the dict with a new cluster id and write to the output FASTA file
    # as representative. For each sequence write to the TSV file to which
    # cluster it belongs.
    for fpath in fpaths:
        fin = open(fpath)
        header, seq = None, ''
        for line in fin:
            if line[0] == '>':
                line = line[1:].rstrip(eol)
                if header is not None and seq != '' and 'X' not in seq:
                    process_seq()
                header = line
                seq = ''
            elif header is not None:
                seq += line.rstrip(eol)
        # Process the last seq at the end of the file.
        if header is not None and seq != '' and 'X' not in seq:
            process_seq()
        fin.close()
    
    # Save cluster sequence counts to a TSV file.
    lines = [ None ] * len(clusts)
    for i, clust in enumerate(clusts.values()):
        lines[i] = f'{clust.clustid}\t{clust.length}{eol}'
    fcount.writelines(lines)
    
    # Close input and output files.
    fcount.close()
    fclust.close()
    frepr.close()

#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

