#!/usr/bin/env python
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Searches an HMM profile database against protein sequences using PyHMMER
# hmmsearch. Domain hits are filtered by minimum query (HMM model) coverage
# (--qcovt) and maximum independent E-value (--iE). Retained hits are annotated
# with group names from a domain summary TSV file; hits not present in that file
# are assigned to group 'OTHER'. Results are written to a TSV file. Arguments:
# --domains  : TSV file listing domains of interest with name and group columns
# --namecol  : column in the domains TSV that holds domain names
# --groupcol : column in the domains TSV that holds domain group names
# --hmmdb    : HMM profile database file
# --seqs     : protein sequence FASTA file to search
# --cpus     : number of CPU threads for PyHMMER
# -E         : sequence-level E-value threshold passed to PyHMMER
# --domE     : domain-level E-value threshold passed to PyHMMER
# --incE     : sequence inclusion E-value threshold passed to PyHMMER
# --incdomE  : domain inclusion E-value threshold passed to PyHMMER
# --qcovt    : minimum query coverage fraction for a domain hit to be retained
# --iE       : maximum independent E-value for a domain hit to be retained
# --output   : output TSV file with filtered and annotated domain hits
#USAGE:
# ./hmmsearch.py --domains  DOMAINS_TSV  --namecol  NAME_COL   --groupcol GROUP_COL \
#                --hmmdb    HMM_DB       --seqs     SEQS_FILE  --cpus     N_CPUS    \
#                -E         SEQ_EVALUE   --domE     DOM_EVALUE --incE     INC_EVALUE \
#                --incdomE  INCDOM_EVAL  --qcovt    QCOV_THRES --iE       I_EVALUE  \
#                --output   OUTPUT_TSV

#-------------------------------------------------------------------------------
import argparse
import pyhmmer
import pandas as pd

#-------------------------------------------------------------------------------
def parse_args():
    '''Parses command line arguments:
       --domains  : TSV file listing domains of interest with name and group columns
       --namecol  : column in the domains TSV that holds domain names
       --groupcol : column in the domains TSV that holds domain group names
       --hmmdb    : HMM profile database file
       --seqs     : protein sequence FASTA file to search
       --cpus     : number of CPU threads for PyHMMER
        -E        : sequence-level E-value threshold passed to PyHMMER
       --domE     : domain-level E-value threshold passed to PyHMMER
       --incE     : sequence inclusion E-value threshold passed to PyHMMER
       --incdomE  : domain inclusion E-value threshold passed to PyHMMER
       --qcovt    : minimum query coverage fraction for a domain hit to be retained
       --iE       : maximum independent E-value for a domain hit to be retained
       --output   : output TSV file with filtered and annotated domain hits
       Returns:
       args : ArgumentParser object
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--domains', type=str, required=True,
        help='Summary TSV file for domains of interest')
    parser.add_argument('--namecol', type=str, required=True,
        help='Column in the domains TSV file that holds domain names')
    parser.add_argument('--groupcol', type=str, required=True,
        help='Column in the domains TSV file that indicates domain groups')
    parser.add_argument('--hmmdb', type=str, required=True,
        help='')
    parser.add_argument('--seqs', type=str, required=True,
        help='')
    parser.add_argument('--cpus', type=int, required=True,
        help='')
    parser.add_argument('-E', type=float, required=True,
        help='')
    parser.add_argument('--domE', type=float, required=True,
        help='')
    parser.add_argument('--incE', type=float, required=True,
        help='')
    parser.add_argument('--incdomE', type=float, required=True,
        help='')
    parser.add_argument('--qcovt', type=float, required=True,
        help='')
    parser.add_argument('--iE', type=float, required=True,
        help='')
    parser.add_argument('--output', type=str, required=True,
        help='')

    args = parser.parse_args()
    return args

def main():
    '''
    '''
    args = parse_args()
    hmms = pyhmmer.plan7.HMMFile(args.hmmdb)
    seqs = pyhmmer.easel.SequenceFile(args.seqs, digital=True)
    
    columns = ('tname srcid start end asmacc clustid qname qacc qlen qcovt ' +
               'E-value c-Evalue i-Evalue hmm_from hmm_to env_from env_to').split()
    
    # TODO: write results directly to a file to save RAM
    data = []
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=args.cpus, E=args.E,
        domE=args.domE, incE=args.incE, incdomE=args.incdomE):
        for hit in hits:
            props = dict(pair.split('=') for pair in hit.description.split(' '))
            for domain in hit.domains:
                alignment = domain.alignment
                qcovt = (alignment.hmm_to - alignment.hmm_from + 1) / alignment.hmm_length
                assert qcovt >= 0, f'Negative domain coverage: {domain.name} {qcovt}'
                if qcovt >= args.qcovt and domain.i_evalue <= args.iE:
                    data.append([hit.name, props['srcid'],
                                 props['start'], props['end'], props['asmacc'],
                                 props['clustid'], hits.query.name,
                                 hits.query.accession, hits.query.M,
                                 qcovt, hit.evalue, domain.c_evalue, domain.i_evalue,
                                 alignment.hmm_from, alignment.hmm_to,
                                 domain.env_from, domain.env_to])
    hmm_df = pd.DataFrame(data=data, columns=columns)
    
    # TODO: push this out as an extra step
    # Load to a DataFrame searched domain data, merge with HMMsearch data.
    dom_df = pd.read_csv(args.domains, index_col=args.namecol, sep='\t')
    hmm_df = hmm_df.merge(dom_df[args.groupcol], left_on='qname', right_index=True, how='left')
    hmm_df[args.groupcol] = hmm_df[args.groupcol].fillna('OTHER')
    
    hmm_df.to_csv(args.output, index=False, sep='\t')
    
#-------------------------------------------------------------------------------
# Entry point.
if __name__ == '__main__':
    main()

