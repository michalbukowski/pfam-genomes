# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# A workflow designed for high-throughput searches for proteins of given domain
# architectures encoded in genomic FASTA nucleotide sequences. The initial input,
# next to genomic sequences, encompassess the following files:
# input/domains.tsv  - a TSV file that must contain at lest 2 named columns:
#     group   - domain arbitrary groups, in the demo file these are CAT (catalytic)
#               and CWT (cell wall targeting) groups of domains
#     pfm_acc - Pfam accession version number for each domain
# input/IUPACDNA.txt - a file with DNA alphabet to be used for reverse-complement
#                      searches (first line - all possible characters, second -
#                      complementary chracters), here ambiguous IUPACDNA that
#                      covers all possible characters in DNA sequecnes
# input/TABLE11.txt  - translation table, here 11 (Bacteria), in a short format
#                      as available on https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# The workflow also requires a HMM file with Pfam-A models (see the next part)

import os
import pandas as pd

# Load environmental variables from config.yaml config file, these are:
# gendir  - a path do the directory where FASTA nucleotide format genomic sequences
#           are located, GZIP compressed, of the following name pattern:
#           {assembly}_genomic.fna.gz, where assembly is assembly accession number
# pfam_db - a path to HMM file with Pfam-A models, may be obtained from
#           https://www.ebi.ac.uk/interpro/download/Pfam/
# archs   - a dictionary with custom architecure names and regex describing
#           architectures of interest
# Each architecture regex is preceeded by ACC| or GRP| prefixes, which are stripped
# from the regex before further processing, that inform what kind of references
# to domains are used:
# ACC - Pfam accession version numbers
# GRP - group names
# values of both must be provided in input/domains.tsv
# Examples:
# 'GRP|CAT-CWT' - search for 'CAT-CWT' architecture: only two out of all searched
#                 domains are present, a domain from CAT group preceeds a domain
#                 from CWT group
# 'ACC|.*PF08460\.13.*' - search for any protein containing PF08460.13 domain
#                         that may be surounded by any other domain
configfile: 'config.yaml'
gen_dir = config['gen_dir']
pfam_db = config['pfam_db']
archs   = config['archs']
# Set the maximal number of threads to the number of CPUs in the system.
max_cores = os.cpu_count()

# Retrieve assembly accession numbers of genomes to be analysed.
assemblies, = glob_wildcards(gen_dir + '/{assembly}_genomic.fna.gz')

# The final output of the main workflow branch for each domain architecure is:
# 1) a protein FASTA file that contains relevant non-redundant set of protein
#    sequences derived from the analysed genomes that agree with the domain architecture
# 2) a GFF3 file with annotation for the aforementioned sequences, the annotations
#    describes location and the kind of the domains as well as N-terminal signal peptides
# Additionally the workflow branch B produces:
# 3) a HTML and PNG files with charts depicting all domain architectures found
#    for all query domains found in searched sequences
rule all:
    input:
        expand('output/final/final_{arch}.faa',  arch=archs.keys()),
        expand('output/final/final_{arch}.gff3', arch=archs.keys()),
        'output/final/architectures.png'

# In the 1st step extract from each genome all posisble open reading frames (ORFs)
# of lenght >= 200 nt, based on provided DNA alphabet and translation table.
# For more information see comments in src/extractorfs.cpp file.
rule extractorfs:
    params:
        minlen = 300
    input:
        alph   = 'input/IUPACDNA.txt',
        tab    = 'input/TABLE11.txt',
        genome = gen_dir + '/{assembly}_genomic.fna.gz'
    output:
        seqs  = 'output/extractorfs/fna/{assembly}.fna',
        trans = 'output/extractorfs/faa/{assembly}.faa'
    log:
        'log/extractorfs/{assembly}.log'
    shell:
        '''gunzip               -c      {input.genome}       \
                                 2>     {log}                \
           |                                                 \
           scripts/extractorfs --alph   {input.alph}         \
                               --tab    {input.tab}          \
                               --asmacc {wildcards.assembly} \
                               --minlen {params.minlen}      \
                               --seqs   {output.seqs}        \
                               --trans  {output.trans}       \
                                 >>     {log} 2>&1
        '''

# In the 2nd step quickly cluster sequences based on their 100% identity to
# prepare a non-redundant set for HMM searches.
rule uniquetrans:
    params:
        mask = rules.extractorfs.output.trans.replace('{assembly}', '*')
    input:
        expand(rules.extractorfs.output.trans, assembly=assemblies)
    output:
        repr   = 'output/unique/all_assembly.faa',
        clust  = 'output/unique/all_assembly_clusts.tsv',
        length = 'output/unique/all_assembly_lengths.tsv'
    log:
        'log/uniquetrans.log'
    shell:
        '''scripts/unique.py --input "{params.mask}" \
                             --output {output.repr}  \
                               > {log} 2>&1
        '''

# In the 3rd step, from Pfam-A HMM file fetch domains that are listed in
# input/domains.tsv. Here only pfam_acc column is used.
rule hmmfetch:
    params:
        col = 'pfm_acc'
    input:
        'input/domains.tsv'
    output:
        'output/hmm/domains.hmm'
    log:
        'log/hmm/hmmfetch.log'
    shell:
        '''scripts/extractpfm.sh --column  {params.col} \
                                 --input   {input}      \
                                   2>      {log}        \
           |                                            \
           hmmfetch               -f                    \
                                  -o       {output}     \
                                           {pfam_db}    \
                                           -            \
                                   >>      {log} 2>&1
        '''

# In the 4th step search for domains retrived from Pfam-A HMM file in the
# non-redundant protein sequence set.
rule hmmsearch:
    params:
        E       = 0.01,
        domE    = 0.01,
        incE    = 0.001,
        incdomE = 0.001
    threads:
        max_cores
    input:
        domains = rules.hmmfetch.output,
        trans   = rules.uniquetrans.output.repr
    output:
        'output/hmm/hmmsearch.txt'
    log:
        'log/hmm/hmmsearch.log'
    shell:
        '''hmmsearch --cpu       {threads}        \
                      -E         {params.E}       \
                     --domE      {params.domE}    \
                     --incE      {params.incE}    \
                     --incdomE   {params.incdomE} \
                     --noali                      \
                     --notextw                    \
                     --acc                        \
                      -o         /dev/null        \
                     --domtblout {output}         \
                                 {input.domains}  \
                                 {input.trans}
        '''

# In the 5th step, preprocess the raw HMMsearch results and save relevant
# columns (leave) in the final TSV file.
rule preprocess:
    params:
        qcovt   = 0.8,
        iEvalue = 0.001,
        cols    = 'tname tacc tlen qname qacc qlen E-value seqscore seqbias '+ \
                  '# of c-Evalue i-Evalue domscore dombias hmm_from hmm_to ' + \
                  'ali_from ali_to env_from env_to acc desc',
        leave   = 'tname srcid start end asmacc clustid qname qacc qlen qcovt ' + \
                  'group E-value c-Evalue i-Evalue hmm_from hmm_to '            + \
                  'env_from env_to'
    input:
        domdata = rules.hmmfetch.input,
        domtbl  = rules.hmmsearch.output
    output:
        'output/hmm/preprocessed.tsv'
    log:
        'log/hmm/preprocess.log'
    shell:
        '''scripts/preprocess.py --qcovt    {params.qcovt}   \
                                 --iEvalue  {params.iEvalue} \
                                 --cols    "{params.cols}"   \
                                 --leave   "{params.leave}"  \
                                 --domdata  {input.domdata}  \
                                 --domtbl   {input.domtbl}   \
                                 --output   {output}         \
                                   > {log} 2>&1
        '''

# In the 6-Bth step, branch to generate charts in HTML that describe
# all domain architectures, also in regard to groups of domains
rule archchart:
    input:
        style   = 'templates/style.css',
        tmpl    = 'templates/tmpl.html',
        colors  = 'input/colors.tsv',
        clstlen = rules.uniquetrans.output.length,
        hmmres  = rules.preprocess.output
    output:
        'output/final/architectures.html'
    log:
        'log/archchart.log'
    shell:
        '''scripts/archchart.py --style   {input.style}   \
                                --tmpl    {input.tmpl}    \
                                --colors  {input.colors}  \
                                --clstlen {input.clstlen} \
                                --hmmres  {input.hmmres}  \
                                --output  {output}        \
                                  > {log} 2>&1
        '''

# In the 7-Bth step, continue the branch to convert charts in HTML format to PNG
rule convertchart:
    params:
       quality = 100,
       width   = 2000,
       zoom    = 3
    input:
        rules.archchart.output
    output:
        'output/final/architectures.png'
    log:
        'log/convertchart.log'
    shell:
        '''wkhtmltoimage --disable-smart-width      \
                         --width {params.width}     \
                         --zoom {params.zoom}       \
                         --quality {params.quality} \
                           {input} {output}         \
                           >> {log} 2>&1
        '''

# In the 6th step filter the results, leave hits of independent E-value (i-Evalue)
# <= 0.001 and domain coverage >= 80% (0.8). Next select sequences with domain
# architecure of interest.
rule filter:
    params:
        arch = lambda wildcards: archs[wildcards.arch]
    input:
        rules.preprocess.output
    output:
        'output/filter/filtered_{arch}.tsv'
    log:
        'log/filter/filter_{arch}.log'
    shell:
        '''scripts/filter.py --arch    '{params.arch}'   \
                             --hmmres   {input}          \
                             --output   {output}         \
                               > {log} 2>&1
        '''

# In the 7th step, extract relevant sequences from ORF set for a given
# assembly accession, the value which is retrieved from sequnce metadata
# (FASTA header) preserved in filtered HMMsearch results (desc column).
rule extracttrans:
    params:
        seqdir = os.path.dirname(rules.extractorfs.output.trans)
    input:
        rules.filter.output
    output:
        'output/final/final_{arch}.faa'
    log:
        'log/final/extracttrans_{arch}.log'
    shell:
        '''scripts/extractfasta.py --seqdir {params.seqdir} \
                                   --hmmres {input}         \
                                   --output {output}        \
                                     > {log} 2>&1
        '''

# In the 8th step, continuing step 6th, use SignalP to detect N-terminal signal
# sequences in the final set of protein sequences for each domain architecture
# of interest. This step requires a separate SignalP installation and an access
# to it via signalp command. If the command is not found, empty output file
# is generated.
rule signalp:
    threads:
        max_cores
    params:
        orgs   = ('gram+', 'gram-'),
        format = 'short',
        plot   = 'none'
    input:
        rules.extracttrans.output
    output:
        'output/signalp/signalp_{arch}.tsv'
    log:
        'log/signalp/signalp_{arch}.log'
    shell:
        '''if [[ ! -s {input} || $(command -v signalp) == '' ]]; then
               touch {output}
               exit
           fi
           rm -f {input}, {log}
           for org in {params.orgs}; do
               signalp -batch    {threads}       \
                       -org    "${{org}}"        \
                       -format   {params.format} \
                       -plot     {params.plot}   \
                       -stdout                   \
                       -fasta    {input}         \
                         >>      {output}        \
                        2>>      {log}
           done
        '''

# In the last, 9th step prepare GFF3 file with annotations for the final protein
# sequence set. Annotations are prepared based on HMMsearch filtered results and
# SignalP results.
rule annotdom:
    input:
        sigres = rules.signalp.output,
        hmmres = rules.filter.output,
        seqs   = rules.extracttrans.output
    output:
        rules.extracttrans.output[0][:rules.extracttrans.output[0].rfind('.')] + '.gff3'
    log:
        'log/final/annotdom_{arch}.log'
    shell:
        '''scripts/annot.py --sigres {input.sigres} \
                            --hmmres {input.hmmres} \
                            --seqs   {input.seqs}   \
                            --output {output}       \
                              > {log} 2>&1
        '''

