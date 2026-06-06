# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

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
# GRP - group names
# NAM - PFAM domain names
# values of both must be provided in input/domains.tsv
# Examples:
# 'GRP|CAT-CWT' - search for 'CAT-CWT' architecture: only two out of all searched
#                 domains are present, a domain from CAT group preceeds a domain
#                 from CWT group
# 'ACC|.*SH3_5.*' - search for any protein containing SH3_5 domain
#                   that may be surounded by any other domain
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
        expand('output/final/final_{arch}.faa',  arch=archs),
        expand('output/final/final_{arch}.gff3', arch=archs),
        'output/final/architectures.html'

# In the 1st step, extract from each genome all posisble open reading frames (ORFs)
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

# In the 2nd step, quickly cluster sequences based on their 100% identity to
# prepare a non-redundant set for HMM searches.
rule uniquetrans:
    conda:
        'envs/pyhmmer.yml'
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

# In the 3rd step, split the unique protein sequences into chunks
# for parallel processing.
checkpoint splitfasta:
    threads:
        max_cores
    input:
        rules.uniquetrans.output.repr
    output:
        directory('output/splitfasta')
    log:
        'log/splitseqs.log'
    shell:
        '''rm -rf   {output}
           mkdir -p {output}
           scripts/splitfasta.py --chunks {threads} \
                                 --input  {input}   \
                                 --outdir {output}  \
                                   > {log} 2>&1
        '''


# In the 4th step, from Pfam-A HMM file fetch domains that are listed in
# input/domains.tsv by their names.
rule hmmfetch:
    conda:
        'envs/pyhmmer.yml'
    params:
        namecol = 'pfm_name'
    input:
        domains = 'input/domains.tsv',
        pfamdb  = pfam_db
    output:
        'output/hmm/domains.hmm'
    log:
        'log/hmm/hmmfetch.log'
    shell:
        '''scripts/hmmfetch.py --domains {input.domains}  \
                               --namecol {params.namecol} \
                               --pfamdb  {input.pfamdb}   \
                               --output  {output}         \
                                 >>      {log} 2>&1
        '''

# In the 5th step, search for domains retrived from Pfam-A HMM file in the
# non-redundant protein sequence set.
rule hmmselsearch:
    conda:
        'envs/pyhmmer.yml'
    params:
        namecol  = rules.hmmfetch.params.namecol,
        groupcol = 'group',
        E        = 0.1,
        domE     = 0.1,
        incE     = 0.001,
        incdomE  = 0.001,
        qcovt    = 0.8,
        iE       = 0.1
    threads:
        1
    input:
        domains = rules.hmmfetch.input.domains,
        hmmdb   = rules.hmmfetch.output,
        trans   = rules.splitfasta.output[0] + '/{chunk}.faa'
    output:
        'output/hmm/hmmselsearch_{chunk}.tsv'
    log:
        'log/hmm/hmmselsearch_{chunk}.log'
    shell:
        '''scripts/hmmsearch.py --domains  {input.domains}   \
                                --namecol  {params.namecol}  \
                                --groupcol {params.groupcol} \
                                --hmmdb    {input.hmmdb}     \
                                --seqs     {input.trans}     \
                                --cpus     {threads}         \
                                 -E        {params.E}        \
                                --domE     {params.domE}     \
                                --incE     {params.incE}     \
                                --incdomE  {params.incdomE}  \
                                --qcovt    {params.qcovt}    \
                                --iE       {params.iE}       \
                                --output   {output}          \
                                   >>      {log} 2>&1
        '''


# In the 6th step, extract the protein seqences for which input
# domain matches were found.
rule extracttrans:
    conda:
        'envs/pyhmmer.yml'
    params:
        seqdir = os.path.dirname(rules.extractorfs.output.trans)
    input:
        rules.hmmselsearch.output
    output:
        'output/hmm/extracttrans_{chunk}.faa'
    log:
        'log/hmm/extracttrans_{chunk}.log'
    shell:
        '''scripts/extractfasta.py --seqdir {params.seqdir} \
                                   --hmmres {input}         \
                                   --output {output}        \
                                     > {log} 2>&1
        '''


# In the 7th step, search the extracted protein sequences for all
# domains from the PFAM database.
rule hmmallsearch:
    conda:
        'envs/pyhmmer.yml'
    params:
        namecol  = rules.hmmfetch.params.namecol,
        groupcol = 'group',
        E        = 0.1,
        domE     = 0.1,
        incE     = 0.001,
        incdomE  = 0.001,
        qcovt    = 0.8,
        iE       = 0.1
    threads:
        1
    input:
        domains = rules.hmmfetch.input.domains,
        hmmdb   = rules.hmmfetch.input.pfamdb,
        trans   = rules.extracttrans.output
    output:
        'output/hmm/hmmallsearch_{chunk}.tsv'
    log:
        'log/hmm/hmmallsearch_{chunk}.log'
    shell:
        '''if [[ ! -s {input.trans} ]]; then
               touch {output}
               exit
           fi
           scripts/hmmsearch.py --domains  {input.domains}   \
                                --namecol  {params.namecol}  \
                                --groupcol {params.groupcol} \
                                --hmmdb    {input.hmmdb}     \
                                --seqs     {input.trans}     \
                                --cpus     {threads}         \
                                 -E        {params.E}        \
                                --domE     {params.domE}     \
                                --incE     {params.incE}     \
                                --incdomE  {params.incdomE}  \
                                --qcovt    {params.qcovt}    \
                                --iE       {params.iE}       \
                                --output   {output}          \
                                   >>      {log} 2>&1
        '''


# In the 8th step, continuing step 6th, use SignalP to detect N-terminal signal
# sequences in the final set of protein sequences for each domain architecture
# of interest. This step requires a separate SignalP installation and an access
# to it via signalp command. If the command is not found, empty output file
# is generated.
rule signalp:
    threads:
        1
    params:
        orgs   = 'gram+',   # ('gram+', 'gram-'),
        format = 'short',
        plot   = 'none'
    input:
        rules.extracttrans.output
    output:
        'output/signalp/signalp_{chunk}.tsv'
    log:
        'log/signalp/signalp_{chunk}.log'
    shell:
        '''if [[ ! -s {input} || $(command -v signalp) == '' ]]; then
               touch {output}
               exit
           fi
           rm -f {output} {log}
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


# In the 9th step, process the signalp results.
rule sigprocess:
    conda:
        'envs/pyhmmer.yml'
    params:
        pthresh = 0.1
    input:
        rules.signalp.output
    output:
        'output/signalp/processed_{chunk}.tsv'
    log:
        'log/signalp/sigprocess_{chunk}.log'
    shell:
        '''scripts/sigprocess.py --pthresh {params.pthresh} \
                                 --sigres  {input}          \
                                 --output  {output}         \
                                   >       {log} 2>&1
        '''


# In the 10th step, remove signal peptides duplicate identifications
# from the hmmallsearch results.
rule hmmclean:
    conda:
        'envs/pyhmmer.yml'
    input:
        sigres = rules.sigprocess.output,
        hmmres = rules.hmmallsearch.output
    output:
        'output/hmm/hmmcleaned_{chunk}.tsv'
    log:
        'log/hmm/hmmclean_{chunk}.log'
    shell:
        '''if [[ ! -s {input.hmmres} ]]; then
               touch {output}
               exit
           fi
           scripts/hmmclean.py --sigres {input.sigres} \
                               --hmmres {input.hmmres} \
                               --output {output}       \
                                 > {log} 2>&1
        '''


# In the 11th step, combine the hmmallsearch and signalp results.
rule join:
    conda:
        'envs/pyhmmer.yml'
    input:
        sigres = rules.sigprocess.output,
        hmmres = rules.hmmclean.output
    output:
        'output/join/joined_{chunk}.tsv'
    log:
        'log/join/join_{chunk}.log'
    shell:
        '''if [[ ! -s {input.hmmres} ]]; then
               touch {output}
               exit
           fi
           scripts/join.py --sigres {input.sigres} \
                           --hmmres {input.hmmres} \
                           --output {output}       \
                             > {log} 2>&1
        '''


# Run the splitfasta checkpoint and generate the list of files expected
# as the output from extracttrans.
def alltrans_agg(wildcards):
    checkpoint_output = checkpoints.splitfasta.get().output[0]
    chunks, = glob_wildcards(os.path.join(checkpoint_output, "{chunk}.faa"))
    return expand(rules.extracttrans.output, chunk=chunks)


# Run the splitfasta checkpoint and generate the list of files expected
# as the output from join. 
def allres_agg(wildcards):
    checkpoint_output = checkpoints.splitfasta.get().output[0]
    chunks, = glob_wildcards(os.path.join(checkpoint_output, "{chunk}.faa"))
    return expand(rules.join.output, chunk=chunks)


# In the 12th step, collect the files from extracttrans and join that were
# generated as parallel output for unique sequences chunks.
rule integrate:
    input:
        alltrans = alltrans_agg,
        allres   = allres_agg
    output:
        alltrans = 'output/integrated.faa',
        allres   = 'output/integrated.tsv'
    log:
        'log/integrate.log'
    shell:
        '''cat {input.alltrans} > {output.alltrans} 2> {log}
           files=({input.allres})
           head -1 ${{files[0]}} > {output.allres} 2>> {log}
           for file in ${{files[@]}}; do
               tail -n +2 ${{file}} >> {output.allres} 2>> {log}
           done
        '''


# In the 13B-th step, generate charts in HTML that describe
# all domain architectures, also in regard to groups of domains.
rule archchart:
    conda:
        'envs/pyhmmer.yml'
    input:
        style   = 'templates/style.css',
        tmpl    = 'templates/tmpl.html',
        colors  = 'input/colors.tsv',
        clstlen = rules.uniquetrans.output.length,
        allres  = rules.integrate.output.allres
    output:
        'output/final/architectures.html'
    log:
        'log/final/archchart.log'
    shell:
        '''scripts/archchart.py --style   {input.style}   \
                                --tmpl    {input.tmpl}    \
                                --colors  {input.colors}  \
                                --clstlen {input.clstlen} \
                                --allres  {input.allres}  \
                                --output  {output}        \
                                  > {log} 2>&1
        '''


# In the 13th step, prepare GFF3 file with annotations for the final protein
# sequence set. Annotations are prepared based on HMMsearch filtered results
# obtained for a search for all domains as well as SignalP results.
rule annotdom:
    conda:
        'envs/pyhmmer.yml'
    input:
        allres = rules.integrate.output.allres,
        seqs   = rules.integrate.output.alltrans
    output:
        rules.integrate.output.alltrans[:rules.integrate.output.alltrans.rfind('.')] + '.gff3'
    log:
        'log/final/annotdom.log'
    shell:
        '''scripts/annot.py --allres {input.allres} \
                            --seqs   {input.seqs}   \
                            --output {output}       \
                              > {log} 2>&1
        '''


# In the 14th step filter the results, leave hits of independent E-value (i-Evalue)
# <= 0.001 and domain coverage >= 80% (0.8). Next select sequences with domain
# architecure of interest.
rule archfilter:
    conda:
        'envs/pyhmmer.yml'
    params:
        arch = lambda wildcards: archs[wildcards.arch]
    input:
        allres = rules.integrate.output.allres
    output:
        'output/archfilter/archfilter_{arch}.tsv'
    log:
        'log/archfilter/archfilter_{arch}.log'
    shell:
        '''scripts/archfilter.py --arch    '{params.arch}' \
                                 --allres   {input.allres} \
                                 --output   {output}       \
                                   >        {log} 2>&1
        '''


# In the 15th step, split the final results per domain architecture
# among those indicated in the config.yml file.
rule splitres:
    input:
        seqs   = rules.integrate.output.alltrans,
        annots = rules.annotdom.output,
        seqids = rules.archfilter.output
    output:
        seqs   = 'output/final/final_{arch}.faa',
        annots = 'output/final/final_{arch}.gff3'
    log:
        'log/final/splitres_{arch}.log'
    shell:
        '''scripts/splitres.py --allseqs   {input.seqs}    \
                               --allannots {input.annots}  \
                               --seqids    {input.seqids}  \
                               --outseqs   {output.seqs}   \
                               --outannots {output.annots} \
                                 > {log} 2>&1
        '''

