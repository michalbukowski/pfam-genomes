import os
import pandas as pd

configfile: 'config.yaml'
gen_dir = config['gen_dir']
pfam_db = config['pfam_db']
archs   = config['archs']
max_cores = os.cpu_count()

assemblies, = glob_wildcards(gen_dir + '/{assembly}_genomic.fna.gz')

rule all:
    input:
        expand('output/final/final_{arch}.gff3', arch=archs.keys())

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
                                 >      {log} 2>&1
        '''

rule uniquetrans:
    params:
        mask = rules.extractorfs.output.trans.replace('{assembly}', '*')
    input:
        expand(rules.extractorfs.output.trans, assembly=assemblies)
    output:
        repr  = 'output/unique/all_assembly.faa',
        clust = 'output/unique/all_assembly.faa.tsv'
    log:
        'log/uniquetrans.log'
    shell:
        '''fpaths=({params.mask})
           for fpath in "${{fpaths[@]}}"; do
               cat "${{fpath}}"
           done                                         \
           |                                            \
           scripts/uniqueseqs.py --output {output.repr} \
                                   >  {log} 2>&1
        '''

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
                                 --domdata {input}      \
                                   2>      {log}        \
           |                                            \
           hmmfetch               -f                    \
                                  -o       {output}     \
                                           {pfam_db}    \
                                           -            \
                                   >       {log} 2>&1
        '''

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

rule filter:
    params:
        qcovt   = 0.8,
        iEvalue = 0.001,
        cols    = 'tname tacc tlen qname qacc qlen E-value seqscore seqbias '+ \
                  '# of c-Evalue i-Evalue domscore dombias hmm_from hmm_to ' + \
                  'ali_from ali_to env_from env_to acc desc',
        leave   = 'tname srcid start end asmacc qname qacc qlen qcovt ' + \
                  'group E-value c-Evalue i-Evalue hmm_from hmm_to '    + \
                  'env_from env_to',
        arch    = lambda wildcards: archs[wildcards.arch]
    input:
        domdata = rules.hmmfetch.input,
        domtbl  = rules.hmmsearch.output
    output:
        'output/filter/filter_{arch}.tsv'
    log:
        'log/filter/filter_{arch}.log'
    shell:
        '''scripts/filter.py --qcovt    {params.qcovt}   \
                             --iEvalue  {params.iEvalue} \
                             --cols    "{params.cols}"   \
                             --leave   "{params.leave}"  \
                             --arch    '{params.arch}'   \
                             --domdata  {input.domdata}  \
                             --domtbl   {input.domtbl}   \
                             --output   {output}         \
                               > {log} 2>&1
        '''

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
        '''if [[ ! -s {input} ]]; then
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
        '''if [[ ! -s {input.sigres} ]]; then
               touch {output}
               exit
           fi
           scripts/annot.py --sigres {input.sigres} \
                            --hmmres {input.hmmres} \
                            --seqs   {input.seqs}   \
                            --output {output}       \
                              > {log} 2>&1
        '''

