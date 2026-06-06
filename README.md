## pfam-genomes

A compact pipeline that allows for searching selected domains (based on [Pfam HMM models](https://www.ebi.ac.uk/interpro/download/pfam)) in proteins encoded in genomic sequences. The pipeline utilises the [PyHMMER](https://pyhmmer.readthedocs.io) library providing bindings to the [HMMER](http://hmmer.org) tools.

### 1. Prerequisites
This pipeline utilises Conda package manager and Conda environments. To set up any Conda environment properly you need to install, for instance, [Miniconda](https://docs.conda.io/en/latest/miniconda.html). The pipeline was tested on Ubuntu 24.04.4 LTS x86-64 using Conda 26.5.2. 

### 2. PFAM database and test genomes
Using `test_data.sh` script, you may download 100 staphylococcal genomic sequences from [NCBI GenBank](https://ftp.ncbi.nlm.nih.gov/genomes) to `genomes` directory. The genome assembly accessions are listes in the `input/test_genomes.txt` file. The script also downloads the [`Pfam-A.hmm.gz`](https://www.ebi.ac.uk/interpro/download/Pfam/) filedatabase. The pipeline utilises a set of 23 CAT (catalytic) and 18 CWT (cell-wall-targeting) domains (`input/domains.tsv`) that occur in bacterial peptidoglycan hydrolases.

The script downloading the test data utilises the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools) tools. Before running it create and activate the `ncbi-cli` environment:
```bash
conda env create --file envs/ncbi-cli.yml
conda activate ncbi-cli
```

Then, run the script in the pipeline directory:
```bash
./test_data.sh
```

### 3. The Snakemake environment
Create the `snakemake` Conda environment necessary to run the pipeline:
```bash
conda env create --file envs/snakemake.yml
```
Next, activate the environment:
```bash
conda activate snakemake
```

### 4. Directory structure and pipeline files
The pipeline utilises the following directory structure:
```
your_pipeline_location/
├── Snakefile
├── config.yaml
├── envs
│   ├── blast.yml
│   ├── ncbi-cli.yml
│   ├── pyhmmer.yml
│   └── snakemake.yml
├── input
│   ├── colors.tsv
│   ├── domains.tsv
│   ├── IUPACDNA.txt
│   ├── TABLE11.txt
│   └── test_genomes.txt
├── src
│   ├── build.sh
│   └── extractorfs.cpp
├── templates
│   ├── style.css
│   └── tmpl.html
├── scripts
│   ├── lib
│   │   └── fasta.py
│   ├── annot.py
│   ├── archchart.py
│   ├── archfilter.py
│   ├── extractfasta.py
│   ├── extractorfs
│   ├── hmmclean.py
│   ├── hmmfetch.py
│   ├── hmmsearch.py
│   ├── join.py
│   ├── sigprocess.py
│   ├── splitfasta.py
│   ├── splitres.py
│   └── unique.py
├── log
└── output
```
In the working directory you can find the `Snakefile` describing the pipeline and `config.yaml`, in which paths to [`Pfam-A.hmm.gz`](https://www.ebi.ac.uk/interpro/download/Pfam) and the directory with genomic sequences as well as domain architectures of interest are provided. In `input/` there are files describing nucleotide sequence alphabet and translation table to be used (`IUPACDNA.txt`, `TABLE11.txt`) as well as groups of domains and domains being searched for next to colors these groups and domains are to be painted with in a visualisation of all possible architectures that are found in searched genomes (`domains.tsv`, `colors.tsv`). In `src/` a source file `extractorfs.cpp` for `scripts/extractorfs` is placed. If `build.sh` is run from that directory, the application will be recompiled and saved as `scripts/extractorfs`. In `templates/` CSS and HTML templates are located. These are used to generate a HTML visualisation of all possible architectures that are found in searched genomes, which is then converted to a PNG file. Necessary scripts and one compiled application are located in `scripts/`. Directories `output/` and `log/` will be created automatically once the pipeline is run. All diagnostic and error messages from tools and scripts used by the pipeline will be redirected to files in the `log/` directory.

### 5. Pipeline architecture
The pipeline described in the Snakefile encompasses the following stages:
1. **extractorfs** -- using `scripts/extractorfs`, provided DNA alphabet (`input/IUPACDNA.txt`) and translation table (`input/TABLE11.txt`), from each genome extracts all possible open reading frames (ORFs) of length >= 300&nbsp;nt.
1. **uniquetrans** -- using `scripts/unique.py`, clusters quickly extracted protein sequences based on their 100% identity to prepare a non-redundant set for HMM searches.
1. **splitfasta** -- splits the unique sequences into chunks, the number of which is equal to the number of cores assigned to the execution of the pipeline.
1. **hmmfetch** -- using `scripts/hmmfetch.py`, fetches domains, from the `Pfam-A.hmm` file, that are listed in `input/domains.tsv` and whose Pfam accession version numbers are provided in the `pfam_acc` column of that file.
1. **hmmselsearch** -- using `scripts/hmmsearch.py`, searches for retrieved domains in the non-redundant protein sequence set chunks and filters the results.
1. **extracttrans** -- using `scripts/extractfasta.py`, extracts sequences indicated by the previous step from the non-redundant set.
1. **hmmallsearch** -- using `scripts/hmmsearch.py`, searches for all domains from the PFAM database in the extracted sequences and filters the results.
1. **signalp** -- using independently installed [SignalP 5.0](https://services.healthtech.dtu.dk/services/SignalP-5.0) tool available as `signalp`, detects N-terminal signal sequences in the extracted sequences. If the tool is not available, the step generates an empty result file and is skipped.
1. **sigprocess** -- using `scripts/sigprocess.py`, filters the SignalP results.
1. **hmmclean** -- using `scripts/hmmclean.py`, removes duplicated signal peptide hits in respect to `signalp` results.
1. **join** -- combines the `hmmsearch` and `signalp` results into one table.
1. **integrate** -- collects result files from the `extracttrans` and `join`, which are run in parallel.
1. **archchart** -- using `scripts/archchart.py`, branches to generate charts in HTML format that describe all domain architectures found, also in regard to groups of domains (column `group` in `input/domains.tsv`).
1. **annotdom** -- using `scripts/annot.py`, prepares a GFF3 file with annotations for the finally selected protein sequences. Annotations are prepared based on `hmmsearch` filtered results and `signalp` results.
1. **archfilter** -- using `scripts/archfilter.py`, continues the main branch and from preprocessed hits select those of independent E-value (i-Evalue) <= 0.001 and domain coverage >= 80% (0.8). Next from target protein sequences select those with domain architecture of interest (`config.yaml`).
1. **split** -- splits the filtered results per architecture listed in the `config.yaml` file.

More detailed description on how the pipeline works you will find in comments in the `Snakefile`, `config.yaml` file and the script/source files.

### 6. Running the pipeline
Providing that you have `snakemake` environment properly set up and activated as well as at least the test data downloaded, you can first look up the list of tasks to be done by running the following command from the pipeline directory, where `Snakefile` is located:
```bash
snakemake --dryrun --quiet rules
```
and then run the pipeline using as many cores as you wish. Indicate that the pipeline utilises `conda` directives in order to create Conda environments for running the pipeline tasks. We strongly recommend using the greedy task scheduler:
```bash
snakemake --cores number_of_cores --scheduler greedy --use-conda
```
Final data, non-redundant protein sequences and corresponding annotations for each domain architecture of interest, as well as HTML visualisations of all domain architectures found, you will find in `output/final` directory.

