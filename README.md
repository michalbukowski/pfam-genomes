## pfam-genomes

A compact pipeline that allows for searching selected domains (based on [Pfam HMM models](https://www.ebi.ac.uk/interpro/download/pfam)) in proteins encoded in genomic sequences.

### 1. Conda environment
To set up the environment properly you need to have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed. Then, from the pipeline directory, you need to create an environment `pfam-gen` described in `envs/pfam-gen.txt`:
```bash
conda create -n pfam-gen -f envs/pfam-gen.txt
```
Next, activate the environment:
```bash
conda activate pfam-gen
```
### 2. Test data
Using `test_data.sh` script, you may download 20 staphylococcal genomic sequences from [NCBI GenBank](https://ftp.ncbi.nlm.nih.gov/genomes) to `genomes` directory, as well as [`Pfam-A.hmm`](https://www.ebi.ac.uk/interpro/download/Pfam/) file to test the pipeline on a set of 10 CAT (catalytic) and 5 CWT (cell-wall-targeting) domains (`input/domains.tsv`) that occurr in bacterial peptidoglycan hydrolases. Simply run the script in the pipeline directory:
```bash
./test_data.sh
```
### 3. Directory structure and pipeline files
The pipeline utilises the following directory structure:
```
your_pipeline_location/
├── Snakefile
├── config.yaml
├── input
│   ├── colors.tsv
│   ├── domains.tsv
│   ├── IUPACDNA.txt
│   └── TABLE11.txt
├── src
│   ├── build.sh
│   └── extractorfs.cpp
├── templates
│   ├── style.css
│   └── tmpl.html
├── scripts
│   ├── annot.py
│   ├── archchart.py
│   ├── extractfasta.py
│   ├── extractorfs
│   ├── extractpfm.sh
│   ├── filter.py
│   ├── lib
│   │   └── fasta.py
│   ├── preprocess.py
│   └── unique.py
├── log
└── output
```
In the working directory you can find the `Snakefile` describing the pipeline and `config.yaml`, in which paths to [`Pfam-A.hmm`](https://www.ebi.ac.uk/interpro/download/Pfam) and the directory with genomic sequences as well as domain architectures of interes are provided. In `input/` there are files describing nucleotide sequence alphabet and translation table to be used (`IUPACDNA.txt`, `TABLE11.txt`), groups of domains and domains being sarched for and colors these groups and domains are to be painted with in a visualisation of all possible architectures that are found in searched genomes (`domains.tsv`, `colors.tsv`). In `src/` a source file `extractorfs.cpp` for `scripts/extractorfs` is placed. If `build.sh` is run from that directory, the application will be recompiled and saved as `scripts/extractorfs`. In `templates/` CSS and HTML templates are located. These are used to generate a HTML visualisation of all possible architectures that are found in searched genomes, which next is converted to a PNG file. Necessary scripts and one compiled application are located in `scripts/`. Directories `output/` and `log/` will be created automatically once the pipeline is run. All diagnostic and error messages from tools and scripts used by the pipeline will be redirected to files in the `log/` directory.

### 4. Pipeline architecture
The pipline described in the Snakefile encompasses the following stages:
1. **extractorfs** -- using `scripts/extractorfs`, provided DNA alphabet (`input/IUPACDNA.txt`) and translation table (`input/TABLE11.txt`), from each genome extract all posisble open reading frames (ORFs) of lenght >= 200&nbsp;nt.
1. **uniquetrans** -- using `scripts/unique.py`, cluster quickly extracted non-redundant protein sequences based on their 100% identity to prepare a non-redundant set for HMM searches.
1. **hmmfetch** -- using `hmmfetch` tool, fetch domains, from `Pfam-A.hmm` file, that are listed in `input/domains.tsv` based on `pfam_acc` column.
1. **hmmsearch** -- using `hmmsearch` tool, search for retrived domains in the non-redundant protein sequence set.
1. **preprocess** -- using `scripts/preprocess.py`, preprocess the raw `hmmsearch` results and save relevant columns.
1. **archchart** -- using `scripts/archchart.py`, branch to generate charts in HTML format that describe all domain architectures found, also in regard to groups of domains (column `group` in `input/domains.tsv`).
1. **convertchart** -- using `wkhtmltoimage` tool, continue the branch to convert charts in HTML format to PNG.
1. **filter** -- using `scripts/filter.py`, continue the main branch and from preprocessed hits select those of independent E-value (i-Evalue) <= 0.001 and domain coverage >= 80% (0.8). Next from target protein sequences select those with domain architecure of interest (`config.yaml`).
1. **extracttrans** -- using `scripts/extractfasta`, extract finally selected sequences from the non-redundant set.
1. **signalp** -- using and independently installed [SignalP 5.0](https://services.healthtech.dtu.dk/services/SignalP-5.0) tool available as `signalp`, detect N-terminal signal sequences in the final set of protein sequences for each domain architecture of interest. If the tool is not available, the step generates an empty result file and is skipped.
1. **annotdom** -- using `scripts/annot.py`, prepare a GFF3 file with annotations for the finally selected protein sequences. Annotations are prepared based on `hmmsearch` filtered results and `signalp` results as well.

More detailed description on how the pipeline works you will find in comments in the `Snakefile`, `config.yaml` file and the script files.

### 5. Running the pipeline
Providing that you have `pfam-gen` environment properly set up and activated as well as at least the test data dowloaded, you can first look up the list of tasks to be done by running the following command from the pipeline directory, where `Snakefile` is located:
```bash
snakemake --dryrun --quiet
```
and then run the pipeline using as many cores as you wish:
```bash
snakemake --cores number_of_cores
```
Final data, non-redundant protein sequences and corresponding annotations for each domain architecture of interest, as well as HTML and PNG visualisatons of all domain architectures found, you will find in `output/final` directory.

