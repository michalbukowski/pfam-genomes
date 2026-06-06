#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Prepares complete test input data for pfam-genomes pipeline by fetching
# test genomes listed in input/test_genomes.txt from NCBI to ./genomes directory
# and Pfam database to ./input directory.

set -euo pipefail

# Path to the directory and a filename of Pfam database od HMM profiles.
PFAM_PATH="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FILE="Pfam-A.hmm.gz"

# File with assembly accession numbers of test genomes.
GENOMES_FILE="input/test_genomes.txt"

# Check the path for saving Pfam database to. If fails, exit with code 1.
if [[ ! -d "input" ]]; then
    echo "Directory input/ does not exist in the current location." \
         "Make sure you downloaded a complete pfam-genomes repository." >&2
    exit 1
fi
if [[ -e "input/${PFAM_FILE}" ]]; then
    echo "The path exists: input/${PFAM_FILE}." \
         "Rename of remove the existing file or directory and start again." >&2
    exit 1
fi

# Check the path for saving genomes to. If fails, exit with code 1.
if [[ -e "genomes" ]]; then
    echo "Cannot create genomes/ output directory. The path exists." \
         "Rename of remove the directory and start again." >&2
    exit 1
fi

# Read genome accessions from the input file. If fails, exit with code 1.
if [[ ! -f "${GENOMES_FILE}" ]]; then
    echo "Genome accessions file not found: ${GENOMES_FILE}" >&2
    exit 1
fi
mapfile -t ACCESSIONS < "${GENOMES_FILE}"
count="${#ACCESSIONS[@]}"

# Warn the user about the size of dowloaded and decomressed data. On rejection,
# exit with code 0.
echo "This script will fetch ${count} test genomes to genomes/ directory and" \
     "Pfam database to input/ directory (~280 MB, ~1.6 GB when decompressed)"
answer=""
while [[ "${answer}" != "yes" && "${answer}" != "no" ]]; do
    echo "Do you want to proceed? (yes/no)"
    read -r answer
done
if [[ "${answer}" == "no" ]]; then
    echo "Exiting..."
    exit 0
fi

# Download Pfam database. If fails, exit with code 1.
echo "Fetching Pfam database..."
wget -q --show-progress -P "input" "${PFAM_PATH}/${PFAM_FILE}"
if (( $? == 0 )); then
    echo "File ${PFAM_FILE} successfully fetched" >&2
else
    echo "Cannot fetch file from: ${PFAM_PATH}/${PFAM_FILE}" >&2
    exit 1
fi
echo "Done"

# Download test genomes using NCBI datasets CLI. If fails, exit with code 1.
echo "Fetching ${count} test genomes..."
mkdir -p "genomes"
tmp_zip="genomes_ncbi_download.zip"
tmp_dir="genomes_ncbi_tmp"

datasets download genome accession \
    --inputfile "${GENOMES_FILE}"  \
    --include   genome             \
    --filename  "${tmp_zip}"

unzip -q "${tmp_zip}" -d "${tmp_dir}"

i=0
for acc in "${ACCESSIONS[@]}"; do
    i=$((i + 1))
    fna_file=$(find "${tmp_dir}/ncbi_dataset/data/${acc}" -name "*.fna" 2>/dev/null | head -1)
    if [[ -z "${fna_file}" ]]; then
        echo "Cannot find genome file for: ${acc}" >&2
        rm -rf "${tmp_dir}" "${tmp_zip}"
        exit 1
    fi
    gzip -c "${fna_file}" > "genomes/${acc}_genomic.fna.gz"
    echo "Genome ${i}/${count} ${acc} successfully fetched" >&2
done

rm -rf "${tmp_dir}" "${tmp_zip}"
echo "Done"
