#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Prepares complete test input data for pfam-genomes pipeline by fetching
# 20 staphylococcal test genomes from NCBI GenBank to ./genomes directory and
# Pfam database to ./input directory.

set -euo pipefail

# Path to the directory and a filename of Pfam database od HMM profiles.
PFAM_PATH="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FILE="Pfam-A.hmm.gz"

# Root FTP location and sublocations for 20 staphylococcal test genomes.
FTP_ROOT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"
FTP_PATHS=(
    "000/276/045/GCA_000276045.1_ASM27604v1"
    "000/313/005/GCA_000313005.1_Staaur_20"
    "000/363/645/GCA_000363645.1_Stap_aure_M1521_V1"
    "000/364/105/GCA_000364105.1_Stap_aure_M1228_V1"
    "000/382/965/GCA_000382965.1_ASM38296v1"
    "000/529/315/GCA_000529315.1_Stap_aure_M0702_V1"
    "000/530/215/GCA_000530215.1_S0460D1"
    "000/538/115/GCA_000538115.1_Stap_aure_SJOS6072_V1"
    "000/538/235/GCA_000538235.1_Stap_aure_HBHO6042_V1"
    "000/541/615/GCA_000541615.1_Stap_aure_OCMM6075_V1"
    "000/546/925/GCA_000546925.1_Stap_aure_M1067_V1"
    "000/557/105/GCA_000557105.1_Stap_aure_M0533_V1"
    "000/561/525/GCA_000561525.1_Stap_aure_H33602_V1"
    "000/570/215/GCA_000570215.1_Stap_aure_F57021_V1"
    "000/570/655/GCA_000570655.1_Stap_aure_T37840_V1"
    "000/594/685/GCA_000594685.1_Stap_aure_T13125_V1"
    "000/644/195/GCA_000644195.1_Stap_aure_VET0388R_V1"
    "000/683/255/GCA_000683255.1_Stap_aure_VET0489R_V1"
    "001/019/375/GCA_001019375.2_ASM101937v2"
    "001/045/625/GCA_001045625.1_ASM104562v1"
)
# Name suffic of files contatining nucleotide genomic sequences.
SUFFIX="_genomic.fna.gz"

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

# Warn the user about the size of dowloaded and decomressed data. On rejection,
# exit with code 0.
echo "This script will fetch 20 test genomes to genomes/ directory (~17 MB) and" \
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

# Decompress Pfam database. If fails, exit with code 1.
echo "Decompressing Pfam database..."
gunzip "input/${PFAM_FILE}"
if (( $? != 0 )); then
    echo "Cannot decompress file: input/${PFAM_FILE}" >&2
    exit 1
fi
echo "Done"

# Download test genomes. If fails at any point, exit with code 1.
count="${#FTP_PATHS[@]}"
echo "Fetching ${count} test genomes.."
mkdir -p "genomes"
for i in "${!FTP_PATHS[@]}"; do
    ftp_path="${FTP_PATHS[$i]}"
    asmname="${ftp_path##*/}"
    if [[ ! "${asmname}" =~ (GCA_[0-9]+\.[0-9]+) ]]; then
        echo "Cannot extract assembly accession from: ${ftp_path}" >&2
        exit 1
    fi
    asmacc="${BASH_REMATCH[1]}"
    inpath="${FTP_ROOT}/${ftp_path}/${asmname}${SUFFIX}"
    outpath="genomes/${asmacc}${SUFFIX}"
    wget -q --show-progress -O "genomes/${asmacc}${SUFFIX}" "${inpath}"
    if (( $? == 0 )); then
        echo "Genome $((i+1))/${count} ${asmacc} successfully fetched" >&2
    else
        echo "Cannot fetch genome from: ${inpath}" >&2
        exit 1
    fi
done
echo "Done"

