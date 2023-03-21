#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Extracts values from a given column from a TSV file and writes them to stdout.
# Arguments:
# --input  : input TSV file with named columns
# --column : name of the column values are to be extracted from
# USAGE:
# ./extractpfm.sh --input INPUT_TSV --column COL_NAME

set -euo pipefail
# Create a list of expected command line arguments.
EXPARGS=("input" "column")
USAGE="USAGE: ./extractpfm.sh --input INPUT_TSV --column COL_NAME"

# Parse command line arguments. Every argument value starting with double dash
# treat as an argument name (strip dashes) and the next value, as long as does
# not start with a double dash, as the argument value. Put these value pairs to
# an associative array args.
declare -A args
raw_args=("${@}")
for i in "${!raw_args[@]}"; do
    if [[ "${raw_args[$i]:0:2}" == "--" && "${raw_args[$i+1]:0:2}" != "--" ]]; then
        args["${raw_args[$i]:2}"]="${raw_args[$i+1]}"
    fi
done
# Check wheter ars array contains values for all required arguments (keys).
for key in "${EXPARGS[@]}"; do
    if [[ "${args[$key]}" == "" ]]; then
        echo $USAGE >&2
        exit 1
    fi
done

# Iterate over input file lines. In the first line find the column of interest to
# extract values of that column from subsequent lines. Write the values to stdout.
index=-1
while IFS=$'\t' read -a line; do
   if (( index == -1 )); then
       for i in "${!line[@]}"; do
           if [[ "${line[$i]}" == "${args['column']}" ]]; then
               index=$i
               break
           fi
       done
       if (( index == -1 )); then
           echo "Cannot find ${args['column']} column" >&2
           exit 1
       fi
   else
       echo "${line[$index]}"
   fi
done < "${args['input']}"

