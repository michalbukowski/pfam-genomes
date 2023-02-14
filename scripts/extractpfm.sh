#!/bin/bash
EXPARGS=("domdata" "column")
USAGE="USAGE: extractpfm.sh --domdata INPUT_TSV --column PFAMACC_COLNAME"

declare -A args
raw_args=("${@}")
for i in "${!raw_args[@]}"; do
    if [[ "${raw_args[$i]:0:2}" == "--" && "${raw_args[$i+1]:0:2}" != "--" ]]; then
        args["${raw_args[$i]:2}"]="${raw_args[$i+1]}"
    fi
done
for key in "${EXPARGS[@]}"; do
    if [[ "${args[$key]}" == "" ]]; then
        echo $USAGE >&2
        exit 1
    fi
done

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
done < "${args['domdata']}"

