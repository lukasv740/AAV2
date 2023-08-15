#!/bin/bash

# User input:
accession_ID=$1                         # Input Serratus accession ID.


# Pathways:
serratus="../Serratus"                          # Main Serratus dir.
acc_dir="${serratus}/${accession_ID}"           # Serratus subdir for the accession ID,
data="${acc_dir}/Data"                          # Data folder for .bam, .bam.ai, .stat and .summary files.
consensus="${acc_dir}/Consensus"                # Consensus folder for the consensus sequences for each SRA.

# Check if correct accession ID has been entered (its directory should exist):
if [ ! -d "${acc_dir}" ]; then
        echo "Incorret accession ID: ${accession_ID}"
        exit 1
fi

# Creating directories, if they do not exist:
mkdir -p ${serratus}
mkdir -p ${acc_dir}
mkdir -p ${data}
mkdir -p ${consensus}



# General script:
while read -r line
do
        samtools mpileup -a -A -d 0 -Q 0 -r ${accession_ID} ${data}/${line}.bam | ivar consensus -p ${line} -t 0.4
        mv ${line}.fa ${consensus}
        mv ${line}.qual.txt ${consensus}
        cat ${consensus}/${line}.fa | paste - - | cut -f2 | tr -d '\n' | grep -o -e A -e T -e C -e G -e N -e - | sort | uniq -c > ${consensus}/${line}.stat
done < ${acc_dir}/${accession_ID}_SRAs.txt

echo -e "\nConsensus sequences for ${accession_ID} have been generated and saved at ${consensus}\n"