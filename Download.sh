#!/bin/bash

# User input:
accession_ID=$1									# Input Serratus accession ID.


# Check if accession number exists on Serratus database:
URL="https://api.serratus.io/matches/nucleotide/download?scoreMin=50&scoreMax=100&identityMin=75&identityMax=100&sequence=${accession_ID}"
if wget --spider "${URL}" 2>/dev/null; then
    echo -e "\nAccession ID ${accession_ID} is correct!\n"
else
    echo -e "\nError! Incorrect accession ID ${accession_ID} provided! \n"
    exit
fi


# Pathways:
serratus="../Serratus"							# Main Serratus dir.
acc_dir="${serratus}/${accession_ID}"			# Serratus subdir for the accession accession ID.
data="${acc_dir}/Data"							# Data folder for .bam, bam.ai, .stat and .summary files.


# Creating directories, if they do not exist:
mkdir -p ${serratus}
mkdir -p ${acc_dir}
mkidr -p ${data}


# General script:
# Gets a summary file for the input accesion number and cuts the first column to get SRA ids:
echo -e "Downloading data\n"
wget --output-document=${acc_dir}/${accession_ID}_matches.csv "https://api.serratus.io/matches/nucleotide/download?scoreMin=50&scoreMax=100&identityMin=75&identityMax=100&sequence=${accession_ID}"	# Downloads csv summary file of specified accession ID.
cut -d "," -f 1 ${acc_dir}/${accession_ID}.csv | tail +2 > ${acc_dir}/${accession_ID}_SRAs.txt	# Cuts the first column of generated csv file without header to get SRA IDS.

# Download BAM files, index them, get sumamry files, get the idxstat file and sort it:
while read -r line
do
	wget --output-document=${data}/${line}.bam "https://s3.amazonaws.com/lovelywater2/bam/${line}.bam"					# Get bam files.
	wget --output-document=${data}/${line}.summary "https://s3.amazonaws.com/lovelywater2/summary2/${line}.summary"		# Get summary file for each bam file.
	samtools index ${data}/${line}.bam															# Index bam files.
	samtools idxstats ${data}/${line}.bam | sort -nr -k3 > ${data}/${line}.sorted.stat			# Creates a sorted stat file that shows number of hits.
	stats_output=$(grep "${accession_ID}" ${data}/${line}.sorted.stat)							# Extracts number of hits for the accession number from stat file.
	echo -e "${line}\t${stats_output}" >> ${acc_dir}/${accession_ID}.stat						# Add stats info into the stat file i.e. /query /ref /lenght /#matched /#unmatched

done < ${acc_dir}/${accession_ID}_SRAs.txt

# Sort generated stat file by the number of hits in descending order:
sort ${acc_dir}/${accession_ID}.stat -nr -k4 > ${acc_dir}/${accession_ID}.sorted.stat
rm ${acc_dir}/${accession_ID}.stat

count="$(cat ${acc_dir}/${accession_ID}.sorted.stat | wc -l)"
echo -e "\nA total of ${count} have been downloaded and saved in ${data}\n"
