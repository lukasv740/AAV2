#!/bin/bash

# User inputs:
query=$1			# Input your query fasta file for blastn.
contigs=$2			# Input your contig fasta file for database creation and selection.


# Pathways (Create if does not exist):
unxplore="../UnXplore"
refs="${unxplore}/Refs"
databases="${unxplore}/Databases"
blastn_results="${unxplore}/Blastn_results"

mkdir -p ${unxplore}
mkdir -p ${refs}
mkdir -p ${databases}
mkdir -p ${blastn_results}


# Test if inputs exist:
if [ ! -f "${query}" ]; then
	echo "${query} does not exist."
	exit
fi
if [ ! -f "${contigs}" ]; then
	echo "${contigs} does not exist."
	exit
fi


# General script:
output="$(basename ${contigs} | cut -d '.' -f1)"	# Gets the prefix name of contigs file.
query_name="$(basename ${query} | cut -d '.' -f1)"	# Gets the prefix name of query file.
database_output="${databases}/${output}"		# Subdir for database location.
blastn_results_output="${blastn_results}/${output}"	# Subdir for blastn results.

# Check if database exists, if not, create one:
if [ ! -d "${database_output}" ]; then
	makeblastdb -dbtype nucl -in ${contigs} -out ${database_output}/${output} -title ${output}
	mkdir -p ${blastn_results_output}
else
	echo -e "\nDatabase ${output} alredy exists.\n"
	mkdir -p ${blastn_results_output}
fi


# Perform blastn for the query:
echo -e "Performing BLASTN.\n"
blastn -db ${database_output}/${output} -query ${query} -outfmt 6 -num_threads 8 -out ${blastn_results_output}/${query_name}.txt

# Display results:
count="$(cat ${blastn_results_output}/${query_name}.txt | wc -l)"
echo -e "Total hits ${count}."
echo -e "Results saved in ${blastn_results_output}/${query_name}.txt \n"
