#!bin/bash

# User input:
accession_ID=$1
input_file=$2
alignmnet_file=$3
output_name=$4
model=$5                        # Select auto if don't know which model to select

# Pathways:
serratus="../Serratus"                          # Main Serratus dir.
acc_dir="${serratus}/${accession_ID}"           # Serratus subdir for the accession ID.
tree_dir="${acc_dir}/Tree"                                              # Tree subdir for output.

# Check if correct accession ID has been entered (its directory should exist):
if [ ! -d "${acc_dir}" ]; then
        echo "Incorret accession ID: ${accession_ID}"
        exit 1
fi

# Creating directories, if they do not exist:
mkdir -p ${tree_dir}


# Script:
mafft --thread 8 --add ${input_file} --keeplength ${alignmnet_file} > ${tree_dir}/${output_name}
iqtree2 -nt 8 -s ${tree_dir}/${output_name} -m ${model} -bb 1000 -czb
