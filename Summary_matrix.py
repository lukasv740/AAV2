# Takes each bam summary file, extracts family names and number of reads, and creates a count matrix.
# Additionally there are output name and filter (number of reads) options available.

# Import libraries:
import os
import re
import pandas as pd
import argparse


# Setting up argparse:
parser = argparse.ArgumentParser(
    prog="Summary matrix generator",
    description="Creates a summary matrix (csv file) from summary files, showing alignment read count for each family")
parser.add_argument("accession_ID", help="Input GenBank accession ID")
parser.add_argument("-output", default="Summary_matrix", help="Choose your output prefix name")
parser.add_argument("-filter", default=100, help="Choose your filter value for reads (min number of reads for each family)")
parser.add_argument("-sra", default=None, help="Choose your file with SRA IDs (relative path from ~/Serratus/NC_xxxx/")
args = parser.parse_args()


# File Pathways:
accession_ID = args.accession_ID
serratus = "../Serratus"
acc_dir = f"{serratus}/{accession_ID}"
data = f"{acc_dir}/Data"


# Check if directories exist and not empty:
if not os.path.exists(acc_dir):
    raise SystemExit(f"Incorrect accession ID: {accession_ID}")
if not os.path.exists(data):
    raise SystemExit(f"Data directory does not exist for {accession_ID}")
if not os.listdir(data):
    raise SystemExit(f"Data directory is empty for {accession_ID}")

# Check if specified SRA file exists:
if args.sra == None:
    sra_file = f'{acc_dir}/{accession_ID}_SRAs.txt'
    if not os.path.exists(sra_file):
        raise SystemExit(f"File {sra_file} does not exist")
else:
    sra_file = f'{acc_dir}/{args.sra}'
    if not os.path.exists(sra_file):
        raise SystemExit(f"File {sra_file} does not exist")


# Script:
# Get each individual file name with each corresponding family name and its aln length:
matrix_dict = {}

# Parse through each file:
# First opens SRA file and parses through SRA IDs:
with open(sra_file, 'r') as file:
    for sample_name in file:
        sample_name = sample_name.strip()
        file_name = f'{sample_name}.summary'

        # Opens each individual summary file, and parses through it, looking for family name and its aln length:
        try:
            with open(f'{data}/{file_name}', 'r') as summary:
                data_dict = {}
                for row in summary:
                    family_name = re.search(r'fam=([^;]+)', row)
                    aln_length = re.search(r'aln=([^;]+)', row)

                    # Append family name to its aln length:
                    if family_name != None and aln_length != None:
                        if int(aln_length.group(1)) >= int(args.filter):
                            data_dict[family_name.group(1)] = aln_length.group(1)

                # Append data to each file in a dictionary:
                matrix_dict[sample_name] = data_dict
        except FileNotFoundError as e:
            print(f"Unexpected error: {e}")
            print("Continue...")
            continue


# Creating an alphabetically sorted summary matrix:
df = pd.DataFrame(matrix_dict).T
df.fillna(0, inplace=True)
df = df.sort_index(axis=1)

# Export it as csv:
summary_matrix_output_path = f'{acc_dir}/{args.output}.csv'
df.to_csv(summary_matrix_output_path)

