# Select Accession ID you want to work on, then select files containing SRA IDs that
# should be intercepted (if more than 1 file is selected). Finally it opens interecepted
# SRA.fa consensus files and filters empty sequences and the ones that do not pass the
# specified filter.

# Import libraries:
import argparse
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Setup argparse:
parser = argparse.ArgumentParser(
    prog="GenBankIDs inteceptor and filter",
    description="Intercepts GenBankIDs, and creates filtere fasta file.")
parser.add_argument("accession_ID", help="Input GenBank accession ID")
parser.add_argument("-files", nargs='+', required=True, help='Input one or more files with SRA IDs for interception')
parser.add_argument('-output', required=True, help='Specify output prefix name')
parser.add_argument('-perc', default=50, type=int, choices=range(0, 101), help='Specify value% to filter sequences with % of Ns')
args = parser.parse_args()


# File pathways:
accession_ID = args.accession_ID
sra_files = args.files
ratio = args.perc/100
serratus = "../Serratus"
acc_dir = f"{serratus}/{accession_ID}"
data = f"{acc_dir}/Data"
neighbours = f"{acc_dir}/Neighbours"
consensus = f"{acc_dir}/Consensus"


# Check if directory exists, if not either create it or raise error:
if not os.path.exists(acc_dir):
    raise SystemExit(f"Incorrect accession ID: {accession_ID}")
if not os.path.exists(data):
    raise SystemExit(f"Data directory does not exist for {accession_ID}")
if not os.path.exists(consensus):
    raise SystemExit(f"Consensus directory does not exist for {accession_ID}")


# Script:
SRA_IDs =[]

# Intercept files if needed, and obtain SRA IDs:
if len(sra_files) == 1:
    filename = sra_files[0]
    file_location = f'{acc_dir}/{filename}'
    try:
        with open(file_location, 'r') as file:
            content = file.read()
            SRA_IDs = content.split()
    except FileNotFoundError as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)
elif len(sra_files) > 1:
    all_SRAs = []
    for filename in sra_files:
        file_location = f'{acc_dir}/{filename}'
        try:
            with open(file_location, 'r') as file:
                content = file.read()
                content = content.split()
                all_SRAs.append(content)
        except FileNotFoundError as e:
            print(f"Unexpected error: {e}")
            sys.exit(1)
    SRA_IDs = set.intersection(*map(set, all_SRAs))
    SRA_IDs = list(SRA_IDs)

print(f'\nIn total {len(SRA_IDs)} fasta sequences will be parsed through the filter.\n')

# Now open each corresponding consensus fasta file and filter zero ones and ones that do not pass the % filter:
fasta_sequences = {}

for SRA_ID in SRA_IDs:
    fasta_file = f'{consensus}/{SRA_ID}.fa'
    try:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if len(record.seq) != 0:
                sequence = str(record.seq)
                n_count = sequence.count('N')
                n_ratio = n_count / len(sequence)
                if n_ratio < ratio:
                    fasta_sequences[record.id] = str(record.seq)
    except FileNotFoundError as e:
            print(f'Unexpected error: {e}')
            print('Continue...')
            continue

# Save filtered sequences into one fasta file:
seq_records = []
for header, sequence in fasta_sequences.items():
    seq_record = SeqRecord(Seq(sequence), id=header, description="")
    seq_records.append(seq_record)

fasta_output = f'{acc_dir}/{args.output}.fasta'
with open(fasta_output, "w") as output_handle:
    SeqIO.write(seq_records, output_handle, "fasta")

print(f'A total of {len(fasta_sequences)} fasta sequences passed the {ratio} N filter')
print(f'Results saved at {fasta_output}\n')
