# The program works by specifying on which GenBank ID we will be working on.
# Then using -query option user specifies what he wants to look for (i.e. AAV2 virus in humand adenovirus F)
# In addition to that user can input a txt file with neighbours (GenBank) IDs for filtering.
# Finally, user can specify output directory name, for more efficient control

# Libraries:
import argparse
import os
import csv
from operator import itemgetter
import sys


# Setting up argparse:
parser = argparse.ArgumentParser(
    prog="Neighbour parser",
    description="Creates a sorted stat file using a query ID + a list of desired neighbours.")
parser.add_argument("accession_ID", help="Input GenBank accession ID")
parser.add_argument("-query", required=True, help="Input your GenBank ID query")
parser.add_argument("-neighbours_file", help="Input a txt file with neighbour list")
parser.add_argument("-output", required=True, help="Specify output prefix name, if test is run multiple times")
args = parser.parse_args()


# File pathways:
accession_ID = args.accession_ID
output = args.output
serratus = "../Serratus"
acc_dir = f"{serratus}/{accession_ID}"
data = f"{acc_dir}/Data"
neighbours = f"{acc_dir}/Neighbours"

# Check if directory exists, if not either create it or raise error:
if not os.path.exists(acc_dir):
    raise SystemExit(f"Incorrect accession ID: {accession_ID}")
if not os.path.exists(data):
    raise SystemExit(f"Data directory does not exist for {accession_ID}")
if not os.path.exists(neighbours):
    os.makedirs(neighbours)


# Script:
# The program works by selecting a NC folder you want to work on, then specifying query ID
# and optionally a list of neighbour IDs. It parses through each sorted bam file summary,
# and filters entries for the query result. Output statistics is forwarded into Neighbours directory,
# generating specified output name .pos.stat, .neg.stat and .general.stat files, as well as
# generating output.txt file with SRA IDs that passed the filter.

# Open file with neighbours GenBank IDs, if provided:
if args.neighbours_file != None:
    try:
        with open(args.neighbours_file, "r") as file:
            content = file.read()
            neighbours_list = content.split()
    except FileNotFoundError as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)
else:
    neighbours_list = []


pos_list = []       # Positive list (samples that passed the filter)
neg_list = []       # Negative list (samples that did not pass the filter)
pos_count = 0
neg_count = 0

# Open file with GenBank IDs:
with open(f"{acc_dir}/{accession_ID}_SRAs.txt", "r") as file_ids:
    content = file_ids.read()
    content = content.split()

    # Open each .sorted.stat file for each SRA ID:
    for SRA_ID in content:
        try:
            with open(f"{data}/{SRA_ID}.sorted.stat", "r") as stat_file:
                stat_reader = csv.reader(stat_file, delimiter="\t")
                for row in stat_reader:
                    # If no neighbours provided, select first row:
                    if neighbours_list == []:
                        stat_row = row
                        stat_row[2] = int(stat_row[2])
                        break

                    # If neighbours provided, select first row that matches one of the neighbours:
                    elif row[0] in neighbours_list or row[0] in args.query:
                        stat_row = row
                        stat_row[2] = int(stat_row[2])
                        break

                # If this row is the same as query, put it into pos filtered list, else into neg:
                # It works because .sorted.stat is sorted in descending order, so it will hit the first highest number:
                if stat_row[0] == args.query:
                    pos_list.append([SRA_ID] + stat_row)
                    pos_count += 1
                else:
                    neg_list.append([SRA_ID] + stat_row)
                    neg_count += 1

        except FileNotFoundError as e:
            print(f"Unexpected error: {e}")
            print("Continue...")
            continue


# Sort lists by the number of hits in descending order:
pos_list = sorted(pos_list, key=itemgetter(3), reverse=True)
neg_list = sorted(neg_list, key=itemgetter(3), reverse=True)


# Write list into corresponding pos and neg files:
# For positive samples:
# Stat file:
with open(f"{neighbours}/{output}_pos.stat", "w", newline="") as pos_out:
    pos_writer = csv.writer(pos_out, delimiter="\t")
    for item in pos_list:
        pos_writer.writerow(item)
# IDs:
with open(f"{neighbours}/{output}_pos_IDs.txt", "w", newline="\n") as pos_id_out:
    for item in pos_list:
        pos_id_out.write(f"{item[0]}\n")



# For negative samples:
# Stat file:
with open(f"{neighbours}/{output}_neg.stat", "w", newline="") as neg_out:
    neg_writer = csv.writer(neg_out, delimiter="\t")
    for item in neg_list:
        neg_writer.writerow(item)
# IDs:
with open(f"{neighbours}/{output}_neg_IDs.txt", "w", newline="\n") as neg_id_out:
    for item in neg_list:
        neg_id_out.write(f"{item[0]}\n")


# For general stats:
with open(f"{neighbours}/{output}_general.stat", "w", newline="\n") as gen_out:
    gen_out.write(f"The query for {accession_ID} was {args.query}\n")
    gen_out.write(f"A total of {pos_count} samples passed the filter.\n")
    gen_out.write(f"A total of {neg_count} samples did not pass the filter.\n\n")
    gen_out.write("Neighbours used:\n")
    for id in neighbours_list:
        gen_out.write(f"{id}\n")

print(f"\nFiles for neighbour filtering have been created and saved at {neighbours}\n")
