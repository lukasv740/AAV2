# Overview
The overall analysis of this study consists of two parts. Analysis of the UnXplore dataset was performed using one custom made bash script, whereas Serratus dataset was analysed using three bash, three python and one R script. 

# UnXplore
The `UnXplore.sh` script requires two positional inputs from the user - a FASTA file with genomic sequency query, and a FASTA file with assembled contigs for a local database creation. It works by first generating a local nucleotide BLAST database using 'makeblastdb' command from NCBI BLAST suite of programs. Once the database is created, the query sequence is seached using 'bastn' command. The results displayed in a text file in a default BLASTN tabular output format 6. Both results and generated databases are saved in `../UnXplore` directory under specific subdirectories. The script creates all required directories automatically.
To run the script the user just simply inputs query FASTA sequence as a first argument, and contig FASTA file as a second:
```
bash UnXplore.sh <query> <database>
```
The required aligned contig datasets are public available in Modha et al., 2022 study. Due to the sheer size (300b long contig file requires over 80GB of space), the dataset was not included.

The paper: Modha, S., Robertson, D.L., Hughes, J. and Orton, R.J. (2022). Quantifying and Cataloguing Unknown Sequences within Human Microbiomes. mSystems. doi:https://doi.org/10.1128/msystems.01468-21.




# Serratus
The overall workflow for analysing Serratus database is outlined in Figure 3. It consists of downloading a Serratus summary file for an entered GenBank ID, followed by an extraction of the SRA sample accession IDs which are used to download corresponding SRA summary and BAM files. These files contain reads that were mapped to the pangenome of viral sequences as well as a text overview of number of reads mapping to different viral taxonomic families. These SRA summary files are used to create a family count summary matrix, which are subsequently analysed and visualised. Other steps include neighbourhood filtering, generation of consensus sequences, low quality sequences filtering, their alignment to study sequences and phylogenetic tree analysis. Documentation how to run each script, how it works, and recommended order is described below:

![](/Plots/Master_project.jpeg)
