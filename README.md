# UnXplore
The `UnXplore.sh` script requires two positional inputs from the user - a FASTA file with genomic sequency query, and a FASTA file with assembled contigs for a local database creation. It works by first generating a local nucleotide BLAST database using 'makeblastdb' command from NCBI BLAST suite of programs. Once the database is created, the query sequence is seached using 'bastn' command. The results displayed in a text file in a default BLASTN tabular output format 6. Both results and generated databases are saved in `../UnXplore` directory under specific subdirectories. The script creates all required directories automatically.
To run the script the user just simply inputs query FASTA sequence as a first argument, and contig FASTA file as a second:
```
bash UnXplore.sh <query> <database>
```




# Serratus
Here are the steps for analysing Serratuts database for AAV2 presence.
