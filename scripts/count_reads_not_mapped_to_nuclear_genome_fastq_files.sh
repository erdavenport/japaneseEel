#!/bin/bash

# This script will count the number of reads for each of the files containing reads not mapped to the nuclear genome:

cd data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/
for i in *.fastq
do
echo -en "$i \t" >> ../../../results/1_general_info/reads_per_sample_not_mapped_to_nuclear_genome_fastq.txt
LINES="$(cat $i | wc -l)"
echo $((LINES/4)) >> ../../../results/1_general_info/reads_per_sample_not_mapped_to_nuclear_genome_fastq.txt
done



echo "DONE!"


