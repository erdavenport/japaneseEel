<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [RADseq of Japanese Eels](#radseq-of-japanese-eels)
  - [Overview](#overview)
- [Prep sequencing data:](#prep-sequencing-data)
  - [Step 1: Trim PE reads using fastx-toolkit](#step-1-trim-pe-reads-using-fastx-toolkit)
  - [Step 2: Combine reads from two sequencing batches together](#step-2-combine-reads-from-two-sequencing-batches-together)
- [Prep data for STACKS](#prep-data-for-stacks)
  - [Step 1: Make accessory files for STACKS](#step-1-make-accessory-files-for-stacks)
  - [Step 2: QC and demultiplex the reads](#step-2-qc-and-demultiplex-the-reads)
  - [Step 3: Download and prep reference genome](#step-3-download-and-prep-reference-genome)
  - [Step 4: Map reads to the nuclear and mitochondrial genome](#step-4-map-reads-to-the-nuclear-and-mitochondrial-genome)
- [Create SNP catalog using STACKS](#create-snp-catalog-using-stacks)
  - [Step 1: Create stacks using a reference genome (ref_map)](#step-1-create-stacks-using-a-reference-genome-ref_map)
  - [Step 2: Run corrections module (rxstacks)](#step-2-run-corrections-module-rxstacks)
  - [Step 3: Reassemble the catalog after running corrections module](#step-3-reassemble-the-catalog-after-running-corrections-module)
  - [Step 4: Run populations module on catalog to calculate stats](#step-4-run-populations-module-on-catalog-to-calculate-stats)
  - [Step 5: Filter out samples with low calls and loci with high depth, rerun populations](#step-5-filter-out-samples-with-low-calls-and-loci-with-high-depth-rerun-populations)
- [Run analyses on final data set](#run-analyses-on-final-data-set)
  - [Create basic stats plots](#create-basic-stats-plots)
  - [Create phylogenetic trees from population level Fst](#create-phylogenetic-trees-from-population-level-fst)
  - [Run admixture](#run-admixture)
  - [Run PCA](#run-pca)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

RADseq of Japanese Eels
=======

This repo contains all of the scripts and log files used to process the RADseq data for the japanese eel project led by Xiaoling Gong.
Note that all analyses start in the base directory. 

Overview
-------
	project
	|- README							# Description of analyses performed
	|
	|- data/							# Any data put into analyses - may be raw or processed (note: not version controlled currently due to size)
	|    |- fastq_files					# original fastq files from Xiaoling (all merged, don't use)
	|    |- fastq_files_updated			# fastq files by lane and run 
	|
	|- log/								# Contains electronic lab notebook files, titled by date
	|
	|- scripts/							# Contains all scripts used to run analyses
	|
	|- results/							# Contains all output from scripts (note: not version controlled currently due to size)
	|    |- 0_fastqc					# FastQC reports of the original sequence files
	+
# Prep sequencing data:

Xiaoling obtained four sequencing files from the company who sequenced the eels (stored on cbsufsrv5 at `data2/japaneseEel/data/fastq_files_updated/`):

* 32manli1_R1.fastq.gz
* 32manli2_R1.fastq.gz
* 36manli2_R1.fastq.gz
* 36manli2_R2.fastq.gz

Files that start with 32 are the single-end libraries. 
Files that start with 36 are the paired-end libraries.
manli1 is library mix 1.
manli2 is library mix 2.
R1 are all forward reads (only read for single-end sequencing).
R2 are the paired end reads (where applicable). 

## Step 1: Trim PE reads using fastx-toolkit

The single-end reads were only sequenced to 100bp, while the paired-end reads were sequenced to 125bp.

Trim off the last 25 bases of the PE reads from `36manli2_R1.fastq.gz` so that all the fragments are the same length (using fastx-toolkit) [from `cbsulogin`]:

```
# Make new directory for output files:
mkdir /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/36manli2_pe_trimmed/

# Unzip fastq.gz file:
gunzip /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/36manli2_R1.fastq.gz

# Submit script to trip reads to 100bp:
qsub scripts/submit_fastx_trim_pe_reads.sh /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/ 36manli2_R1.fastq /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/36manli2_pe_trimmed/

# Zip up fastq file:
gzip /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/36manli2_R1.fastq
```

## Step 2: Combine reads from two sequencing batches together 

Library 2 (li2) sequenced twice, once single-end and once paired-end.
Combining the forward reads two libraries together, using the trimmed PE reads [from `cbsufsrv5`]:

```
# Combine li2 files:
zcat data/fastq_files_updated/36manli2_pe_trimmed/36manli2_R1.fastq.gz data/fastq_files_updated/32manli2_R1.fastq.gz > ../data/fastq_files_updated/32_36_manli2_R1_merged.fastq

# zip up the merged file:
gzip data/fastq_files_updated/32_36_manli2_R1_merged.fastq
``` 

# Prep data for STACKS
## Step 1: Make accessory files for STACKS

First, make a barcode file for stacks [from `cbsulm06`]:

```
scripts/make_bcfiles_for_STACKS.R
```

Next, make a population map file for STACKS [from `cbsulm06`]:

```
scripts/make_population_map_file_for_STACKS.R
```

## Step 2: QC and demultiplex the reads

Finally, run `process_radtags` to QC the reads and demultiplex the samples [from `cbsulm06`]:

Process 32manli1 file:

```
/programs/stacks-1.48/bin/process_radtags \
	-f /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/32manli1_R1.fastq.gz \
	-i 'gzfastq' \
	-o /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out/ \
	-b data/sample_info/bcfile_manli1_for_stacks.txt \
	-c \
	-q \
	-D \
	-e 'kpnI' \
	--inline_null
```

Process merged manli2 file:

```
/programs/stacks-1.48/bin/process_radtags \
	-f /fs/cbsufsrv5/data2/japaneseEel/data/fastq_files_updated/32_36_manli2_R1_merged.fastq.gz \
	-i 'gzfastq' \
	-o /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out/ \
	-b data/sample_info/bcfile_manli2_for_stacks.txt \
	-c \
	-q \
	-D \
	-e 'kpnI' \
	--inline_null
```

## Step 3: Download and prep reference genome

Downloaded the japanese eel nuclear and mtDNA genomes from (<http://www.ncbi.nlm.nih.gov/genome/?term=txid7937[Organism:noexp]>):

Nuclear source:  
<ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000470695.1_japanese_eel_genome_v1_25_oct_2011_japonica_c401b400k25m200_sspacepremiumk3a02n24_extra.final.scaffolds/GCA_000470695.1_japanese_eel_genome_v1_25_oct_2011_japonica_c401b400k25m200_sspacepremiumk3a02n24_extra.final.scaffolds_genomic.fna.gz>

mtDNA source (send to fasta file and then gzipped):  
<http://www.ncbi.nlm.nih.gov/nuccore/595077858?report=fasta>

The mitochondrial genome is included in the nuclear source file, but we want to consider nuclear and mito DNA separately.

Use BBmap to remove that mito sequence [from `cbsulm06`]:

```
/programs/bbmap-35.66/filterbyname.sh \
	in=/fs/cbsufsrv5/data2/japaneseEel/genomes/a_japonica/GCA_000470695.1_japanese_eel_genome_v1_25_oct_2011_japonica_c401b400k25m200_sspacepremiumk3a02n24_extra.final.scaffolds_genomic.fna.gz \
	out=a_japonica_nuclear_genome.fasta \
	names=CM002536.1 \
	include=f
```

Unzip genome so Bowtie can build it:

```
gunzip /data/genomes/a_japonica/nuclear/a_japonica_nuclear_genome.fasta.gz 
```

Use Bowtie to make an index of the reference genome:

```
mkdir bowtie_indexes/ 
/programs/bowtie-1.1.2/bowtie-build a_japonica_nuclear_genome.fasta /data2/japaneseEel/data/genomes/a_japonica/nuclear/bowtie_indexes/
```

## Step 4: Map reads to the nuclear and mitochondrial genome

Unzip Stacks processed fastq files (in `data/STACKS_processed/0_process_radtags_out/` on `cbsufsrv5`)

```
for i in *.gz
do
	echo $i
	gunzip $i
done
```

From the login node, submit jobs to map reads [from `cbsulogin`]:

```
# Make directory for output
mkdir /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/

# Map reads
for i in /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out/*.fq
do
myarg=`basename $i`
qsub scripts/map_fastq_using_bowtie_2_mismatches_unique_only.sh \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out/ \
	$myarg \
	/fs/cbsufsrv5/data2/japaneseEel/data/genomes/a_japonica/nuclear/bowtie_indexes/ \
	a_japonica_nuclear \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/
done
```

Rezip up the fastq files (from `data/STACKS_processed/0_process_radtags_out/` on `cbsufsrv5`)

```	
for i in *.fq
do
	echo $i
	gzip $i
done
```

Pull unmapped reads from sam files [from `cbsulogin`]:

Make a directory for the unmapped read files:

```
mkdir /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/
```

Run scripts to pull reads:

```
for i in /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/*.sam
do
myarg=`basename $i`
qsub scripts/submit_convert_unmapped_sam_to_fastq.sh \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/ \
	$myarg \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/
done
```

Rename to remove the unmapped bit (from `data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/` on `cbsufsrv5`):

```
for i in *.fastq
do
	mv $i $(echo $i | sed s/_unmapped//)
done
```

Make directory for general processing information

```
mkdir results/1_general_info/
```

Double check that those files have the same number of reads as we expect from the sam files:

```
scripts/count_reads_not_mapped_to_nuclear_genome_fastq_files.sh
```

Mapping those reads to mtDNA [from `cbsulogin`]:

Make a directory for output:

```
mkdir /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/3_aligned_to_mitochondrial_genome/
```
Map the reads to mtDNA:

```
for i in /fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/*.fastq
do
myarg=`basename $i`
qsub /fs/cbsufsrv5/data2/japaneseEel/scripts/map_fastq_using_bowtie_2_mismatches_unique_only.sh \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/2_unmapped_reads_from_nuclear_genome_bowtie/ \
	$myarg \
	/fs/cbsufsrv5/data2/japaneseEel/data/genomes/a_japonica/mtDNA/bowtie_indexes/ \
	a_japonica_mtDNA \
	/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/3_aligned_to_mitochondrial_genome/
done
```

Remove the .fastq from the middle of the file name for those mtDNA mapped reads (in `/data2/japaneseEel/data/STACKS_processed/3_aligned_to_mitochondrial_genome/` on `cbsufsrv5`):

```
for i in *
do
	mv $i $(echo $i | sed s/.fastq//)
done
```


Count the number of reads mapping and unmapped after bowtie:

```
# Make directory for results:
mkdir results/2_processing_info/

scripts/summarize_mapped_reads_after_bowtie_nuclear.R
```

Summarize the number of reads mapped to mtDNA:

```
scripts/summarize_mapped_reads_after_bowtie_mtDNA.R
```


# Create SNP catalog using STACKS

## Step 1: Create stacks using a reference genome (ref_map)

## Step 2: Run corrections module (rxstacks)

## Step 3: Reassemble the catalog after running corrections module

## Step 4: Run populations module on catalog to calculate stats

## Step 5: Filter out samples with low calls and loci with high depth, rerun populations

# Run analyses on final data set

## Create basic stats plots
## Create phylogenetic trees from population level Fst
## Run admixture
## Run PCA





