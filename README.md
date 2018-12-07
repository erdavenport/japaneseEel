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
## Step 1: Trim PE reads using fastx-toolkit

## Step 2: Combine reads from two sequencing batches together 

# Prep data for STACKS
## Step 1: Make accessory files for STACKS

## Step 2: QC and demultiplex the reads

## Step 3: Download and prep reference genome

## Step 4: Map reads to the nuclear and mitochondrial genome

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





