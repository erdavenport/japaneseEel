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
  - [Step 6: Run populations module on "meta-populations"](#step-6-run-populations-module-on-meta-populations)
    - [Final data sets](#final-data-sets)
- [Run analyses on final data set](#run-analyses-on-final-data-set)
  - [Create basic stats plots](#create-basic-stats-plots)
  - [Pull data for Table 2](#pull-data-for-table-2)
  - [Create phylogenetic trees from population level Fst](#create-phylogenetic-trees-from-population-level-fst)
  - [Run admixture](#run-admixture)
    - [Admixture - all individuals](#admixture---all-individuals)
    - [Admixture - no outgroup](#admixture---no-outgroup)
  - [Run PCA](#run-pca)
  - [Calculate significant differences in diversity statistics using ANOVA](#calculate-significant-differences-in-diversity-statistics-using-anova)
  - [Marker vs. Sampling size simulations](#marker-vs-sampling-size-simulations)
  - [Calculate stats for manuscript](#calculate-stats-for-manuscript)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

RADseq of Japanese Eels
=======

This repo contains all of the scripts and log files used to process the RADseq data for the paper ["Lack of spatial and temporal genetic structure of Japanese eel (_Anguilla japonica_) populations"](https://link.springer.com/article/10.1007/s10592-019-01146-8), led by Xiaoling Gong.
A [PDF](https://link.springer.com/article/10.1007/s10592-019-01146-8) of the article and the [raw data](https://doi.org/10.5061/dryad.4v286) are available. 
Note that all analyses start in the base directory. 

Overview
-------
	project
	|- README				# Description of analyses performed
	|
	|- data/				# Any data put into analyses - may be raw or processed (note: not version controlled currently due to size)
	|    |- fastq_files_updated		# Contains raw fastq files
	|    |- genomes				# Contains downloaded A. japonica reference genome
	|    |- sample_info			# contains sample information files 
	|    |- STACKS_processed		# Contains all of the files, processed by STACKS (not version controlled due to size)
	|
	|- log/					# Contains electronic lab notebook files, titled by date
	|
	|- scripts/				# Contains all scripts used to run analyses
	|
	|- results/				# Contains all output from scripts (note: not version controlled currently due to size)
	|    |- 1_general_info			# Contains sample and read info
	|    |- 2_processing_info		# Contains depth and alignment info
	|    |- 3_optimizing_depth		# Plots and stats from STACKS depth optimization
	|    |- 4_Fst				# Fst tables and plots for each depth
	|    |- 5_admixture			# Admixture results and plots for each depth
	|    |- 6_PCA				# PCA plots for each depth
	|    |- 7_ANOVA				# ANOVA p-values m3 only
	+

A table listing the Figure/Table in the final paper, script, and generated file (with location) is listed below:

| Item                   | file                                                                                            | script                         |
|------------------------|-------------------------------------------------------------------------------------------------|--------------------------------|
| Figure 1a              | results/1\_general\_info/map\_sample\_locations.pdf                                                 | make\_map\_for\_eelseq\_project.R  |
| Figure 1b              | results/1\_general\_info/timeline\_sample\_collection.pdf                                           | make\_map\_for\_eelseq\_project.R  |
| Figure 2               | results/3\_optimizing\_depth/barplot\_variant\_sites\_per\_pop\_across\_min\_stack\_depths.pdf            | plot\_loci\_per\_population.R     |
| Figure 3               | results/3\_optimizing\_depth/barplot\_diversity\_by\_population\_m3.pdf                               | plot\_loci\_per\_population.R     |
| Figure 4a              | results/3\_optimizing\_depth/m3\_plot\_time\_distance\_by\_genetic\_distance\_yangzte.png                | plot\_Fst\_over\_time\_and\_space.R |
| Figure 4b              | results/3\_optimizing\_depth/m3\_plot\_geographic\_distance\_by\_genetic\_distance\_no\_out\_group.pdf     | plot\_Fst\_over\_time\_and\_space.R |
| Figure 5               | results/4\_Fst/m3/plot\_individual\_genetic\_distance\_rooted\_tree\_122018ERD.pdf                     | Fst\_analyses\_for\_eelseq.R      |
| Figure 6               | results/6\_PCA/m3/PCA.no.outgroup.122018ERD.pdf                                                  | plot_PCA.R                     |
| Supplementary Figure 1 | results/5\_admixture/m3/all\_individuals/plot\_admixture\_K2.pdf                                    | plot\_admixture\_barchart.R      |
| Supplementary Table 1  | results/7\_ANOVA/m3/table\_Tukey\_ANOA\_p\_values\_11\_japanese\_eel\_pops\_whiteout\_P\_less\_than\_0.01.txt | ANOVA\_script.R                 |
| Supplementary Table 2  | results/4\_Fst/m3/                                                                               | fst\_for\_supplementary\_table\_2  |
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

Make output directories:

```
mkdir -p data/STACKS_processed/4_depth_optimization/m3/
mkdir -p data/STACKS_processed/4_depth_optimization/m6/
mkdir -p data/STACKS_processed/4_depth_optimization/m10/
```

ref_map options:
  
* `-O` = the population map file  
* `--samples` = path to where the .sam files are located (names are read from the population map file)  
* `-b` = the batch ID  
* `-B` = the database name (must end in _radtags)  
* `-D` = description of the data on this run to be displayed (parameters)  
* `--create_db` = create a database if it doesn't exist (--overw_db to overwrite it)  
* `-m` = minimum stack depth  
* `-T` = number of threads to execute  
* `-X` "populations:--fstats" = calculate F statistics for each population  
* `-o` = path to store output files  
* `-e` = path to stacks executables  
 
Run `ref_map` on cbsulm06, varying minimum stack depth (choices: 3, 6, 10) [from `cbsulm06`]:

```
# Read depth of 3 to call stack:
/programs/stacks-1.48/bin/ref_map.pl \
	--samples data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/ \
	-O data/sample_info/population_file_for_stacks.txt \
	-b 1 \
	-B japaneseEel_radtags \
	-D "Population RADseq Samples: m:3" \
	--create_db \
	-m 3 \
	-T 4 \
	-X "populations:--fstats" \
	-o data/STACKS_processed/4_depth_optimization/m3/ \
	-e /programs/stacks-1.48/bin/ 

# Read depth of 6 to call stack:
/programs/stacks-1.48/bin/ref_map.pl \
	--samples data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/ \
	-O data/sample_info/population_file_for_stacks.txt \
	-b 3 \
	-B japaneseEel_radtags \
	-D "Population RADseq Samples: m:6" \
	-m 6 \
	-T 4 \
	-X "populations:--fstats" \
	-o data/STACKS_processed/4_depth_optimization/m6/ \
	-e /programs/stacks-1.48/bin/ 

# Read depth of 10 to call stack:
/programs/stacks-1.48/bin/ref_map.pl \
	--samples data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/ \
	-O data/sample_info/population_file_for_stacks.txt \
	-b 4 \
	-B japaneseEel_radtags \
	-D "Population RADseq Samples: m:10" \
	-m 10 \
	-T 4 \
	-X "populations:--fstats" \
	-o data/STACKS_processed/4_depth_optimization/m10/ \
	-e /programs/stacks-1.48/bin/ 
```

## Step 2: Run corrections module (rxstacks)

Make output folders for rxstacks

```
mkdir -p data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/
mkdir -p data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/
mkdir -p data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/
```

`rxstacks` parameters:

* `-b` = batch id
* `-P` = path to stacks output files to be analyzed
* `-t` = number of threads 
* `--lnl_lim` = minimum log likelihood required to keep a catalog locus
* `--conf_lim` = proportion of loci in population that must be confounded relative to the catalog locus [between 0 and 1]
* `--prune_haplo` = prune out non-biological haplotypes unlikely to occur in the population
* `--model_type` = 'snp' (default), 'bounded', or 'fixed'
* `--alpha` = chi square significance level required to call a het or homozygote [either 0.1 (default), 0.05, 0.01, 0.001]
* `--bound_high` = upper bound for error rate (between 0 an 1 (default))
* `--verbose` = extended logging (forces single-threaded execution)
* `--lnl_dist` = print distribution of mean log likelihoods for catalog loci 

Run `rxstacks` on the three depths:

```
# Read depth of 3 to call stack:
/programs/stacks-1.48/bin/rxstacks -b 1 \
	-P data/STACKS_processed/4_depth_optimization/m3/ \
	-o data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	-t 2 \
	--lnl_lim -10 \
	--conf_lim 0.25 \
	--prune_haplo \
	--model_type bounded \
	--bound_high 0.1 \
   	--lnl_dist \
  	--verbose \
   	--conf_filter
    
# Read depth of 6 to call stack:
/programs/stacks-1.48/bin/rxstacks -b 3 \
	-P data/STACKS_processed/4_depth_optimization/m6/ \
	-o data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	-t 2 \
	--lnl_lim -10 \
	--conf_lim 0.25 \
	--prune_haplo \
	--model_type bounded \
	--bound_high 0.1 \
    --lnl_dist \
    --verbose \
    --conf_filter
    
# Read depth of 10 to call stack:
/programs/stacks-1.48/bin/rxstacks -b 4 \
	-P data/STACKS_processed/4_depth_optimization/m10/ \
	-o data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	-t 2 \
	--lnl_lim -10 \
	--conf_lim 0.25 \
	--prune_haplo \
	--model_type bounded \
	--bound_high 0.1 \
   	--lnl_dist \
   	--verbose \
   	--conf_filter
``` 	

## Step 3: Reassemble the catalog after running corrections module

Create scripts to rerun `cstacks` and `sstacks` after running `rxstacks`:

```
# Read depth of 3 to call stack:
scripts/create_stacks_running_scripts.R \
	--pop_map=data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	--b=1 \
	--p=4 \
	--io_path=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	--outfile=scripts/rerun_cstacks_and_sstacks_optimization_m3.sh 

# Read depth of 6 to call stack:
scripts/create_stacks_running_scripts.R \
	--pop_map=data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	--b=3 \
	--p=4 \
	--io_path=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	--outfile=scripts/rerun_cstacks_and_sstacks_optimization_m6.sh 

# Read depth of 10 to call stack:
scripts/create_stacks_running_scripts.R \
	--pop_map=data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	--b=4 \
	--p=4 \
	--io_path=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	--outfile=scripts/rerun_cstacks_and_sstacks_optimization_m10.sh 

```

Run rerun scripts:

```
chmod +x scripts/rerun*.sh
scripts/rerun_cstacks_and_sstacks_optimization_m3.sh
scripts/rerun_cstacks_and_sstacks_optimization_m6.sh 
scripts/rerun_cstacks_and_sstacks_optimization_m10.sh
```
   
## Step 4: Run populations module on catalog to calculate stats

`populations` parameters:

* `-P` = path to stacks output files to be analyzed
* `-b` = batch ID
* `-M` = path to population map
* `-t` = numer of threads to run in parallele sections of code
* `-s` = output file to SQL database
* `-p` = minimum number of populations a locus must be present in to process a locus.
* `-r` = minimum percentage of individuals in a population required to process a locus for that population
* `-m` = minimum stack depth required for individuals at a locus
* `--write_single_snp` = restrict data analysis to only the first SNP per locus
* `--fstats` = SNP and haplotype-based F-statistics
* `--genomic` = output each nucleotide position (fixed or polymorphic) in all population members to a file
* `--vcf` = output SNPs in VCF
* `--plink` = output genotypes in plink format
* `--verbose` = turn on logging

Run the `populations` module to calculate Fst, etc:

```
# Stack depth = 3:
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	-b 1 \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 3 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--verbose 

# Stack depth = 6:
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	-b 3 \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 6 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--verbose 
  
# Stack depth = 10:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	-b 4 \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 10 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--verbose 
```

## Step 5: Filter out samples with low calls and loci with high depth, rerun populations

First, use vcftools to get a list of SNPs called per individual:

```
# m3
/programs/vcftools-v0.1.14/bin/vcf-stats data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.vcf -p data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.vcf.stats

# m6
/programs/vcftools-v0.1.14/bin/vcf-stats data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3.vcf -p data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3.vcf.stats

# m10
/programs/vcftools-v0.1.14/bin/vcf-stats data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/batch_4.vcf -p data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/batch_4.vcf.stats	
```

Then, check number of SNPs called per sample

Create output directory:

```
mkdir -p results/3_optimizing_depth/
```

```
# m3
scripts/plot_variants_per_sample.R \
	--vcf_stats_counts_file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.vcf.stats.counts \
	--desc=m3 \
	--outpath=results/3_optimizing_depth/

# m6
scripts/plot_variants_per_sample.R \
	--vcf_stats_counts_file=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3.vcf.stats.counts \
	--desc=m6 \
	--outpath=results/3_optimizing_depth/

# m10
scripts/plot_variants_per_sample.R \
	--vcf_stats_counts_file=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/batch_4.vcf.stats.counts \
	--desc=m10 \
	--outpath=results/3_optimizing_depth/
```

For a depth of 3:

![m3](results/3_optimizing_depth/m3_barplot_SNPs_per_sample_in_vcf_before_coverage_filtering_121118ERD.png)

For a depth of 6:

![m6](results/3_optimizing_depth/m6_barplot_SNPs_per_sample_in_vcf_before_coverage_filtering_121118ERD.png)

For a depth of 10:

![m10](results/3_optimizing_depth/m10_barplot_SNPs_per_sample_in_vcf_before_coverage_filtering_121118ERD.png)

Second, examine mean coverage per loci in each depth.
Filter out any locus that is more than 2 x SD higher than the mean coverage:

```
# m3
scripts/filter_SNPs_and_create_whitelist_for_populations.R \
	--vcf.file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.vcf \
	--white.list.out=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	--plot.out=results/3_optimizing_depth/ \
	--desc=m3
	
# m6
scripts/filter_SNPs_and_create_whitelist_for_populations.R \
	--vcf.file=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3.vcf \
	--white.list.out=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	--plot.out=results/3_optimizing_depth/ \
	--desc=m6
	
# m10
scripts/filter_SNPs_and_create_whitelist_for_populations.R \
	--vcf.file=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/batch_4.vcf \
	--white.list.out=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	--plot.out=results/3_optimizing_depth/ \
	--desc=m10 
```

Mean coverage before and after filtering (m3):

![m3](results/3_optimizing_depth/hist_coverage_before_and_after_filtering_m3_121918ERD.png)

m6:

![m6](results/3_optimizing_depth/hist_coverage_before_and_after_filtering_m6_121918ERD.png)

m10:

![m10](results/3_optimizing_depth/hist_coverage_before_and_after_filtering_m10_121918ERD.png)

Third, rerun `populations` using the whitelist to include only loci that are:

* i) biallelic (loci in VCF from previous step are biallelic only) and  
* ii) not of very high coverage (>2SD + mean coverage)  

Make output directories:

```
mkdir data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/
mkdir data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/
mkdir data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/
```

Rerun populations:

```
# Stack depth = 3:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	-b 1 \
	-O data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-W data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 3 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
	
# Stack depth = 6:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	-b 3 \
	-O data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-W data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 6 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
	
# Stack depth = 10:  
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	-b 4 \
	-O data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107.txt \
	-W data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r .6 \
	-m 10 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
```

## Step 6: Run populations module on "meta-populations"

One concern is that 5 individuals per population might be on the low side. 
To address this, samples were divided into four metapopulations, and statistics on those metapopulations calculated using populations. 

The four metapopulations are as follows:

1. “South China Sea” includes Guangdong2009/09
2. “East China Sea” includes Fujian2009/09, Yangtze2005/05, Yangtze2006/06, Yangtze2007/07, Yangtze2008/08, Yangzte2009/09
3. “Yellow Sea” includes Jiangsu2009/09
4. “Pacific Ocean” includes Chibaken2001/01 and Kagawa2001/01)

First, make output directories for results:

```
# Make output directories:
mkdir data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/four_meta_populations/
mkdir data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/four_meta_populations/
mkdir data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/four_meta_populations/
```

Next, run populations for each stack depth:

Stack depth = 3:  

```
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/ \
	-b 1 \
	-O data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 3 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
```

Stack depth = 6:  

```
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ \
	-b 3 \
	-O data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 6 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
```

Stack depth = 10:  

```
/programs/stacks-1.48/bin/populations \
	-P data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/ \
	-b 4 \
	-O data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/four_meta_populations/ \
	-M data/sample_info/population_file_for_stacks_no_JJ-107_four_meta_populations.txt \
	-W data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/whitelist_loci_after_coverage_filtering.txt \
	-t 4 \
	-s \
	-p 1 \
	-r 0.6 \
	-m 10 \
	--write_single_snp \
	--fstats \
	--genomic \
	--vcf \
	--plink \
	--phylip \
	--verbose
```

### Final data sets

The final datasets for each read depth are in the following directories on `cbsulm06`:

* m3: `/workdir/japaneseEel/data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/`  
* m6: `/workdir/japaneseEel/data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/`  
* m10:  `/workdir/japaneseEel/data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/`  

Stacks requires each batch of processing to have a batch number.
m3 was batch\_1, m6 was batch\_3, and m10 was batch\_4. 
In each of those folders, SNP called across all individuals are saved in multiple different formats. 
In particular, the `batch_*.vcf` and `batch_*.plink.raw` files will be useful for further analyses. 
SNPs saved in these files are:

* i) bi-allelic
* ii) in 60% of individuals in at least one population
* iii) filtered for coverage too deep (indicating a duplication/mapping issue)
* iv) filtered so that only one SNP per stack remains (to limit LD, first SNP position retained). 


Additionally, output files from `populations` live in this folder as well, including Fst calculations between all pairs of populations.  

# Run analyses on final data set

Note that all results from analyses are stored in the `results/` directory. 
If analyses were run for each of the three read depths, typically results are stored in sub-directories according to read depth (m3, m6, or m10). 
For population-wide analyses, it seems like having more SNPs called per population gives more reliable results (m3), although results between m3 and m6 are largely consistent. 
Some caution should be noted for examining m10 results.
Several samples have very few SNPs called using this threshold, which affects across population comparisons. 
In particular, pairwise comparisons between some samples could not be made because they simply don't share overlapping SNPs. 
Therefore, I would consider m3 results to be the most robust and would use those in the final analysis. 

## Create basic stats plots

Make maps of where and when samples were taken:

```
scripts/make_map_for_eelseq_project.R 
```

![map](results/1_general_info/map_sample_locations.png)

![timeline](results/1_general_info/timeline_sample_collection.png)

Make plot of the number of loci per population at each stack depth and plot of the nucleotide diversity in each population:

```
scripts/plot_loci_per_population.R \
	--m3=results/3_optimizing_depth/batch_1.sumstats_summary_tail.tsv \
	--m6=results/3_optimizing_depth/batch_3.sumstats_summary_tail.tsv \
	--m10=results/3_optimizing_depth/batch_4.sumstats_summary_tail.tsv \
	--pop_info=data/sample_info/eelseq_sample_info_degrees_removed_corrected_with_names_coordinates_fixed.txt \
	--outpath=results/3_optimizing_depth/
```

![loci per pop](results/3_optimizing_depth/barplot_variant_sites_per_pop_across_min_stack_depths.png)

![pi per pop](results/3_optimizing_depth/barplot_diversity_by_population_m3.png)

Look at genetic distance versus geographic distance and temporal distance.
This script uses a Mantel test to assess significance:

```
scripts/plot_Fst_over_time_and_space.R \
	--sumstats_summary=results/3_optimizing_depth/batch_1.sumstats_summary_tail.tsv \
	--fst_summary=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.fst_summary.tsv \
	--desc=m3 \
	--outpath=results/3_optimizing_depth/ 
```

There is a trend between genetic distance (Fst/(1 - Fst)) and geographic distance (lm p = 0.05907), but it is opposite of what we would expect to see: 

![geographic vs genetic distance](results/3_optimizing_depth/m3_plot_geographic_distance_by_genetic_distance_no_out_group.png)

There is not a relationship between genetic distance and temporal distance (lm p = 0.2369):

![temporal vs. genetic distance](results/3_optimizing_depth/m3_plot_time_distance_by_genetic_distance_yangzte.png)

## Pull data for Table 2

The first row of Table 2 comes from the log files from `process_radtags`.
Here's information from each of the log files:

From `/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out_32_36_manli2/process_radtags.fastq_files_updated.log`:

```
File    Retained Reads  Low Quality     Ambiguous Barcodes      Ambiguous RAD-Tag       Total
32_36_manli2_R1_merged.fastq.gz 55963714        253415  7841343 7186317 71244789

Total Sequences 71244789
Ambiguous Barcodes      7841343
Low Quality     253415
Ambiguous RAD-Tag       7186317
Retained Reads  55963714
```

From `/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/0_process_radtags_out/process_radtags.fastq_files_updated.log`:

```
File    Retained Reads  Low Quality     Ambiguous Barcodes      Ambiguous RAD-Tag       Total
32manli1_R1.fastq.gz    68699159        423170  9247278 10371910        88741517

Total Sequences 88741517
Ambiguous Barcodes      9247278
Low Quality     423170
Ambiguous RAD-Tag       10371910
Retained Reads  68699159
```

Sum those up:
Total Sequences = 159,986,306
Ambiguous Barcodes = 17,088,621 
Low Quality = 676,585 
Ambiguous RAD-Tag = 17,558,227 
Retained Reads = 124,662,873 

The second row of Table 2 comes from the Bowtie .log files. 
Run this script to summarize the .log files:

```
scripts/analyze_bowtie_logs.R \
	--inpath=/fs/cbsufsrv5/data2/japaneseEel/data/STACKS_processed/1_aligned_to_nuclear_genome_bowtie/ \
	--outpath=results/1_general_info/
```

The final row of Table 2 comes from counting the lines in the catalog files generated by Stacks and the lines of the final VCF file (-11):

```
wc -l data/STACKS_processed/4_depth_optimization/m3/batch_1.catalog.tags.tsv
# 214,210

wc -l data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.catalog.tags.tsv
# 213,610

wc -l data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.vcf
# 106663 - 11 = 106,652
 
```

## Create phylogenetic trees from population level Fst

First, recode plink files to be in dosage format: 

```
# m3:
/programs/plink-1.07-x86_64/plink --file data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink --recodeA --noweb --out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink

# m6:
/programs/plink-1.07-x86_64/plink --file data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink --recodeA --noweb --out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink

# m10:
/programs/plink-1.07-x86_64/plink --file data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink --recodeA --noweb --out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink
```

Make directories to store analysis output:

```
mkdir -p results/4_Fst/m3/
mkdir -p results/4_Fst/m6/
mkdir -p results/4_Fst/m10/
```

Next, calculate the genetic distance between all pairs of individuals: 

```
# m3
scripts/create_individual_distance_matrix.R \
	--raw_file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.raw \
	--cores=12 \
	--outpath=results/4_Fst/m3/

# m6	
scripts/create_individual_distance_matrix.R \
	--raw_file=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.raw \
	--cores=12 \
	--outpath=results/4_Fst/m6/

# m10	
scripts/create_individual_distance_matrix.R \
	--raw_file=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.raw \
	--cores=12 \
	--outpath=results/4_Fst/m10/
```

Create trees based on pairwise Fst and individual genetic distances:

```
# m3
scripts/Fst_analyses_for_eelseq.R \
	--fst_summary=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/batch_1.fst_summary.tsv \
	--desc=m3 \
	--ind_dist=results/4_Fst/m3/table_pairwise_individual_genetic_distances.txt \
	--outpath=results/4_Fst/m3/

# m6
scripts/Fst_analyses_for_eelseq.R \
	--fst_summary=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3.fst_summary.tsv \
	--desc=m6 \
	--ind_dist=results/4_Fst/m6/table_pairwise_individual_genetic_distances.txt \
	--outpath=results/4_Fst/m6/
	
# m10	
scripts/Fst_analyses_for_eelseq.R \
	--fst_summary=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/batch_4.fst_summary.tsv \
	--desc=m10 \
	--ind_dist=results/4_Fst/m10/table_pairwise_individual_genetic_distances.txt \
	--outpath=results/4_Fst/m10/
```
	
Using the stack depth=3 data (m3), we can see that the outgroup clusters separately and that there is very little divergence between any of the eel species.
In addition, samples don't necessarily cluster by population/sampling date. 

![rooted tree, created using UPGMA (m3 depth)](results/4_Fst/m3/plot_Fst_rooted_tree_all_pops_122018ERD.png) 
"rooted tree, created using UPGMA (m3 depth)"

Above, rooted tree created using UPGMA (m3 depth).

![unrooted tree, created using neighbor-joining (m3 depth)](results/4_Fst/m3/plot_Fst_unrooted_tree_all_pops_122018ERD.png) 
"unrooted tree, created using neighbor-joining (m3 depth)"

Above, unrooted tree created using neighbor-joining (m3 depth).

Additionally, examined pairwise genetic distances between all individuals.
To do this, for all SNPs that a pair shared in common, the number of differences were summed and divided by 2/number of shared SNPs. 
This created a distance matrix that was then used to build the trees below. 
At the individual level, no apparent clustering by either population or year of sampling.  

![unrooted tree, individual genetic distances created using neighbor-joining (m3 depth)](results/4_Fst/m3/plot_individual_genetic_distance_rooted_tree_122018ERD.png) "unrooted tree, created using neighbor-joining (m3 depth")
Above, rooted tree of individual genetic distances (m3 depth)

![rooted tree, individual genetic distances created using neighbor-joining (m3 depth)](results/4_Fst/m3/plot_individual_genetic_distance_unrooted_tree_122018ERD.png) "unrooted tree, created using neighbor-joining (m3 depth")

Above, unrooted tree of individual genetic distances (m3 depth)

## Run admixture

Run admixture to discover population substructure/admixing. 


### Admixture - all individuals 


First, need to remove crazy eel contig names as chromosome or admixture won't like the file:

```
scripts/convert_plink_chr_names_for_admixture.R --file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.map
scripts/convert_plink_chr_names_for_admixture.R --file=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.map
scripts/convert_plink_chr_names_for_admixture.R --file=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.map
```

Second, copy ped file with new name:

```
cp data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.ped data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.ped
cp data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.ped data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.ped
cp data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.ped data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.ped
```

Third, use plink to convert to bed/bim/fam:

```	
/programs/plink-1.9-x86_64-beta3.30/plink \
	--file data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture \
	--make-bed \
	--allow-extra-chr \
	--out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture
	
/programs/plink-1.9-x86_64-beta3.30/plink \
	--file data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture \
	--make-bed \
	--allow-extra-chr \
	--out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture
	
/programs/plink-1.9-x86_64-beta3.30/plink \
	--file data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture \
	--make-bed \
	--allow-extra-chr \
	--out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture
```

Switch to the correct working directory and run admixture with the following parameters:

* `-j4` = use four processors 
* `-s` = set seed (1)
* `--cv` = do cross validation
* `-C` = set min delta to hit before declaring convergence (if float) or max iterations (int)
* `$K` = number of populations (K)

For m3:

```
# switch to running directory:
cd results/5_admixture/m3/all_individuals/

# run admixture using different values of K:
for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

For m6:

```
cd results/5_admixture/m6/all_individuals/

for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

m10: 

```
cd results/5_admixture/m10/all_individuals/

for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

Plot optimal number of K for each depth:

```
for i in m3 m6 m10
do
scripts/plot_admixture_CVs.R \
	--cv_file=results/5_admixture/$i/all_individuals/admixture_cross_validations.txt  \
	--outfile=results/5_admixture/$i/all_individuals/plot_admixture_cross_validations_${i}.pdf
done
```

Optimal number of K for m3:
![m3](results/5_admixture/m3/all_individuals/plot_admixture_cross_validations_m3.png)

Optimal number of K for m6:
![m6](results/5_admixture/m6/all_individuals/plot_admixture_cross_validations_m6.png)

Optimal number of K for m10:
![m10](results/5_admixture/m10/all_individuals/plot_admixture_cross_validations_m10.png)

The optimal number of populations for m3 and m6 is 2, which makes sense given the two species. 
The optimal number of populations for m10 is 4 (although 5 was not tested). 

Create admixture barplots for optimal K for each depth. 

m3:

```
scripts/plot_admixture_barchart.R \
	--Q_file=results/5_admixture/m3/all_individuals/batch_1.plink.for.admixture.2.Q \
	--fam_file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.fam \
	--K=2 \
	--outfile=results/5_admixture/m3/all_individuals/plot_admixture_K2.pdf
```

![m3 K2 admixture plot](results/5_admixture/m3/all_individuals/plot_admixture_K2.png)

m6

```
scripts/plot_admixture_barchart.R \
	--Q_file=results/5_admixture/m6/all_individuals/batch_3.plink.for.admixture.2.Q \
	--fam_file=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.fam \
	--K=2 \
	--outfile=results/5_admixture/m6/all_individuals/plot_admixture_K2.pdf
```

![m6 K2 admixture plot](results/5_admixture/m6/all_individuals/plot_admixture_K2.png)

m10

```
scripts/plot_admixture_barchart.R \
	--Q_file=results/5_admixture/m10/all_individuals/batch_4.plink.for.admixture.4.Q \
	--fam_file=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.fam \
	--K=4 \
	--outfile=results/5_admixture/m10/all_individuals/plot_admixture_K4.pdf
```

![m10 K4 admixture plot](results/5_admixture/m10/all_individuals/plot_admixture_K4.png)

### Admixture - no outgroup


First, use plink to exclude the Hainan individuals

```
/programs/plink-1.9-x86_64-beta3.30/plink --file data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt --out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture_no_outgroup --make-bed
/programs/plink-1.9-x86_64-beta3.30/plink --file data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt --out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture_no_outgroup --make-bed
/programs/plink-1.9-x86_64-beta3.30/plink --file data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt --out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture_no_outgroup --make-bed
```

Second, eliminate SNPs with all missing genotypes from analysis

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture_no_outgroup --geno 0.95 --make-bed --out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture_no_outgroup_genotype_filtered
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture_no_outgroup --geno 0.95 --make-bed --out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture_no_outgroup_genotype_filtered
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture_no_outgroup --geno 0.95 --make-bed --out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture_no_outgroup_genotype_filtered
```

Third, make output directories:

```
mkdir -p results/5_admixture/m3/no_outgroup/
mkdir -p results/5_admixture/m6/no_outgroup/
mkdir -p results/5_admixture/m10/no_outgroup/
```

Fourth, run admixture. Admixture program will only save output to working directory, so create working directories for each iteration:

* `-j4` = use four processors 
* `-s` = set seed (1)
* `--cv` = do cross validation
* `-C` = set min delta to hit before declaring convergence (if float) or max iterations (int)
* `$K` = number of populations (K)

Switch to the correct working directory and run admixture with the following parameters:

For m3:

```
cd results/5_admixture/m3/no_outgroup/

for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture_no_outgroup_genotype_filtered.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

For m6:

```
cd results/5_admixture/m6/no_outgroup/

for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture_no_outgroup_genotype_filtered.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

m10 

```
cd results/5_admixture/m10/no_outgroup/

for K in 1 2 3 4 6 8 10 
do
/programs/admixture_linux-1.23/admixture -j4 -s 1 -C 0.01 --cv ../../../../data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture_no_outgroup_genotype_filtered.bed $K | tee log${K}.out
done

#...Pull together a file of cross-validation errors
grep -h CV log*.out > admixture_cross_validations.txt

#...switch back to base directory
cd ../../../..
```

Plot the cross validation rate for the different values of K:

```
for i in m3 m6 m10
do
scripts/plot_admixture_CVs.R \
	--cv_file=results/5_admixture/$i/no_outgroup/admixture_cross_validations.txt  \
	--outfile=results/5_admixture/$i/no_outgroup/plot_admixture_cross_validations_${i}.pdf
done
```

Cross validation error rate for m3:

![CV m3 no outgroup](results/5_admixture/m3/no_outgroup/plot_admixture_cross_validations_m3.png)

Cross validation error rate for m6:

![CV m6 no outgroup](results/5_admixture/m6/no_outgroup/plot_admixture_cross_validations_m6.png)

Cross validation error rate for m10:
![CV m10 no outgroup](results/5_admixture/m10/no_outgroup/plot_admixture_cross_validations_m10.png)

For all stack depths (3, 6, and 10), the lowest cross-validation error rate is seen at K = 1, lending support for panmixia. 

## Run PCA

Use Plink to run PCA for all individuals: 

m3:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture --pca 59 --out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.pca 
```

m6:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture --pca 59 --out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.pca 
```
	
m10:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture --pca 59 --out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.pca 
```

Perform PCA, excluding Hainan population (the outgroup)

m3:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture --pca 54 --out data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.pca.no.hainan --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt
```

m6:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture --pca 54 --out data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.pca.no.hainan --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt
```

m10:

```
/programs/plink-1.9-x86_64-beta3.30/plink --bfile data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture --pca 54 --out data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.pca.no.hainan --remove results/1_general_info/hainan_individuals_to_remove_in_plink.txt
```

Plot first two principal components:

```
mkdir -p results/6_PCA/m3/
mkdir -p results/6_PCA/m6/
mkdir -p results/6_PCA/m10/
```

m3:

```
scripts/plot_PCA.R \
	--full_plink_in=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.pca \
	--part_plink_in=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.plink.for.admixture.pca.no.hainan \
	--outpath=results/6_PCA/m3/	

```

![m3 PCA all samples](results/6_PCA/m3/PCA.all.samples.122018ERD.png)

![m3 PCA no outgroup](results/6_PCA/m3/PCA.no.outgroup.122018ERD.png)

m6:

```
scripts/plot_PCA.R \
	--full_plink_in=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.pca \
	--part_plink_in=data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/coverage_filtered/batch_3.plink.for.admixture.pca.no.hainan \
	--outpath=results/6_PCA/m6/
```

![m6 PCA all samples](results/6_PCA/m6/PCA.all.samples.122018ERD.png)

![m6 PCA no outgroup](results/6_PCA/m6/PCA.no.outgroup.122018ERD.png)

m10:

```
scripts/plot_PCA.R \
	--full_plink_in=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.pca \
	--part_plink_in=data/STACKS_processed/4_depth_optimization/m10/rxstacks_corrected/coverage_filtered/batch_4.plink.for.admixture.pca.no.hainan \
	--outpath=results/6_PCA/m10/
```

![m10 PCA all samples](results/6_PCA/m10/PCA.all.samples.122018ERD.png)

![m10 PCA no outgroup](results/6_PCA/m10/PCA.no.outgroup.122018ERD.png)

## Calculate significant differences in diversity statistics using ANOVA

The following script uses an ANOVA followed by TukeyHSD to identify significant differences between each pair of populations in the study: 

```
scripts/ANOVA_script.R \
	--sumstats_file=data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.sumstats.tsv \
	--outpath=results/7_ANOVA/m3/
```

## Marker vs. Sampling size simulations

The following script (from Andy Clark) demonstrates that Fst is accurately calculated using large numbers of variants, even with small population sizes (n = 5, as in this manuscript):

```
scripts/eelsim.txt
```

## Calculate stats for manuscript

```
scripts/stats_in_paper.R
```

This produces a text table at `results/1_general_info/stats_for_paper.txt` that lists all of the statistics reported in the text of the manuscript. 
 


