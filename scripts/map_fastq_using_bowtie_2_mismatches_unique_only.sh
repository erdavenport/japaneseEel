#!/bin/bash

#$ -S /bin/bash
#$ -q regular.q 
#$ -j n
#$ -pe bscb 2
#$ -N map_me
#$ -p -200

# script inputs:
# $1 path to fastq file to be mapped (assuming mount!)
# $2 fastq file name
# $3 path to bowtie index files
# $4 prefix of bowtie index files  (a_japonica)
# $5 where to write out sam and log files

# date
d1=$(date +%s)
echo $HOSTNAME
echo $2

# Make a directory on the local drive to run computations:
mkdir -p /SSD/$USER/$JOB_ID/

# Mount fsvr5 to get access to EELseq files:
/programs/bin/labutils/mount_server cbsufsrv5 /data2

# Copy input files to local folder:
cp $1/$2 /SSD/$USER/$JOB_ID/ 		# fastq file
cp -r $3 /SSD/$USER/$JOB_ID/bowtie_index/ # barcode file

# Save prefix of file to variable to write out:
prefix=$(echo $2 | sed 's/.fq//')
echo $prefix

# Switch to that directory to run fastx:
cd /SSD/$USER/$JOB_ID/

# Map reads using bowtie
/programs/bowtie-1.1.2/bowtie -S -p 4 -v 2 -m 1 --seed 1 \
	bowtie_index/$4 \
	$2 > \
	${prefix}.sam 2> \
	${prefix}.log


# Move out written files to specified location:
cp *.sam $5
cp *.log $5

# Remove the directory on the local machine:
rm -r /SSD/$USER/$JOB_ID

# Print out total run time to out file:
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)