#!/bin/bash

#$ -S /bin/bash
#$ -q regular.q 
#$ -j n
#$ -N trim_pe_me
#$ -p -200

# script inputs:
# $1 path to fastq file assoc files (assuming mount!)
# $2 fastq file name
# $3 path to save trimmed fastqs to

# date
d1=$(date +%s)
echo $2

# Make a directory on the local drive to run computations:
mkdir -p /SSD/$USER/$JOB_ID/

# Mount fsvr5 to get access to EELseq files:
/programs/bin/labutils/mount_server cbsufsrv5 /data2

# Copy input files to local folder:
cp $1/$2 /SSD/$USER/$JOB_ID/ 		# fastq file

# Switch to that directory to run fastx:
cd /SSD/$USER/$JOB_ID/
mkdir out

# Run barcode splitter:
/programs/bin/fastx/fastx_trimmer \
-f 1 \
-l 100 \
-Q33 \
-i $2 \
-o out/$2

# Move out written files to specified location:
cp out/* $3

# Remove the directory on the local machine:
rm -r /SSD/$USER/$JOB_ID

# Print out total run time to out file:
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)