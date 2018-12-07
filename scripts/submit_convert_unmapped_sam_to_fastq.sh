#!/bin/bash

#$ -S /bin/bash
#$ -q regular.q 
#$ -j n
#$ -N get_unmapped_reads
#$ -p -200

# script inputs:
# $1 path to sam file to be parsed (assuming mount!)
# $2 sam file name
# $3 path to write out fastq file of unmapped reads

# date
d1=$(date +%s)
echo $HOSTNAME
echo $2

# Make a directory on the local drive to run computations:
mkdir -p /SSD/$USER/$JOB_ID/

# Mount fsvr5 to get access to EELseq files:
/programs/bin/labutils/mount_server cbsufsrv5 /data2

# Copy input files to local folder:
cp $1/$2 /SSD/$USER/$JOB_ID/myfile.sam		# sam file

# Save prefix of file to variable to write out:
prefix=$(echo $2 | sed 's/.sam//')
echo $prefix

# Switch to that directory to run fastx:
cd /SSD/$USER/$JOB_ID/

# Convert unmapped reads in sam file back to fastq:
/programs/samtools-1.3/bin/samtools view -b -f 4 myfile.sam > myfile.bam
/programs/samtools-1.3/bin/samtools bam2fq myfile.bam > ${prefix}_unmapped.fastq

# Move out written files to specified location:
cp *.fastq $3

# Remove the directory on the local machine:
rm -r /SSD/$USER/$JOB_ID

# Print out total run time to out file:
d2=$(date +%s)
sec=$(( ( $d2 - $d1 ) ))
hour=$(echo - | awk '{ print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)