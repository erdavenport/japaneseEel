#!/usr/bin/env Rscript

######  
suppressMessages(library("docopt"))
"
Usage:
	convert_plink_chr_names_for_admixture.R --file=<file> 

Description: This script will change the eel chromosome names to numbers so that admixture can run. Note: the new file will be save in the same directory with plink.map replaced with plink.for.admixture.map

Options:
	--file=<file>				plink map file
" -> doc
######

###### PARAMETERS ##########
# Set the parameters:
today <- Sys.Date()											# Set the date that will go on the end of the files generated by this script
today <- format(today, format="%m%d%y")
#############################


##### Load arguments:
opts <- docopt(doc)
file <- opts$file


##### Read file, change chromosome names all to "9", and save new file
# Read in plink map file
plink <- read.table(file, sep="\t", header=FALSE)

# Change chromosome names to 9:
plink$V1 <- "9"

# create new output name
out <- gsub(".plink.map", ".plink.for.admixture.map", file) 

write.table(plink, out, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

print("plink file for admixture created")