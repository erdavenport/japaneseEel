#!/usr/bin/env Rscript

##### 
# This script contains functions that are used across scripts in the eelseq study. 
# It can be loaded within a script that calls these functions by running: source("common_functions_for_eelseq_analyses.R")
#####



##### Sample information table
samples <- read.table("data/sample_info/eelseq_sample_info_degrees_removed_corrected_with_names_coordinates_fixed.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

samples$latitude <- as.numeric(gsub("N", "", samples$latitude))

pops <- unique(samples[,-which(colnames(samples) %in% c("sample_ID"))])
pops$lat_help <- c(4, 4, 4, 4, 4, 7, 6, 5, 3, 1, 2, 8)
pops$popName <- paste0(pops$pop_name, " - ", pops$lat_help, " - ", pops$year)
pops$popNameY <- paste0(pops$year, " - ", pops$lat_help, " - ", pops$pop_name)

samples$popName <- pops$popName[match(samples$population, pops$population)]
samples$popNameY <- pops$popNameY[match(samples$population, pops$population)]



##### Function for population colors/shapes:
whoIsMyPop <- function(popName) {
	p <- samples[which(samples$pop_name == popName),]
	mypop <- c(popName, unique(p$population), p$sample_ID)
	return(mypop)
}

pop.cols <- function(x) {
	if (x %in% whoIsMyPop("Chiba-ken")) {
		return("firebrick")
	} else if (x %in% whoIsMyPop("Kagawa")) {
		return("firebrick1")
	} else if (x %in% whoIsMyPop("Jiangsu province")) {
		return("blue4")
	} else if (x %in% whoIsMyPop("Yangtze River estuary")) {
		return("blue") 
	} else if (x %in% whoIsMyPop("Zhejiang province")) {
		return("royalblue3")
	} else if (x %in% whoIsMyPop("Fujian province")) {
		return("deepskyblue")
	} else if (x %in% whoIsMyPop("Guangdong province")) {
		return("cyan2")
	} else if (x %in% whoIsMyPop("Hainan province")) {
		return("lightblue1")
	}
}

pop.pch <- function(x) {
	if (x %in% whoIsMyPop("Chiba-ken")) {
		return(0)
	} else if (x %in% whoIsMyPop("Kagawa")) {
		return(1)
	} else if (x %in% whoIsMyPop("Jiangsu province")) {
		return(19)
	} else if (x %in% whoIsMyPop("Yangtze River estuary")) {
		return(18) 
	} else if (x %in% whoIsMyPop("Zhejiang province")) {
		return(17)
	} else if (x %in% whoIsMyPop("Fujian province")) {
		return(20)
	} else if (x %in% whoIsMyPop("Guangdong province")) {
		return(15)
	} else if (x %in% whoIsMyPop("Hainan province")) {
		return(13)
	}
}

# Original pop cols function (from before 11/4/16)
pop.cols.ORIGINAL <- function(x) {
	if (x %in% c("Chiba-ken", "Q")) {
		return("firebrick")
	} else if (x %in% c("Kagawa", "H")) {
		return("firebrick1")
	} else if (x %in% c("Jiangsu province", "DF")) {
		return("blue4")
	} else if (x %in% c("Yangtze River estuary", "CXD", "DH", "LCG07", "JDS", "LCG09")) {
		return("blue") 
	} else if (x %in% c("Zhejiang province", "CX")) {
		return("royalblue3")
	} else if (x %in% c("Fujian province", "FQ")) {
		return("deepskyblue")
	} else if (x %in% c("Guangdong province", "XH")) {
		return("cyan2")
	} else if (x %in% c("Hainan province", "HM")) {
		return("lightblue1")
	}
}