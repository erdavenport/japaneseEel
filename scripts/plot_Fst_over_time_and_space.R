#!/usr/bin/env Rscript

######  
suppressMessages(library("docopt"))
"
Usage:
	plot_Fst_over_time_and_space.R --sumstats_summary=<sumstats_summary> --fst_summary=<fst_summary> --desc=<desc> --outpath=<outpath> 

Description: This script will generate a shell script that will run stacks using the parameters specified on command line. 

Options:
	--sumstats_summary=<sumstats_summary> 	from populations module
	--fst_summary=<fst_summary>				from populations module
	--desc=<desc>							decription of analysis (m3)
	--outpath=<outpath> 					path to save output files
" -> doc
######

###### PARAMETERS ##########
# Set the parameters:
today <- Sys.Date()											# Set the date that will go on the end of the files generated by this script
today <- format(today, format="%m%d%y")
#############################



##### Load libraries:
suppressMessages(library("dplyr"))
library("reshape2")
library("ggplot2")
suppressMessages(library("geosphere"))
library("testit")
library("ade4")
library("ape")



##### Load arguments:
opts <- docopt(doc)
#print(opts)
sumstats_file <- opts$sumstats_summary
fst_file <- opts$fst_summary
mydesc <- opts$desc
outpath <- opts$outpath



#sumstats_file <- "results/3_optimizing_depth/batch_1.sumstats_summary_tail.tsv"
#fst_file <- "data/STACKS_processed/4_depth_optimization/m3/rxstacks_corrected/coverage_filtered/batch_1.fst_summary.tsv"
#mydesc <- "m3"
#outpath <- "results/3_optimizing_depth/"



##### Load data:
samples <- read.table("data/sample_info/eelseq_sample_info_degrees_removed_corrected_with_names_coordinates_fixed.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
samples$latitude <- as.numeric(gsub("N", "", samples$latitude))
sumstats <- read.table(file=sumstats_file, sep="\t", header=TRUE)

fst <- read.table(file=fst_file, sep="\t", header=TRUE)
row.names(fst) <- fst$X
fst <- fst[,-1]



##### Pull out the location data:
reps <- samples[seq(1,60,5),]
reps$source <- paste(reps$pop_name, reps$year)
pops <- reps %>% 
	filter(population %in% c("LCG09", "XH", "FQ", "CX", "DF", "Q", "H", "HM"))

# Reorder levels for plot
pops$pop_name <- factor(pops$pop_name, levels = c("Chiba-ken", "Kagawa", "Jiangsu province", "Yangtze River estuary", "Zhejiang province", "Fujian province", "Guangdong province", "Hainan province"))



##### Create matrix of distances:
dists <- distm(reps[,c("longitude", "latitude")])
colnames(dists) <- reps$population
rownames(dists) <- reps$population
dists[lower.tri(dists)] <- NA
diag(dists) <- NA



##### Plot genetic distance (Fst/(1 - Fst)) x distance (exclude outgroup)
no_out_fst <- fst[,-which(colnames(fst) %in% "HM")]
no_out_dists <- dists[-which(rownames(dists) %in% "HM"), -which(colnames(dists) %in% "HM")]
no_out_fst_vec <- as.vector(no_out_fst)
no_out_dists_vec <- as.vector(no_out_dists)
myyprime <- no_out_fst_vec[!is.na(no_out_fst_vec)]
myx <- no_out_dists_vec[!is.na(no_out_dists_vec)]/1000

# Calculate genetic distance
myy <- myyprime/(1-myyprime)

# Find linear relationship between geographic and genetic distance
mylm <- lm(myy ~ myx)

# Write summary of lm to file:
sink(paste0(outpath, mydesc, "_genetic_vs_geographic_distance_lm_summary.txt"))
summary(mylm)
sink()

pdf(paste0(outpath, mydesc, "_plot_geographic_distance_by_genetic_distance_no_out_group.pdf"), width=6, height=6)
plot(myx, myy, ylab="Fst/(1 - Fst)", xlab="geographic distance (km)", pch=16, col="gray40")
abline(mylm, col="red")
hi <- dev.off()



##### Plot genetic distance versus time for the Yangzee samples:
distsp <- dists[-12,]
assert(colnames(distsp) == colnames(fst))
assert(rownames(distsp) == rownames(fst))

# melt:
md <- melt(distsp)
mf <- melt(as.matrix(fst))

# combine:
assert(md$Var1 == mf$Var1)
assert(md$Var2 == mf$Var2)
mf$dist <- md$value

# remove NAs
mfnona <- mf[-which(is.na(mf$value)), ]

# only keep Yangtze river data points:
yangtze_info <- reps[which(reps$pop_name == "Yangtze River estuary"), ]

yangtzemf <- mfnona %>%
	filter(Var1 %in% yangtze_info$population & Var2 %in% yangtze_info$population)
	
# Add year info to table
yangtzemf$Var1Year <- yangtze_info$year[match(yangtzemf$Var1, yangtze_info$population)]
yangtzemf$Var2Year <- yangtze_info$year[match(yangtzemf$Var2, yangtze_info$population)]

# find time difference:
yangtzemf$time <- abs(yangtzemf$Var1Year - yangtzemf$Var2Year)

# make year label:
yangtzemf$label <- paste0(yangtzemf$Var1Year, "-", yangtzemf$Var2Year)

# calculate genetic distance:
yangtzemf$gendist <- yangtzemf$value/(1 - yangtzemf$value)

# Find relationship
mlm2 <- lm(yangtzemf$gendist ~ yangtzemf$time)

# Write summary of lm to file:
sink(paste0(outpath, mydesc, "_genetic_vs_temporal_distance_lm_summary.txt"))
summary(mlm2)
sink()

# plot:
pdf(paste0(outpath, mydesc, "_plot_time_distance_by_genetic_distance_yangzte.pdf"), width=6, height=6)
plot(yangtzemf$time, yangtzemf$gendist, ylab="Fst/(1 - Fst)", xlab="temporal distance (years)", pch=16, col="gray40")
abline(mlm2, col="red")
hi <- dev.off()



##### Mantel tests

### For just the Yangtze temporal samples, see if there's a relationship:

# Make wide matricies for geographic and time distance:
gendist <- acast(yangtzemf, Var1 ~ Var2, value.var = "gendist")
tiiiime <- acast(yangtzemf, Var1 ~ Var2, value.var = "time")

# Make square with all values:
gendist2 <- cbind(NA, gendist)
colnames(gendist2)[1] <- "CXD"
tiiiime2 <- cbind(NA, tiiiime)
colnames(tiiiime2)[1] <- "CXD"

gendist3 <- rbind(gendist2, NA)
rownames(gendist3)[5] <- "LCG09"
tiiiime3 <- rbind(tiiiime2, NA)
rownames(tiiiime3)[5] <- "LCG09"

# Fill in diagonals:
diag(gendist3) <- 0
diag(tiiiime3) <- 0

# Fill in bottom tri:
gendist3[lower.tri(gendist3)] <- t(gendist3)[lower.tri(gendist3)]
tiiiime3[lower.tri(tiiiime3)] <- t(tiiiime3)[lower.tri(tiiiime3)]

y_mantel <- mantel.test(gendist3, tiiiime3, graph = FALSE)

yangze_temporal_samples_p <- y_mantel$p

### Genetic vs. Geographic distances
# Keep only the Yangze River 2009 and remove the outgroup:
rmme <- c("CXD", "DH", "LCG07", "JDS", "HM")

fstp <- fst[-which(colnames(fst) %in% rmme), -which(rownames(fst) %in% rmme)]
distsp <- dists[-which(colnames(dists) %in% rmme), -which(rownames(fst) %in% rmme)]

# Make fstp genetic distance:
genp <- fstp/(1 - fstp)

# Add bottom row:
genp <- rbind(genp, NA)
distsp <- rbind(distsp, NA)

# Fill in diagonals:
diag(genp) <- 0
diag(distsp) <- 0

# Fill in bottom of triangle:
genp[lower.tri(genp)] <- t(genp)[lower.tri(genp)]
distsp[lower.tri(distsp)] <- t(distsp)[lower.tri(distsp)]

g_mantel <- mantel.test(genp, distsp, graph = FALSE)

geographic_samples_p <- g_mantel$p

mantelps <- cbind(c("temporal_p", "geographic_p"), c(yangze_temporal_samples_p, geographic_samples_p))

write.table(mantelps, paste0(outpath, "table_mantel_test_p_values.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



print("DONE!")