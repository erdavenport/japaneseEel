#!/usr/bin/env Rscript

######  
suppressMessages(library("docopt"))
"
Usage:
	Fst_analyses_for_eelseq.R --fst_summary=<fst_summary> --ind_dist=<ind_dist> --desc=<desc> --outpath=<outpath> 

Description: This script will create UPGMA and NJ trees from pairwise 1) Fst and 2) individual genetic differences. 

Options:
	--fst_summary=<fst_summary>				from populations module
	--ind_dist=<ind_dist>					table of pairwise genetic distances by individual
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
suppressMessages(library("ape"))
library(phangorn)
library(ggplot2)
suppressMessages(library(ggtree))
suppressMessages(library(phylobase))
library(ggrepel)
library(gridExtra)
library(testit)

# Load common functions:
source("scripts/common_functions_for_eelseq_analyses.R")



##### Load arguments:
opts <- docopt(doc)
#print(opts)
fst_file <- opts$fst_summary
ind_dist <- opts$ind_dist
mydesc <- opts$desc
outpath <- opts$outpath


#fst_file <- "../data/STACKS_processed/7_depth_optimization/m3/rxstacks_corrected/batch_2.fst_summary.tsv"
#ind_dist <- "../results/6_Fst/m3/table_pairwise_individual_genetic_distances.txt"
#mydesc <- "m3"
#outpath <- "../results/6_Fst/m3/"



##### Generate output folder if it isn't there already:
if (!file.exists(outpath)) {
	print(paste0("creating ",outpath," in filesystem"))
	dir.create(file.path(outpath))
}



##### Load data:
samples <- read.table("data/sample_info/eelseq_sample_info_degrees_removed_corrected_with_names_coordinates_fixed.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
fst <- read.table(file=fst_file, sep="\t", header=TRUE)
dists <- read.table(file=ind_dist, sep="\t", header=TRUE)



##### Reformat sample data:
samples$latitude <- as.numeric(gsub("N", "", samples$latitude))
pops <- unique(samples[,-which(colnames(samples) %in% c("sample_ID"))])
pops$lat_help <- c(4, 4, 4, 4, 4, 7, 6, 5, 3, 1, 2, 8)
pops$popName <- paste0(pops$pop_name, " - ", pops$lat_help, " - ", pops$year)
pops$popName2 <- paste0(pops$year, " - ", pops$pop_name)



##### Reformat Fst data:
row.names(fst) <- fst$X
fst <- fst[,-1]

fst2 <- rbind(fst, rep(NA, dim(fst)[2]))
rownames(fst2) <- colnames(fst2)
fst2[lower.tri(fst2)] <- t(fst2)[lower.tri(fst2)]
diag(fst2) <- 0
ind <- match(rownames(fst2), pops$population)
#rownames(fst2) <- pops$popName[ind]




##### All samples
# UPGMA and NJ tree on everything:
treeUPGMA <- upgma(as.dist(fst2))
treeNJ <- NJ(as.dist(fst2))
	
##reformat to phylo4 for plotting
# UPGMA
g1U = as(treeUPGMA, 'phylo4')
dU = data.frame(color=sapply(rownames(fst2), pop.cols))
rownames(dU) = treeUPGMA$tip.label
g2U = phylo4d(g1U, dU)

# NJ
g1NJ <- as(treeNJ, 'phylo4')
dNJ <- data.frame(color=sapply(rownames(fst2), pop.cols))
rownames(dNJ) <- treeNJ$tip.label
g2NJ <- phylo4d(g1NJ, dNJ)

## Plot using ggtree
# Plot rooted UPGMA:
labs <- pops$popName
treelabs <- paste0("   ", c(labs, rep(NA, length(labels(g2U, "internal")))))

plot1 <- ggtree(g2U, ladderize = TRUE) +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_tiplab(aes(label=treelabs), size=3.5, color="gray20") +
	ggtitle("Rooted tree - Fst - (UPGMA)") +
	xlim(0, 0.40)
	
# Plot unrooted neighbor-joining tree:
plot2 <- ggtree(g2NJ, ladderize = TRUE, layout = "unrooted") +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_text_repel(aes(label=c(pops$popName, rep(NA, 10))), force=5, point.padding = unit(1.6, 'lines'), box.padding = unit(0.5, 'lines'), segment.size=0.75, max.iter = 3e3, arrow = arrow(length = unit(0.01, 'npc')), size=3.5, segment.color="gray") +
	ggtitle("Unrooted tree - Fst - (neighbor-joining)") +
	xlim(-0.5, 0.35)

# Save ggtree plots
ggsave(paste0(outpath, "plot_Fst_rooted_tree_all_pops_",today,"ERD.pdf"), plot=plot1, width=6, height=6)
ggsave(paste0(outpath, "plot_Fst_unrooted_tree_all_pops_",today,"ERD.pdf"), plot=plot2, width=6, height=6)








##### Remove the outgroup:
fst2nog <- fst2[-12, -12]

# UPGMA and NJ tree on everything:
treeUPGMAnog <- upgma(as.dist(fst2nog))
treeNJnog <- NJ(as.dist(fst2nog))

##reformat to phylo4 for plotting
# UPGMA
g1Unog = as(treeUPGMAnog, 'phylo4')
dUnog = data.frame(color=sapply(rownames(fst2nog), pop.cols))
rownames(dUnog) = treeUPGMAnog$tip.label
g2Unog = phylo4d(g1Unog, dUnog)

# NJ
g1NJnog <- as(treeNJnog, 'phylo4')
dNJnog <- data.frame(color=sapply(rownames(fst2nog), pop.cols))
rownames(dNJnog) <- treeNJnog$tip.label
g2NJnog <- phylo4d(g1NJnog, dNJnog)

## Plot using ggtree
# Plot rooted UPGMA:
labs <- pops$popName[-12]
treelabs <- paste0("   ", c(labs, rep(NA, length(labels(g2Unog, "internal")))))

plot1 <- ggtree(g2Unog, ladderize = TRUE) +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_tiplab(aes(label=treelabs), size=3.5, color="gray20") +
	ggtitle("Rooted tree - Fst - (UPGMA)") +
	xlim(0, 0.15)
	
# Plot unrooted neighbor-joining tree:
plot2 <- ggtree(g2NJnog, ladderize = TRUE, layout = "unrooted") +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_text_repel(aes(label=c(pops$popName[-12], rep(NA, 9))), force=5, point.padding = unit(1.6, 'lines'), box.padding = unit(0.5, 'lines'), segment.size=0.75, max.iter = 3e3, arrow = arrow(length = unit(0.01, 'npc')), size=3.5, segment.color="gray") +
	ggtitle("Unrooted tree - Fst - (neighbor-joining)") 

# Save ggtree plots
ggsave(paste0(outpath, "plot_Fst_rooted_tree_no_outgroup_",today,"ERD.pdf"), plot=plot1, width=6, height=6)
ggsave(paste0(outpath, "plot_Fst_unrooted_tree_no_outgroup_",today,"ERD.pdf"), plot=plot2, width=6, height=6)



##### Only 2009 samples:
keeps <- grep("2009", pops$year)
fst2009 <- fst2[keeps, keeps]

# UPGMA and NJ tree on everything:
treeUPGMA2009 <- upgma(as.dist(fst2009))
treeNJ2009 <- NJ(as.dist(fst2009))

##reformat to phylo4 for plotting
# UPGMA
g1U2009 = as(treeUPGMA2009, 'phylo4')
dU2009 = data.frame(color=sapply(rownames(fst2009), pop.cols))
rownames(dU2009) = treeUPGMA2009$tip.label
g2U2009 = phylo4d(g1U2009, dU2009)

# NJ
g1NJ2009 <- as(treeNJ2009, 'phylo4')
dNJ2009 <- data.frame(color=sapply(rownames(fst2009), pop.cols))
rownames(dNJ2009) <- treeNJ2009$tip.label
g2NJ2009 <- phylo4d(g1NJ2009, dNJ2009)

## Plot using ggtree
# Plot rooted UPGMA:
labs <- pops$popName[keeps]
treelabs <- paste0("   ", c(labs, rep(NA, length(labels(g2U2009, "internal")))))

plot1 <- ggtree(g2U2009, ladderize = TRUE) +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_tiplab(aes(label=treelabs), size=3.5, color="gray20") +
	ggtitle("Rooted tree - Fst - (UPGMA)") +
	xlim(0, 0.15)
	
# Plot unrooted neighbor-joining tree:
plot2 <- ggtree(g2NJ2009, ladderize = TRUE, layout = "unrooted") +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	geom_text_repel(aes(label=c(pops$popName[keeps], rep(NA, 3))), force=5, point.padding = unit(1.6, 'lines'), box.padding = unit(0.5, 'lines'), segment.size=0.75, max.iter = 3e3, arrow = arrow(length = unit(0.01, 'npc')), size=3.5, segment.color="gray") +
	ggtitle("Unrooted tree - Fst - (neighbor-joining)") 

# Save ggtree plots
ggsave(paste0(outpath, "plot_Fst_rooted_tree_2009_only_",today,"ERD.pdf"), plot=plot1, width=6, height=6)
ggsave(paste0(outpath, "plot_Fst_unrooted_tree_2009_only_",today,"ERD.pdf"), plot=plot2, width=6, height=6)




##### Individual trees:
# Convert individual genetic distances to distance matrix:
dm <- reshape(dists[,-3], direction="wide", idvar="V2", timevar="V1")  # reshape to matrix
dm <- as.matrix(rbind(rep(NA, ncol(dm)), dm)) # add missing row
dm <- cbind(dm, rep(NA, nrow(dm))) # add missing column
colnames(dm) <- gsub("V4.", "", colnames(dm)) # remove V4 from rownames
dm[1,1] <- colnames(dm)[2] # add missing name to rownames
rownames(dm) <- dm[,1] # make the first column the rownames
dm <- dm[,-1] # remove that first column
colnames(dm)[ncol(dm)] <- rownames(dm)[nrow(dm)] # Make the last rowname the last colname
assert(colnames(dm) == rownames(dm)) # Make sure the rownames and colnames match
class(dm) <- "numeric" # convert to numeric matrix
dm[upper.tri(dm)] <- t(dm)[upper.tri(dm)] # Make matrix symmetrical
diag(dm) <- 0 # set self distance to be zero


# UPGMA and NJ tree on everything:
treeUPGMAind <- upgma(as.dist(dm))
treeNJind <- NJ(as.dist(dm))

dm2 <- dm
rownames(dm2) <- paste0("s", rownames(dm2)) # phylo4d does not want sample names starting with a number.
treeUPGMAind <- upgma(as.dist(dm2))
treeNJind <- NJ(as.dist(dm2))


##reformat to phylo4 for plotting
# UPGMA
g1Uind = as(treeUPGMAind, 'phylo4')
indPops <- samples$pop_name[match(rownames(dm), samples$sample_ID)]
dUind = data.frame(color=sapply(indPops, pop.cols))
rownames(dUind) = treeUPGMAind$tip.label
g2Uind = phylo4d(g1Uind, dUind)

# NJ
g1NJind <- as(treeNJind, 'phylo4')
dNJind <- data.frame(color=sapply(indPops, pop.cols))
rownames(dNJind) <- treeNJind$tip.label
g2NJind <- phylo4d(g1NJind, dNJind)


## Plot using ggtree

# for labels:
labs <- pops$popName[match(samples$population[match(colnames(dm2), samples$sample_ID)], pops$population)]
treelabs <- paste0("   ", c(labs, rep(NA, length(labels(g2Uind, "internal")))))

# Plot rooted UPGMA:
plot1 <- ggtree(g2Uind, ladderize = TRUE) +
	geom_tippoint(aes(color=I(color)), shape=16, size=3) +
	geom_tiplab(aes(label=treelabs), size=2, color="gray20") +
	ggtitle("Rooted tree - genetic distance - (UPGMA)") +
	xlim(0, 0.15)
	
# Plot unrooted neighbor-joining tree:
plot2 <- ggtree(g2NJind, ladderize = TRUE, layout = "unrooted") +
	geom_tippoint(aes(color=I(color)), shape=16, size=5) +
	ggtitle("Unrooted tree - genetic distance - (neighbor-joining)") 

# Save ggtree plots
ggsave(paste0(outpath, "plot_individual_genetic_distance_rooted_tree_",today,"ERD.pdf"), plot=plot1, width=6, height=6)
ggsave(paste0(outpath, "plot_individual_genetic_distance_unrooted_tree_",today,"ERD.pdf"), plot=plot2, width=6, height=6)


##### rooted UPGMA genetic distance plot with shapes and black and white:
print("note: not saving UPGMA genetic distance plot in black and white")
# UPGMA
g1Uind = as(treeUPGMAind, 'phylo4')
indPops <- samples$pop_name[match(rownames(dm), samples$sample_ID)]
dUind = data.frame(shape=sapply(indPops, pop.pch))
rownames(dUind) = treeUPGMAind$tip.label
g2Uind = phylo4d(g1Uind, dUind)

# for labels:
labs <- pops$popName2[match(samples$population[match(colnames(dm2), samples$sample_ID)], pops$population)]
treelabs <- data.frame(treelabs = paste0("   ", c(labs, rep(NA, length(labels(g2Uind, "internal"))))))

# Plot rooted UPGMA:
plot1 <- ggtree(g2Uind, ladderize = TRUE) +
	geom_tippoint(aes(shape=I(shape)), size=2.5) +
	geom_tiplab(aes(label=treelabs, shape = I(shape)), size=2, color="gray20") +
	xlim(0, 0.15) 
	
# Save ggtree plot
# ggsave(paste0(outpath, "plot_individual_genetic_distance_rooted_tree_shapes_not_colors_",today,"ERD.pdf"), plot=plot1, width=6, height=6)




print("DONE!")

