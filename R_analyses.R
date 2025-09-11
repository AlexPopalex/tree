###############################
### Display the phylogenies ###
###############################


# Load phytools
library(phytools)

# Load metadata for each individual
meta <- read.table("full_table_tree_paper_species.csv", h=T, sep=";")

# Add colour information for each individual based on lineage
meta$col <- "grey"
meta$col[which(meta$sp == "rossius")] <- "#2E368F"
meta$col[which(meta$sp == "redtenbacheri")] <- "#0599CE"
meta$col[which(meta$sp == "Vesuve")] <- "#4166AA"
meta$col[which(meta$sp == "grandii")] <- "#BE1E2D"
meta$col[which(meta$sp == "benazzii")] <- "#F7941D"
meta$col[which(meta$sp == "maretimi")] <- "#F9ED32"
meta$col[which(meta$sp == "atticus")] <- "#71BF44"
meta$col[which(meta$sp == "whitei")] <- "#92278F"
meta$col[which(meta$sp == "hybrido1")] <- "#EC008C"
meta$col[which(meta$sp == "hybrido2")] <- "#F47F72"
meta$col[which(meta$sp == "lynceorum")] <- "#8B5E3C"


# Read the trees
# Whole-genome
genome.tree <- read.tree("genomes_rerooted.treefile")

# Maternal (rossius) RADseq
mat.tree.all <- read.tree("maternal_rerooted.treefile")

# Paternal (grandii and atticus) RADseq
pat.tree.all <- read.tree("paternal_rerooted.treefile")


# Remove individuals not in the metadata
# (artificial hybrids used for validation)
mat.tree <- drop.tip(phy = mat.tree.all, tip = mat.tree.all$tip.label[mat.tree.all$tip.label %in% meta$ID == F], trim.internal = T)

# A bit more complex for the paternal one
# because we want to keep the lynceorum with "_as" at the end
pat.tree <- drop.tip(phy = pat.tree.all, tip = pat.tree.all$tip.label[sapply(strsplit(pat.tree.all$tip.label, "_"), `[[`, 1) %in% meta$ID == F], trim.internal = T)


### Figure 2

par(mfrow=c(1,3))

# Maternal RADseq tree; Figure 2A
plot.phylo(mat.tree, show.tip.label = F)
tiplabels(pch=15, col = meta$col[match(mat.tree$tip.label, meta$ID)], offset = 0.0002, cex=0.5)
add.scale.bar()

# Paternal RADseq tree; Figure 2B
plot.phylo(pat.tree, show.tip.label = F)
add.scale.bar()
# We need a trick to put colours on this tree because some individuals have "_as" at the end
tiplabels(pch=15, col = meta$col[match(sapply(strsplit(pat.tree$tip.label, "_"), `[[`, 1), meta$ID)], cex=0.5, offset = 0.005)

# Whole-genome tree; Figure 2C
plot.phylo(genome.tree, show.tip.label = T)
# tiplabels(pch=15, col = meta$col[match(genome.tree$tip.label, meta$ID)])
add.scale.bar()


### Subsets of genome tree (parts of Figure 2C)

rossius_tips <- c("NWhaplome", "SEhaplome", "Blm_rsri", "Blmge_rsri", "Bwi_rsri")
grandii_tips <- c("Blm_gigi", "Blmge_gigi", "Bwi_gigi")
atticus_tips <- c("Blmge_as", "Bas", "Blm_as")

ros_genome.tree <- drop.tip(phy = genome.tree, tip = genome.tree$tip.label[genome.tree$tip.label %in% rossius_tips == F])
gra_genome.tree <- drop.tip(phy = genome.tree, tip = genome.tree$tip.label[genome.tree$tip.label %in% grandii_tips == F])
att_genome.tree <- drop.tip(phy = genome.tree, tip = genome.tree$tip.label[genome.tree$tip.label %in% atticus_tips == F])

plot.phylo(ros_genome.tree, show.tip.label = T)
add.scale.bar()

plot.phylo(gra_genome.tree, show.tip.label = T)
add.scale.bar()

plot.phylo(att_genome.tree, show.tip.label = T)
add.scale.bar()

## Figure display and species models were made with Affinity Design


# Create a new species category for the second lineage of lynceorum
# (a posteriori)
meta$sp2 <- meta$sp
meta$sp2[which(meta$sp == "lynceorum" & meta$pop %in% c("Gela", "SantoPietro"))] <- "lynceorum2"

###############################
### plotting the mtDNA tree ###
###############################

## For Figure S6
mt.tree <- read.newick("mtDNA_tree.newick")
plot.phylo(mt.tree, type = "unrooted")
add.scale.bar()

## Final display and species signs added with Affinity Design


#####################################
### Transition to parthenogenesis ###
#####################################

## In a loop for each individual (listed in windows/list.txt)
## All plots were then scored independently by 3 authors (AB, GL, TS),
## whose independent scores converged to the same conclusions.
## Representative individuals were then displayed in Figure 3C

for (ind in read.table("windows/list.txt", h=F)$V1) {
    # Read per-window coverage of the grandii genome
  gra <- read.table(paste0("windows/", ind, "_covplot_prim_Bgigi_nl.txt"), h=F, col.names = c("window_number", "chr", "start_pos", "end_pos", "coverage"))
  
  # Read per-window coverage of the rossius genome
  ros <- read.table(paste0("windows/", ind, "_covplot_prim_Brsri_corr_nl.txt"), h=F, col.names = c("window_number", "chr", "start_pos", "end_pos", "coverage"))
  
  # Set alternating colours for adjacent chromosomes
    gra$col="#BE1E2D"
  gra$col[which(gra$chr %in% c(2,4,6,8,10,12,14,16))] <- "#ED8E95"
  ros$col="#0599CE"
  ros$col[which(ros$chr %in% c(2,4,6,8,10,12,14,16))] <- "#7BD1F1"
  
  ## Plot coverage 
  # First an empty plot
  png(paste0("windowsplots/", ind, "coverage_gra_ros.png"), width = 1600, height = 400)
  plot(gra$coverage ~ gra$window_number, type="l", xaxt="n", col=NULL, pch=NULL, cex=0, las=1, ylab="Average coverage", xlab="", main=ind, ylim=c(0, max(c(gra$coverage, ros$coverage))))
  # then for each chromosome
  for(chrom in 1:17) {
    # plot the value for grandii
    points(gra$coverage[which(gra$chr == chrom)] ~ gra$window_number[which(gra$chr == chrom)], type="l", xaxt="n", col=gra$col[which(gra$chr == chrom)], lwd=3)
    # and for rossius
    points(ros$coverage[which(ros$chr == chrom)] ~ ros$window_number[which(ros$chr == chrom)], type="l", xaxt="n", col=ros$col[which(ros$chr == chrom)], lwd=3)
  }
  dev.off()
}


###############################
### Transition to triploidy ###
###############################

## Essentially the same thing as above
## but without the loop since it's only one individual
## Used in Figure 3D

# Read per-window coverage of the grandii genome
gra <- read.table("B24_covplot_prim_Bgigi_nl.txt", h=F, col.names = c("window_number", "chr", "start_pos", "end_pos", "coverage"))

# Read per-window coverage of the rossius genome
ros <- read.table("B24_covplot_prim_Brsri_corr_nl.txt", h=F, col.names = c("window_number", "chr", "start_pos", "end_pos", "coverage"))

# Set alternating colours for adjacent chromosomes
gra$col="#BE1E2D"
gra$col[which(gra$chr %in% c(2,4,6,8,10,12,14,16))] <- "#ED8E95"
ros$col="#0599CE"
ros$col[which(ros$chr %in% c(2,4,6,8,10,12,14,16))] <- "#7BD1F1"

## Plot coverage 
# First an empty plot
plot(gra$coverage ~ gra$window_number, type="l", xaxt="n", col=NULL, pch=NULL, cex=0, las=1, ylab="Average coverage", xlab="", main="", ylim=c(0, max(c(gra$coverage, ros$coverage))))
# then for each chromosome
for(chrom in 1:17) {
  # plot the value for grandii
  points(gra$coverage[which(gra$chr == chrom)] ~ gra$window_number[which(gra$chr == chrom)], type="l", xaxt="n", col=gra$col[which(gra$chr == chrom)], lwd=3)
  # and for rossius
  points(ros$coverage[which(ros$chr == chrom)] ~ ros$window_number[which(ros$chr == chrom)], type="l", xaxt="n", col=ros$col[which(ros$chr == chrom)], lwd=3)
}
