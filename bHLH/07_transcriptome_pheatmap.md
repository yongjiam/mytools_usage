## plot morex transcriptome pheatmap
```R
library(pheatmap)
library(ape)
help(pheatmap)
setwd("/Volumes/Elements5T/Programs/Students/Tong_bhlh_transcriptome_heatmap/transcriptome")
df1 = read.table("seded_raw_TPM_log2.csv",sep="\t",header=T) ## raw data


df_num1 = as.matrix(df1[,2:ncol(df1)])


rownames(df_num1) = df1$X

heatmap1 <- pheatmap(df_num1,main = "seded_raw_TPM",
                     border_color = "black",
                     cluster_cols = FALSE,
                     fontsize = 13,
                     fontsize_row = 8,
                     fontsize_col = 15,
                     angle_col = 45)

my_tree <- as.phylo(heatmap1$tree_row)
my_tips <- my_tree$tip.label
write.csv(my_tips, file = "testing_treetips_order.csv")
write.tree(phy=my_tree, file="seded_raw_TPM_log2.pheatmap.newick") # look for the file in your working directory

save_pheatmap_pdf <- function(x, filename, width=12, height=20) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


save_pheatmap_pdf(heatmap1, "seded_raw_TPM_log2.pheatmap.pdf")
```
## plot cophyloplot for gene phylogeny and transcriptome clustering tree
```R
########### cophyloplot to compare 2 trees
#setwd("/Volumes/Elements5T/Programs/Students/Tong_bhlh_transcriptome_heatmap/transcriptome")
setwd("/Volumes/Elements5T/Backups/BHLH/update_analyses")
library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)

library("RColorBrewer")
library(randomcoloR)

##################create a color legend
n <- 23
palette <- distinctColorPalette(n)
display.brewer.pal(n = 8, name = 'RdBu')

# Generate some data
# Add a legend
legends <- c("NA","SF01","SF02","SF03","SF04","SF05","SF07","SF08","SF09","SF10","SF11","SF12","SF13","SF14","SF15","SF16","SF24","SF25","SF26","SF27","SF28","SF31","SF_like")
cols <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99")

pdf("updated_cophyloplot_legends.pdf",width = 12, height = 9)
x<-1:10; y1=x*x; y2=2*y1
plot(x, y1, type="b", pch=19, col="red", xlab="x", ylab="y")
legend(1, 95, legend=legends, fill = cols,cex=0.5)
dev.off()

##########################
#mx <- read.tree("pheatmap_tree_gene_ids.fas.nwk") ## gene family tree
#mx <- read.tree("bHLH.nwk") ## original full tree
#mx <- read.tree("./morexV2_newHLH/seded_merged_dom_pep.fas.nwk") ## newly built full tree
#mx <- read.tree("./morexV2_newHLH/seded_merged_dom_pep_reversed.fas.nwk") ## newly built full tree
#mx <- read.tree("./morexV2_newHLH/seded_reversedID2_tiplabels_in_order_dom_pep.fas.nwk") ## newly built ID2 only tree
mx <- read.tree("./morexV2_newHLH/seded_ID2_tiplabels_in_order_dom_pep.fas.nwk") ## newly built ID2 only tree
my <- read.tree("pheatmap_mean_log2_tree.newick") ## pheatmap cluster tree
#my <- read.tree("test.nwk") ## export pheatmap tree using figtree
#my <- read.tree("pheatmap_mean_log2_tree_no_branch.nwk") ## pheatmap tree without branch length

ID1 <- mx$tip.label
ID2 <- my$tip.label
updated_ID2 <- ID2[grepl("HLH",ID2)]

association <- cbind(ID2, ID2)
association

LINE_COLOR <- readLines("./cophyloplot_cols_inputs.txt")
# use fix(plotCophylo3) to modify cex= to change font size
pdf("ID2only_cophyloplot.pdf", width = 25, height = 25)
cophyloplot(mx, my, assoc = association,
                    length.line = 0,
                    gap = 25,
                    col = LINE_COLOR,
                    edge.width = 8,
                    lwd = 4,
                    cex = 8,
                    font = 1,
                    space = 300)
# Add a legend
legends <- c("NA","SF01","SF02","SF03","SF04","SF05","SF07","SF08","SF09","SF10","SF11","SF12","SF13","SF14","SF15","SF16","SF24","SF25","SF26","SF27","SF28","SF31","SF_like")
cols <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99")
legend("topleft", legend=legends, fill = cols, cex=2)
dev.off()

pdf("updated_cophyloplot_treeonly.pdf", width = 25, height = 30)
cophyloplot(mx, my,
            length.line = 0,
            gap = 25,
            col = LINE_COLOR,
            lwd = 4,
            font = 2,
            space = 200)
# Add a legend
legends <- c("NA","SF01","SF02","SF03","SF04","SF05","SF07","SF08","SF09","SF10","SF11","SF12","SF13","SF14","SF15","SF16","SF24","SF25","SF26","SF27","SF28","SF31","SF_like")
cols <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99")
legend("topleft", legend=legends, fill = cols, cex=1.5)
dev.off()
```
