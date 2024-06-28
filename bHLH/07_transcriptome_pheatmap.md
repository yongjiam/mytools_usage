## plot transcriptome pheatmap
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
