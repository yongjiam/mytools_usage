## CNV heatmap for bHLH gene family
```R
library(pheatmap)

help(pheatmap)
setwd("/Users/yongjia/Desktop/Programs/HMMER/pawsey_HMMER/bHLH_cnv_pheatmap")
#df1 = read.table("pangenome_reference_gene_cnv_npNan.csv",sep=",",header=T)
df1 = read.table("pangenome_reference_gene_cnv_presenceNum_annotated_shortID.csv",sep=",",header=T)

df_num1 = as.matrix(df1[,2:ncol(df1)])
rownames(df_num1) = df1$reference_gene

heatmap1 <- pheatmap(df_num1,
                     na_col = "white",
                     #scale = 'row',
                     main = "barley bHLH gene cnv",
                     border_color = "black",
                     #cluster_cols = FALSE,
                     fontsize = 13,
                     fontsize_row = 8,
                     fontsize_col = 15,
                     angle_col = 45)
save_pheatmap_pdf <- function(x, filename, width=12, height=23) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap1, "test.pdf")
```
