## generate kaks gene pairs, plot kaks output
https://github.com/yongjiam/mytools_usage/blob/main/bHLH/Filter_pan_genes-Copy1.ipynb

## plot kaks with phylogeny in R using ggtree
```R
setwd("/Volumes/Elements5T/Programs/Students/Tong_bhlh_transcriptome_heatmap/")
library(ape)
library(ggtree)
library(ggplot2)
library(treeio)

### kaks_data_plotting
p <- read.newick("sed_Hv2_bHLH_gene_id_domseq.pep.fas.treefile")

KAKS <- read.csv("sed_tong_bhlh_pan_gene_pairs_updated.kaks.mean.treeSeq.csv", sep=",", stringsAsFactor=F)

KA <- read.csv("Ka.csv",sep = "\t")

pl <- ggtree(p, size = 1, branch.length = 'none') 
pl <- pl +
  ggtree::geom_facet(panel = "Ka", data = KAKS, geom = ggstance::geom_barh,
                     aes(x = Ka),
                     stat = "identity", width = .3) +
  ggtree::geom_facet(panel = "Ks", data = KAKS, geom = ggstance::geom_barh,
                     aes(x = Ks),
                     stat = "identity", width = .3) +
  ggtree::geom_facet(panel = "Ka/Ks", data = KAKS, geom = ggstance::geom_barh,
                     aes(x = Ka.Ks, color = Ka.Ks),
                     stat = "identity", width = .3) +
  scale_colour_gradient2(low = "white", high = "red") +
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  #scale_fill_gradient(low = "white", high = "red", na.value = NA) +
  geom_tiplab(size=2.5) +
  ggtree::xlim_tree(40) +
  theme_tree2() +
  theme(strip.text.x = element_text(size = 20, colour = "black", angle = 0))
  
ggsave("sed_kaks_tree2.pdf", plot = pl, width = 15, height = 15)
```
