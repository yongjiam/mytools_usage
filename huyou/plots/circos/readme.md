## visulize the annotation and synteny
## install circlize R package
https://github.com/jokergoo/circlize
```
install.packages("circlize")
```
## prepare hap1.ideogram
```
bioawk -c fastx '{print $name"\t"1"\t"length($seq)}' hap1.fasta|sed "s/_RagTag//"|grep -v "chrUn" > hap1.ideogram
```
## prepare gene bed files
```
## only retain gene lines
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_gemoma_final_annotation.gff|cut -d ';' -f1> hap1_primary.gff
sed -i '' 's/_RagTag//g' hap1_primary.gff
```
## prepare te bed, from earlgrey output
```
## run tRNAscan-SE
conda activate bio
srun --export=all -n 1 -c 64   tRNAscan-SE -j hap2tRNA.gff3 --thread 64 hap2.softmasked.fasta

## get output file
hap2.filteredRepeats.gff
```
## prepare other noncoding RNA bed file, from inferal output
https://docs.rfam.org/en/latest/genome-annotation.html
```
## run inferal 
conda activate bio
srun --export=all -n 1 -c 64 cmscan -Z 638.800824 \
	--cut_ga --rfam --nohmmonly --tblout hap2_mcscan.tblout \
	--fmt 2 --cpu 64 --clanin Rfam.clanin Rfam.cm hap2.softmasked.fasta > hap2.cmscan
srun --export=all -n 1 -c 64 grep -v " = " hap2_mcscan.tblout > hap2_mcscan.deoverlapped.tblout

## convert tblout to bed
awk '!/^#/{print $4"\t"$10"\t"$11"\t"$12"\t"$5"\t"$2}' inferalRfam/hap1_mcscan.deoverlapped.tblout > hap1inferal.bed

## remove _RagTag in chromosome ids
sed -i '' 's/_RagTag//g' *.bed

## some lines start > end, swap
awk '{if ($2 > $3) {temp = $2; $2 = $3; $3 = temp} print $0}' hap2inferal.bed > tmp && mv tmp hap2inferal.bed
```
## synteny bed file
```
### gemoma output replace seqids with Names
awk '/^>/ {for (i=1; i<=NF; i++) if ($i ~ /Name=/) print ">"$i; next} {print}' predicted_proteins.fasta > simple_predicted_proteins.fasta
sed -i 's/Name=//' simple_predicted_proteins.fasta
bioawk -c fastx '{if ($name ~ /\.1/) print ">"$name"\n"$seq}' simple_predicted_proteins.fasta > hap1_primary_pep.fasta ## only use primary trans

### whatever reason, some protein ids not present in gff, use the smaller gene ids from gff gene rows
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_final_annotation.gff|cut -d ';' -f1|cut -d "=" -f2 > hap1_mcscan.gff.id
## add ".1"
sed -i 's/$/.1/' hap1_mcscan.gff.id
## extract protein sequences present in gff file
makeblastdb -in hap1_primary_pep.fasta -dbtype prot -parse_seqids -out hap1_primary_pep
blastdbcmd -db hap1_primary_pep -entry_batch hap1_mcscan.gff.id -out hap1.fasta
makeblastdb -in hap1.fasta -dbtype prot -parse_seqids -out hap1

## prepare mcscan input
blastp -db hap -evalue 1e-10 -query hap1.fasta -outfmt 6 -num_threads 30 -max_target_seqs 5 -out hap1.blast
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_final_annotation.gff|cut -d ';' -f1|sed 's/Name=//'> hap1_primary.gff
awk '{print $1, $4".1", $2, $3}' hap1_primary.gff > hap1.gff  ## has to be tab-delimited file

### install and run mcscanx
git clone https://github.com/wyp1125/MCScanX
cd MCScanX
TMPDIR=/media/hhd1/yjia/tools/MCScanX/tmp make
cp MCScanX MCScanX_h duplicate_gene_classifier ../mybin/
MCScanX ./hap1

### create collinear blocks
bash extract_collinear_region.sh > hap1.collinearity.blocks
bash create_block_bed.sh
```

## color ideogram based on ancestor origin determined by kPart
```
## processed kpart results, exported by tab-delimited text files: hap1_seded_100kb_kmer_dist_classified.txt, hap2_seded_100kb_kmer_dist_classified.txt

## separate start and end position
awk 'NR>1{print $1, $3, $NF}' hap1_seded_100kb_kmer_dist_classified.txt|tr '-' '\t' > hap1_kpart.bed
awk 'NR>1{print $1, $3, $NF}' hap2_seded_100kb_kmer_dist_classified.txt|tr '-' '\t' > hap2_kpart.bed

## replace space with tab
sed -i '' 's/ /\t/g' hap1_kpart.bed
sed -i '' 's/ /\t/g' hap2_kpart.bed

## replace ancestor with colors, sed_colors
s/mandarin/orange/
s/pummelo/purple/
s/citron/green/
s/Undetermined/white/

sed -i '' -f sed_colors hap1_kpart.bed
sed -i '' -f sed_colors hap2_kpart.bed

## merge neighbouring regions
bash merged_regions.sh

cat merged_hap1_kpart.bed merged_hap2_kpart.bed > merged_hap1hap2_kpart.bed

```
## input for circlize.ipynb
https://jokergoo.github.io/circlize_book/book/graphics.html#labels ## circlize
https://bookdown.org/hneth/ds4psy/D-3-apx-colors-basics.html ## basic R colors
```R
library(circlize)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(ComplexHeatmap) ## legend
library(stringr)

# Step 2: Read and Process GFF3 Data
# Replace 'gene_annotation.gff3' and 'te_annotation.gff3' with your file paths
gene_gff3 <- import("hap2_updated_gene.gff") ## only gene rows, grep -v "^#" hap1_gemoma_final_annotation.gff > tmp
te_gff3 <- import("hap2.filteredRepeats.gff") ## sed -i '' 's/_RagTag//g' *.gff
trna_gff3 <- import("hap2_updated_tRNA.gff3") ## sed -i '' 's/_RagTag//g' *.gff

# Filter for relevant feature types (e.g., "gene" and "transposable_element")
genes <- gene_gff3[gene_gff3$type == "gene"] 
tes <- te_gff3[te_gff3$type != "Simple_repeat"]
trna <- trna_gff3[trna_gff3$type == "tRNA"] 

gene_bed0 <- genes %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")
te_bed0 <- tes %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")
trna_bed0 <- trna %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")

# Rename columns to match BED format
colnames(gene_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")
colnames(te_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")
colnames(trna_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")

# read inferal annotation
ncRNA_bed0 <- read.table("hap2updated_inferal.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(ncRNA_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")

##remove chrUn
hap2gene_bed <- gene_bed0 %>%
  filter(chrom != "chrUn")
hap2te_bed <- te_bed0 %>%
  filter(chrom != "chrUn")
hap2trna_bed <- trna_bed0 %>%
  filter(chrom != "chrUn")

hap2ncRNA_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn")
hap2sno_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn" & str_detect(name, "sno") )
hap2MIR_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn" & str_detect(name, "MIR") )

## add ".2" to chromosome names
#hap2gene_bed$chrom <- paste0(hap2gene_bed$chrom, ".2")
#hap2te_bed$chrom <- paste0(hap2te_bed$chrom, ".2")
#hap2trna_bed$chrom <- paste0(hap2trna_bed$chrom, ".2")
#hap2sno_bed$chrom <- paste0(hap2sno_bed$chrom, ".2")
#hap2MIR_bed$chrom <- paste0(hap2MIR_bed$chrom, ".2")


# Step 2: Read and Process GFF3 Data
# Replace 'gene_annotation.gff3' and 'te_annotation.gff3' with your file paths
gene_gff3 <- import("hap1_updated_gene.gff") ## only gene rows, grep -v "^#" hap1_gemoma_final_annotation.gff > tmp
te_gff3 <- import("hap1.filteredRepeats.gff") ## sed -i '' 's/_RagTag//g' *.gff
trna_gff3 <- import("hap1_updated_tRNA.gff3") ## sed -i '' 's/_RagTag//g' *.gff


# Filter for relevant feature types (e.g., "gene" and "transposable_element")
genes <- gene_gff3[gene_gff3$type == "gene"] 
tes <- te_gff3[te_gff3$type != "Simple_repeat"]
trna <- trna_gff3[trna_gff3$type == "tRNA"] 

gene_bed0 <- genes %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")
te_bed0 <- tes %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")
trna_bed0 <- trna %>%
  as.data.frame() %>%
  select(seqnames, start, end, strand, score = NULL, name = NULL) %>%
  mutate(score = ".", name = ".")

# Rename columns to match BED format
colnames(gene_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")
colnames(te_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")
colnames(trna_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")

# read inferal annotation
ncRNA_bed0 <- read.table("hap1updated_inferal.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(ncRNA_bed0) <- c("chrom", "chromStart", "chromEnd", "strand", "score", "name")

##remove chrUn
hap1gene_bed <- gene_bed0 %>%
  filter(chrom != "chrUn")
hap1te_bed <- te_bed0 %>%
  filter(chrom != "chrUn")
hap1trna_bed <- trna_bed0 %>%
  filter(chrom != "chrUn")

ncRNA_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn")
hap1sno_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn" & str_detect(name, "sno") )
hap1MIR_bed <- ncRNA_bed0 %>%
  filter(chrom != "chrUn" & str_detect(name, "MIR") )

## add ".2" to chromosome names
#hap1gene_bed$chrom <- paste0(hap1gene_bed$chrom, ".1")
#hap1te_bed$chrom <- paste0(hap1te_bed$chrom, ".1")
#hap1trna_bed$chrom <- paste0(hap1trna_bed$chrom, ".1")
#hap1sno_bed$chrom <- paste0(hap1sno_bed$chrom, ".1")
#hap1MIR_bed$chrom <- paste0(hap1MIR_bed$chrom, ".1")

# Combine dataframes row-wise
gene_bed <- rbind(hap1gene_bed, hap2gene_bed)
te_bed <- rbind(hap1te_bed, hap2te_bed)
trna_bed <- rbind(hap1trna_bed, hap2trna_bed)
sno_bed <- rbind(hap1sno_bed, hap2sno_bed)
MIR_bed <- rbind(hap1MIR_bed, hap2MIR_bed)

# Sort dataframe by 'chrom'
gene_bed_sorted <- gene_bed[order(gene_bed$chrom), ]
te_bed_sorted <- te_bed[order(te_bed$chrom), ]
trna_bed_sorted <- trna_bed[order(trna_bed$chrom), ]
sno_bed_sorted <- sno_bed[order(sno_bed$chrom), ]
MIR_bed_sorted <- MIR_bed[order(MIR_bed$chrom), ]

## Read the custom cytoband file
## bioawk -c fastx '{print $name"\t"1"\t"length($seq)}' hap1.fasta|sed "s/_RagTag//" > hap1.ideogram
custom_cytoband1 <- read.table("hap1_updated.ideogram", header = FALSE, stringsAsFactors = FALSE)
custom_cytoband2 <- read.table("hap2_updated.ideogram", header = FALSE, stringsAsFactors = FALSE)

## Rename columns for clarity
colnames(custom_cytoband1) <- c("chrom", "start", "end")
#custom_cytoband1$chrom <- paste0(custom_cytoband1$chrom, ".1")

colnames(custom_cytoband2) <- c("chrom", "start", "end")
#custom_cytoband2$chrom <- paste0(custom_cytoband2$chrom, ".2")

# merge 
custom_cytoband <- rbind(custom_cytoband1, custom_cytoband2)
# sort
custom_cytoband_sorted <- custom_cytoband[order(custom_cytoband$chrom), ]
custom_cytoband_sorted

hap2bed1 <- read.table("nonChrUn_hap2block_bed1", header = FALSE, stringsAsFactors = FALSE)
hap2bed2 <- read.table("nonChrUn_hap2block_bed2", header = FALSE, stringsAsFactors = FALSE)
## Rename columns for clarity
colnames(hap2bed1) <- c("chrom", "start", "end")
colnames(hap2bed2) <- c("chrom", "start", "end")
#hap2bed1$chrom <- paste0(hap2bed1$chrom, ".2")
#hap2bed2$chrom <- paste0(hap2bed2$chrom, ".2")

hap1bed1 <- read.table("nonChrUn_hap1block_bed1", header = FALSE, stringsAsFactors = FALSE)
hap1bed2 <- read.table("nonChrUn_hap1block_bed2", header = FALSE, stringsAsFactors = FALSE)
## Rename columns for clarity
colnames(hap1bed1) <- c("chrom", "start", "end")
colnames(hap1bed2) <- c("chrom", "start", "end")
#hap1bed1$chrom <- paste0(hap1bed1$chrom, ".1")
#hap1bed2$chrom <- paste0(hap1bed2$chrom, ".1")

# merge
bed1 <- rbind(hap1bed1,hap2bed1)
bed2 <- rbind(hap1bed2,hap2bed2)

# Step 1: Get unique names and generate colors
unique_names1 <- unique(hap1bed1$chrom)
colors1 <- rainbow(length(unique_names1))  # Generate colors

unique_names2 <- unique(hap2bed1$chrom)
colors2 <- rainbow(length(unique_names2))  # Generate colors

# Create a dictionary mapping names to colors
name_color_dict1 <- setNames(colors1, unique_names1)
name_color_dict2 <- setNames(colors2, unique_names2)
name_color_dict <- c(name_color_dict1,name_color_dict2)

# Step 2: Perform replacement and output color list
bed1$color <- name_color_dict[bed1$chrom]  # Map names to colors
color_list <- bed1$color  # Output color list
head(bed1)
```
## sort by chromosome
```
png("merged_circos_sort_chrom.png", width = 1000, height = 1000)
# Initialize the circos plot with the custom ideogram
# Set the start angle
circos.par(start.degree = 90, gap.after = 0.5, track.margin = c(0,0))
circos.initializeWithIdeogram(custom_cytoband_sorted, plotType = NULL)

circos.track(ylim = c(0, 1), track.height = 0.05, panel.fun = function(x, y) {
    breaks = seq(0, 1e9, by = 5e6)
    circos.genomicAxis(major.at = breaks, labels = paste0(breaks/1e6, ""), labels.cex = 1)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, cex=1.5)
})

# Plot gene density
circos.genomicDensity(gene_bed_sorted, col = "blue", track.height = 0.1, window.size = 1e6, 
                      bg.lwd = 0.3 )

# Plot te density
circos.genomicDensity(te_bed_sorted, col = "red", track.height = 0.1, window.size = 1e6,
                    bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(trna_bed_sorted, col = "yellowgreen", track.height = 0.1, window.size = 1e6,
                     bg.lwd = 0.3)
# circos.genomicDensity(tRNA_bed, col = "yellowgreen", track.height = 0.1, window.size = 1e6)

# Plot trna density
# circos.genomicDensity(rRNA_bed, col = "yellow2", track.height = 0.1, window.size = 1e6,
#                      bg.col = "whitesmoke", bg.lwd = 0.5)

# Plot trna density
circos.genomicDensity(sno_bed_sorted, col = "tan2", track.height = 0.1, window.size = 1e6,
                     bg.lwd = 0.3)

# Plot microRNA/MIR density
circos.genomicDensity(MIR_bed_sorted, col = "darkorange", track.height = 0.1, window.size = 1e6,
                     bg.lwd = 0.3)

# plot mcscan synteny
circos.genomicLink(bed1, bed2, col = color_list, border = NA)

# text(0, 0, "Hap2", cex = 2)

circos.clear()
# Close the PNG device
dev.off()
```
## unsorted chromosome
```
# Add a new column with the last four characters
custom_cytoband <- custom_cytoband_sorted %>%
    mutate(last_four = substr(chrom, nchar(chrom) - 3, nchar(chrom)))

# Sort the data frame based on the new column
custom_cytoband <- custom_cytoband %>%
    arrange(last_four) %>%
    select(-last_four)  # Remove the helper column if not needed

# View the sorted data
print(custom_cytoband)

png("merged_circos_nosort_chrom.png", width = 1000, height = 1000)
# Initialize the circos plot with the custom ideogram
# Set the start angle
circos.par(start.degree = 90, gap.after = 0.5, track.margin = c(0,0))
circos.initializeWithIdeogram(custom_cytoband, plotType = NULL, sort.chr = FALSE)

circos.track(ylim = c(0, 1), track.height = 0.05, panel.fun = function(x, y) {
    breaks = seq(0, 1e9, by = 5e6)
    circos.genomicAxis(major.at = breaks, labels = paste0(breaks/1e6, ""), labels.cex = 1)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, cex=1.5)
})

# Plot gene density
circos.genomicDensity(gene_bed_sorted, col = "blue", track.height = 0.1, window.size = 1e6, 
                      bg.lwd = 0.3 )

# Plot te density
circos.genomicDensity(te_bed_sorted, col = "red", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(trna_bed_sorted, col = "yellowgreen", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)
# circos.genomicDensity(tRNA_bed, col = "yellowgreen", track.height = 0.1, window.size = 1e6)

# Plot trna density
# circos.genomicDensity(rRNA_bed, col = "yellow2", track.height = 0.1, window.size = 1e6,
#                      bg.col = "whitesmoke", bg.lwd = 0.5)

# Plot trna density
circos.genomicDensity(sno_bed_sorted, col = "tan2", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(MIR_bed_sorted, col = "darkorange", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 1)

# plot mcscan synteny
circos.genomicLink(bed1, bed2, col = color_list, border=NA) ## border = NA

# text(0, 0, "Hap2", cex = 2)

circos.clear()
# Close the PNG device
dev.off()
```
## display ancestors
```
duplicated_colors1 <- rep(colors1, each = 2)
duplicated_colors1

custom_cytoband$color = duplicated_colors1
custom_cytoband

## read color bed data
bed_data <- read.table("merged_hap1hap2_kpart.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("chrom", "start", "end", "color")

## sort by chromosomes
# Add a new column with the last four characters
bed_data <- bed_data %>%
    mutate(last_four = substr(chrom, nchar(chrom) - 3, nchar(chrom)))

# Sort the data frame based on the new column
bed_data <- bed_data %>%
    arrange(last_four) %>%
    select(-last_four)  # Remove the helper column if not needed

# View the sorted data
head(bed_data)

pdf("merged_circos_nosort_chrom_ancestors.pdf", width = 1000, height = 1000)
# Set circos parameters
circos.par(start.degree = 90, gap.after = 0.5, track.margin = c(0, 0))

# Initialize with the ideogram
circos.initializeWithIdeogram(bed_data, plotType = NULL, sort.chr = FALSE)

# Plot the ideogram with specified colors and add the genomic axis
circos.genomicTrackPlotRegion(bed_data, 
                              ylim = c(0, 1), 
                              panel.fun = function(region, value, ...) {
                                  circos.genomicRect(region, value, col = value$color, border = NA, ...)
                                  # Add the genomic axis
                                  breaks = seq(0, max(bed_data$end), by = 1e7)
                                  circos.genomicAxis(major.at = breaks, labels = paste0(breaks / 1e6, ""), labels.cex = 1)
                                  circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, cex=1.5)
                              },
                              track.height = 0.05)

# Plot gene density
circos.genomicDensity(gene_bed_sorted, col = "blue", track.height = 0.1, window.size = 1e6, 
                      bg.lwd = 0.3 )

# Plot te density
circos.genomicDensity(te_bed_sorted, col = "red", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(trna_bed_sorted, col = "yellowgreen", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)
# circos.genomicDensity(tRNA_bed, col = "yellowgreen", track.height = 0.1, window.size = 1e6)

# Plot trna density
# circos.genomicDensity(rRNA_bed, col = "yellow2", track.height = 0.1, window.size = 1e6,
#                      bg.col = "whitesmoke", bg.lwd = 0.5)

# Plot trna density
circos.genomicDensity(sno_bed_sorted, col = "tan2", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(MIR_bed_sorted, col = "darkorange", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 1)

# plot mcscan synteny
circos.genomicLink(bed1, bed2, col = color_list, border=NA) ## border = NA

# text(0, 0, "Hap2", cex = 2)

circos.clear()
# Close the PNG device
dev.off()
```
## move chrom label outside
```
svg("merged_circos_nosort_chrom_ancestors_updated.svg", width = 1000, height = 1000)
# Set circos parameters
circos.par(start.degree = 90, gap.after = 0.5, track.margin = c(0, 0))

# Initialize with the ideogram
circos.initializeWithIdeogram(bed_data, plotType = NULL, sort.chr = FALSE)

# Plot the ideogram with specified colors and add the genomic axis
circos.genomicTrackPlotRegion(bed_data, 
                              ylim = c(0, 1), 
                              panel.fun = function(region, value, ...) {
                                  circos.genomicRect(region, value, col = value$color, border = NA, ...)
                                  # Add the genomic axis
                                  breaks = seq(0, max(bed_data$end), by = 1e7)
                                  circos.genomicAxis(major.at = breaks, labels = paste0(breaks / 1e6, ""), labels.cex = 1)
                                  circos.text(CELL_META$xcenter, 3, CELL_META$sector.index, cex=1.2)
                              },
                              track.height = 0.06)

# Plot gene density
circos.genomicDensity(gene_bed_sorted, col = "blue", track.height = 0.1, window.size = 1e6, 
                      bg.lwd = 0.3 )

# Plot te density
circos.genomicDensity(te_bed_sorted, col = "red", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(trna_bed_sorted, col = "yellowgreen", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)
# circos.genomicDensity(tRNA_bed, col = "yellowgreen", track.height = 0.1, window.size = 1e6)

# Plot trna density
# circos.genomicDensity(rRNA_bed, col = "yellow2", track.height = 0.1, window.size = 1e6,
#                      bg.col = "whitesmoke", bg.lwd = 0.5)

# Plot trna density
circos.genomicDensity(sno_bed_sorted, col = "tan2", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 0.3)

# Plot trna density
circos.genomicDensity(MIR_bed_sorted, col = "plum4", track.height = 0.1, window.size = 1e6,
                      bg.lwd = 1)

# plot mcscan synteny
circos.genomicLink(bed1, bed2, col = color_list, border=NA) ## border = NA

# text(0, 0, "Hap2", cex = 2)

circos.clear()
# Close the PNG device
dev.off()
```
