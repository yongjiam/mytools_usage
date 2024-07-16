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
awk '{print $1, $4".1", $2, $3}' hap1_primary.gff > hap1.gff

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
## input for circlize.ipynb
https://jokergoo.github.io/circlize_book/book/graphics.html#labels ## circlize
https://bookdown.org/hneth/ds4psy/D-3-apx-colors-basics.html ## basic R colors
```
library(circlize)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(ComplexHeatmap) ## legend

hap1_gemoma_final_annotation.gff
hap1.filteredRepeats.gff

hap1.ideogram

block_bed1
block_bed2
```
