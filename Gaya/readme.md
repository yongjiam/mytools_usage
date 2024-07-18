## blast search for query homologous protein sequences
```
## /Volumes/Elements5T/other_species/chickpea_NCBI_genome
# top 5 hits for each sequence
blastp -db protein -evalue 1e-30 -query updated_gaya.fasta -outfmt 6 -num_threads 8 -max_target_seqs 5 -out updated_gaya.blastp

## get protein ID and gene ID match file
cut -f2 gaya.blastp|while read R;do echo $R " "$(grep $R GCF_000331145.1_ASM33114v1_translated_cds.faa|cut -d ' ' -f2);done > protein_gene_id_match
sed -i ''  's/gene=//' protein_gene_id_match
sed -i ''  's/\[//' protein_gene_id_match
sed -i ''  's/\]//' protein_gene_id_match
cat protein_gene_id_match|sort|uniq > tmp && mv tmp protein_gene_id_match_uniq

cat protein_gene_id_match_uniq|while read R1 R2;do grep -m1 $R2 all_genes.gff;done > gaya_gene.gff
paste protein_gene_id_match_uniq gaya_gene.gff > homolog_gene.txt
```
## create chickpea ideogram
```
bioawk -c fastx '{print $name"\t"1"\t"length($seq)}' GCF_000331145.1_ASM33114v1_genomic.fna|head -n8 > chickpea.ideogram

## change chromosome name
cat chickpea.ideogram|cut -f1|awk '{print "s/"$1"/chr"NR"/"}' > sed_chr
sed -i '' -f sed_chr chickpea.ideogram
sed -i '' -f sed_chr homolog_gene.txt
```

## extract snp for target genes
https://cegresources.icrisat.org/cicerseq/?page_id=3444
```
## download snp in hmp format from link above
## install tassel5.0 to convert hmp to vcf
git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git

gzip -d *.gz
sed -i '/^rs#/ s/ /-/g' Cultivated_*.hmp.txt ## accession ID containing spaces
ls *.hmp.txt|while read R;do /data/tools/tassel-5-standalone/run_pipeline.pl -Xmx100g -importGuess $R -export $(echo $R|cut -d '.' -f1)".vcf" -exportType VCF;done

bgzip *.vcf
ls *.vcf.gz|while read R;do bcftools index --threads 30 $R;done
bcftools concat --threads 30 -Oz -o Cultivated_merged.vcf.gz Cultivated_Ca*.vcf.gz


```
