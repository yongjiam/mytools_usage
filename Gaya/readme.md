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
