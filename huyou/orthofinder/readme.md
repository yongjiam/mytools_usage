## gemoma output replace seqids with Names
```
awk '/^>/ {for (i=1; i<=NF; i++) if ($i ~ /Name=/) print ">"$i; next} {print}' predicted_proteins.fasta > simple_predicted_proteins.fasta
sed -i 's/Name=//' simple_predicted_proteins.fasta
bioawk -c fastx '{if ($name ~ /\.1/) print ">"$name"\n"$seq}' simple_predicted_proteins.fasta > hap1_primary_pep.fasta ## only use primary trans
```
## run orthofinder with all downloaded species, and hap1, hap2
/data/huyou/orthofinder/selected/OrthoFinder/Results_Jul18/Orthogroups/Orthogroups_SingleCopyOrthologues.txt, Orthogroups.tsv

## extract the single copy gene OG and OG gene matrix
```
## loop single copy SWO gene ID and OG id, grep gene id in SWO.gff, output OG and chromosome
head -n1 Orthogroups.tsv > header
grep -f Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > Orthogroups_SingleCopyOrthologues.tsv
awk '{print $1 "\t" $18}' Orthogroups_SingleCopyOrthologues.tsv|while read R1 R2;do echo $(grep -m1 $R2 SWO.v3.0.gene.model.gff3|cut -f1)" "$R1" "$R2;done > SWO_singlecopy_gene_and_OG.txt

## create OG id file for each chromosome
cut -d ' ' -f1,2 SWO_singlecopy_gene_and_OG.txt|while read R1 R2;do (echo $R2 >> $R1"_OG_id");done

## extract gene id for each chromosome
ls chr*_OG_id|while read R;do (grep -f $R Orthogroups_SingleCopyOrthologues.tsv > $R".tsv");done

## add header row
ls chr*OG_id.tsv|while read R;do (cat header $R > "header_"$R);done

## extract gene id for each chromosome and variety
ls header_chr*.tsv|while read R; do awk -v CHR="query_"$R 'NR==1 {for (i=1; i<=NF; i++) {col[i]=$i;}} NR>1 {for (i=1; i<=NF; i++) {print $i >> CHR"_"col[i]".txt";}}' $R;done

## Hap1.fasta and Hap2.fasta gene ID names are overlapping and need to be modified
sed -i '/^>/ s/HY/H1Y/' Hap1.fasta
sed -i '/^>/ s/HY/H2Y/' Hap2.fasta

## extract the fasta for each chromosome and species
cat *.fasta > merged.pep
makeblastdb -in merged.pep -dbtype prot -parse_seqids -out merged
ls query_header_chr*|while read R;do EN=$(echo $R|sed 's/query_header_//;s/OG_id.tsv_//'); mv $R $EN; blastdbcmd -db ../merged_pep/merged -entry_batch $EN -out $EN".fa";done

## separate the sequences files by chromosome and run orthofinder separately
seq 1 9|while read R;do mkdir "chr"$R; mv "chr"$R*.txt.fa "chr"$R;done
ls --color=never -d */|while read R;do orthofinder -t 30 $R;done

## download species tree for each chromosome
seq 1 9|while read R;do scp -i ~/.ssh/mynimbuskey.pem "ubuntu@146.118.64.65:/data/huyou/orthofinder/merged_chromosomes/chr$R/OrthoFinder/Results_Jul24/Species_Tree/SpeciesTree_rooted.txt" "chr$R-SpeciesTree_rooted.txt";done

## modify taxa name to show species names
sed -i '' 's/.txt//g' chr*SpeciesTree_rooted.txt
sed -i '' -f sed_species_name chr* 
```
## based on orthofinder results, swap chr8 between Hap1 and Hap2, redo earlGrey and gemoma
```
## /media/hhd1/yjia/data/huyou/assembled_genomes/swap_chr8
sed '/^>/ s/_RagTag//' hap1.fasta > hap1_updated.fasta
sed '/^>/ s/_RagTag//' hap2.fasta > hap2_updated.fasta

sed -i 's/^>chr/>H1chr/' hap1_updated.fasta
sed -i 's/^>chr/>H2chr/' hap2_updated.fasta

sed -i 's/^>H1chr8/>H2chr8/' hap1_updated.fasta
sed -i 's/^>H2chr8/>H1chr8/' hap2_updated.fasta

cat hap1_updated.fasta hap2_updated.fasta > merged.fasta
bioawk -c fastx '{if ($name ~ /H1/) print ">"$name"\n"$seq}' merged.fasta  > hap1_updated.fasta
bioawk -c fastx '{if ($name ~ /H2/) print ">"$name"\n"$seq}' merged.fasta  > hap2_updated.fasta

## sort sequence name
mamba install bioconda::seqkit
seqkit sort -n hap1_updated.fasta -o tmp && mv tmp hap1_updated.fasta
seqkit sort -n hap2_updated.fasta -o tmp && mv tmp hap2_updated.fasta
```
