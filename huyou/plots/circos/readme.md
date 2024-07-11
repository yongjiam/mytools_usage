## visulize the annotation and synteny
## install circlize R package
https://github.com/jokergoo/circlize
```
install.packages("circlize")
```
## prepare input file
```
## hap1.ideogram
bioawk -c fastx '{print $name"\t"1"\t"length($seq)}' hap1.fasta|sed "s/_RagTag//"

## gene density, only gene rows,
grep -v "^#" hap1_gemoma_final_annotation.gff > tmp && mv tmp hap1_gemoma_final_annotation.gff
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_final_annotation.gff|cut -d ';' -f1> hap1_primary.gff

## te density, from earlgrey output
hap1.filteredRepeats.gff

## synteny bed file
### gemoma output replace seqids with Names
awk '/^>/ {for (i=1; i<=NF; i++) if ($i ~ /Name=/) print ">"$i; next} {print}' predicted_proteins.fasta > simple_predicted_proteins.fasta
sed -i 's/Name=//' simple_predicted_proteins.fasta
bioawk -c fastx '{if ($name ~ /\.1/) print ">"$name"\n"$seq}' simple_predicted_proteins.fasta > hap1_primary_pep.fasta ## only use primary trans

### whatever reason, some protein ids not present in gff, use the smaller gene ids from gff gene rows
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_final_annotation.gff|cut -d ';' -f1|cut -d "=" -f2 > hap1_mcscan.gff.id
sed -i 's/$/.1/' hap1_mcscan.gff.id
makeblastdb -in hap1_primary_pep.fasta -dbtype prot -parse_seqids -out hap1_primary_pep
blastdbcmd -db hap1_primary_pep -entry_batch hap1_mcscan.gff.id -out hap1.fasta
makeblastdb -in hap1.fasta -dbtype prot -parse_seqids -out hap1

### prepare mcscan input
blastp -db hap -evalue 1e-10 -query hap1.fasta -outfmt 6 -num_threads 30 -max_target_seqs 5 -out hap1.blast
awk '!/^#/{if ($3 == "gene") print $1"\t"$4"\t"$5"\t"$9}' hap1_final_annotation.gff|cut -d ';' -f1> hap1_primary.gff
