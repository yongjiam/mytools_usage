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


```
## use SWO as reference, extract the OG matrix for each chromosome 1-9

## for each chromosome, extract the sequence ids for each species

## run orthofinder again for each chromosome, get the phylogeny tree
