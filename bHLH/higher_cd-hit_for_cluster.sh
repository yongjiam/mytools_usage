#!/bin/bash

echo "input cluster name"
read Cname
Cname0=$(echo $Cname | sed 's/ //') # remove space in cluster name

echo "input gene number"
read Gnum

echo "input cd-hit identify threshold"
read Cidnt

grep -A $Gnum "$Cname" 20vars_HLH_merged_30aa_uniqid_list.cd-hit.c95.clstr | grep -v "^>" | awk '{print substr($3,2,length($3)-4)}' > $Cname0"_gene_id.txt"

#get gene sequences
cat $Cname0"_gene_id.txt" | while read R1;do grep -A1 $R1 20vars_HLH_merged_30aa_uniqid_list_newName.fasta;done > $Cname0"_gene_id.fasta"

#higher cd-hit
cd-hit -i $Cname0"_gene_id.fasta" -o $Cname0"_gene_id.cd-hit" -d 0 -c $Cidnt
