awk '!/^#/{print $3"\t"$4}' hap1.collinearity.blocks|paste - - > hap1.collinearity.blocks.oneline
cat hap1.collinearity.blocks.oneline|while read R1 R2 R3 R4;do echo $(grep -m1 $R1 hap1.gff|cut -f1,3)" "$(grep -m1 $R3 hap1.gff|cut -f3);done > block_bed1
cat hap1.collinearity.blocks.oneline|while read R1 R2 R3 R4;do echo $(grep -m1 $R2 hap1.gff|cut -f1,3)" "$(grep -m1 $R4 hap1.gff|cut -f3);done > block_bed2
tr ' ' '\t' < block_bed1 > tmp && mv tmp block_bed1
tr ' ' '\t' < block_bed2 > tmp && mv tmp block_bed2
