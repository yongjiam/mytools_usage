## filter for black white lines
paste white_black.ids white_black.ids  > tmp && mv tmp white_black.ids
plink --bfile 2071_all_bed --maf 0.01 --geno 0.5 --keep white_black.ids --make-bed --out filtered_samples

## plink binary to vcf
plink --bfile filtered_samples --recode vcf --out filtered_samples
## vcf to hapmap
run_pipeline.pl -Xmx2G -importGuess filtered_samples.vcf -export filtered_samples -exportType HapmapDiploid
cut -f2-3 white_black.phenotype.txt > tmp && mv tmp white_black.phenotype.txt

## phenotype data
IID	white-black
YL0413	1
YL0414	1
YL0417	1

## genotype data
