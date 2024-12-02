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

awk '{print $3"\t"$6"\t"$7}' homolog_gene_less.txt |sort > tmp && mv tmp gaya_genes.bed
sed -i '' 's/ca/CA/' gaya_genes.bed
bcftools view -R gaya_gene.bed Cultivated_merged.vcf.gz -o gaya_gene_snp.vcf
```
## annotate snp 
```
## produce protein.fa and cds.fa from genome and gff
## /Volumes/Elements5T/other_species/chickpea_NCBI_genome
gffread  -x cds.fa -y protein.fa -g GCF_000331145.1_ASM33114v1_genomic.fna GCF_000331145.1_ASM33114v1_genomic.gff
scp -i ~/.ssh/mynimbuskey.pem cds.fa protein.fa GCF_000331145.1_ASM33114v1_genomic.fna GCF_000331145.1_ASM33114v1_genomic.gff ubuntu@146.118.64.65:/data/tools/snpEff/data/chickpea/
mv GCF_000331145.1_ASM33114v1_genomic.fna sequences.fa
mv GCF_000331145.1_ASM33114v1_genomic.gff genes.gff

## build snpeff database
## /data/tools/snpEff/
echo "chickpea.genome : chickpea" >> snpEff.config
java -jar snpEff.jar build -gff3 -v chickpea
 == Protein check:	chickpea	OK: 35574	Not found: 36316	Errors: 105	Error percentage: 0.2942907592701589%

##
awk -F "/" '{print $1"/"$3"/"$2"/"}' sed_chr > sed_chr_reverse
sed -i -f sed_chr_reverse gaya_gene_snp.vcf
java -Xmx8g -jar /data/tools/snpEff/snpEff.jar chickpea gaya_gene_snp.vcf > gaya_gene_snp.annotated.vcf
sed -i '/^##contig=/d' gaya_gene_snp.annotated.vcf
```
### chromosome rename
```
## sed_chr
s/NC_021160.1/ca1/
s/NC_021161.1/ca2/
s/NC_021162.1/ca3/
s/NC_021163.1/ca4/
s/NC_021164.1/ca5/
s/NC_021165.1/ca6/
s/NC_021166.1/ca7/
s/NC_021167.1/ca8/
```

## indel variant calling
```
## data source
https://www.nature.com/articles/s41586-021-04066-1#data-availability
PRJNA657888
41586_2021_4066_MOESM3_ESM.xlsx ## passport, /Users/yongjia/Desktop/workstation/Students/Gayathri/pangenomes

/scratch/pawsey0399/yjia/chickpea
## passport, get the 4-7 columns uniq
S.No	Genotype	Alternate name	Cicer species	Market class	Biological status	Geographic region	Latitude	Longitude	Elevation	DoI
1	Annigiri 1	-	Cicer arietinum	Desi	Cultivar	South Asia	n.a.	n.a.	n.a.	10.18730/MFSRZ
2	Avorodhi	-	Cicer arietinum	Desi	Cultivar	South Asia	n.a.	n.a.	n.a.	n.a.
3	BG 1003	-	Cicer arietinum	Kabuli	Cultivar	South Asia	n.a.	n.a.	n.a.	n.a.
4	BG 1053	-	Cicer arietinum	Kabuli	Cultivar	South Asia	n.a.	n.a.	n.a.	n.a.
5	ICC 19665	Blanco Lechoso	Cicer arietinum	Kabuli	Breeding line	Mediterranean	37.42	-6.11028	152	10.18730/MXETZ
6	ICC 4935	C 235	Cicer arietinum	Desi	Cultivar	South Asia	n.a.	n.a.	201	10.18730/MF2GV
7	CDC 512-51	-	Cicer arietinum	Desi	Breeding line	Americas	n.a.	n.a.	n.a.	n.a.
8	Chaffa	Chaffa	Cicer arietinum	Desi	Cultivar	South Asia	n.a.	n.a.	n.a.	10.18730/MF2FT![image](https://github.com/user-attachments/assets/65019bf9-913b-48ee-985f-6aba1cbe84ba)

## get uniq combinations
cut -f4-7 passport.tsv |sort|uniq > sort_uniq ## 77 in total, get 3 each for variant calling
sed -i "s/\t/=/g;s/ /_/g" sort_uniq

cat sort_uniq |while read R;do grep -m3 $R passport.tsv ;done > selected_passport_grep3
sed "s/\t/=/g;s/ /_/g" filereport_read_run_PRJNA657888_tsv.txt > tmp
cut -d '=' -f2 selected_passport_grep3|while read R;do grep "="$R"=" tmp;done > selected_passport_grep3_tsv.txt
sed -i "s/=/\t/g" selected_passport_grep3_tsv.txt
awk -F "\t" '{print $7}' selected_passport_grep3_tsv.txt|sed "s/\;/\\n/" > selected189_ftp_links

## download fastq
sed 's/^/wget /' selected189_ftp_links > selected189_ftp_links.sh
cat selected189_ftp_links.sh|while read R;do echo $R > $(echo $R|cut -d '/' -f7|cut -d '.' -f1)".sh";done

#### template.conf
#!/bin/bash --login

#SBATCH --job-name=SAMPLE
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

conda activate base
srun --export=all -n 1 -c 1 bash SAMPLE

ls --color=never SRR*.sh|while read R;do (sed "s/SAMPLE/$R/" template.conf > $R".conf");done
ls --color=never SRR*.conf|while read R;do sbatch $R;done
```
