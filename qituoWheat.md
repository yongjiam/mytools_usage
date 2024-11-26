## Fielder ref genome:
https://shigen.nig.ac.jp/wheat/komugi/genome/download.jsp
![image](https://github.com/yongjiam/mytools_usage/assets/88641886/c9eae711-1edc-4707-bbd7-19769a3cbf64)

## genome index
```
bwa index 201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz
samtools faidx 201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta
PJAR=/scratch/pawsey0399/yjia/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar
srun --export=all -n 1 -c 64 java -jar $PJAR CreateSequenceDictionary \
   R=201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz \
   O=201216_Fielder_pseudomolecules_V1+unanchored_contigs.dict
```
## mapping and variant calling
```
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
module load gatk4/4.2.5.0--hdfd78af_0
GENOME=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz
export REF=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta
srun --export=all -n 1 -c 32  bwa mem -t 32 -R '@RG\tID:SAMPLE\tSM:SAMPLE\tPL:ILLUMINA' $GENOME SAMPLE_good_1.fq.gz SAMPLE_good_2.fastq.gz | samtools view -@ 16 -Sb - | samtools sort -@ 16 -o SAMPLE_sort.bam
srun --export=all -n 1 -c 32 samtools index -@ 32 -c SAMPLE_sort.bam
srun --export=all -n 1 -c 2 gatk HaplotypeCaller -R $REF -I SAMPLE_sort.bam -ERC GVCF -O SAMPLE_sort.bam.g.vcf
srun --export=all -n 1 -c 2 bgzip -@ 2 SAMPLE_sort.bam.g.vcf ## compress gvcf file
```
## mark duplicates
```
## template.conf
srun --export=all -n 1 -c 10 gatk MarkDuplicates -I sample.bam -M sample.metrics.txt -O mark_sample.bam
srun --export=all -n 1 -c 10 samtools index -@ 10 sample.bam

ls --color=never *.bam|paste -|cut -d '.' -f1|while read R;do (sed "s/sample/$R/g" template.conf > $R".conf");done

```
```
cat qituo_gwas/sample_ids.txt |while read R;do (sed "s/SAMPLE/$R/g" bwa.conf > qituo_gwas/$R".conf");done
```
```
while IFS= read -r filename; do [[ ! -e $filename"_sort.bam" ]] && echo "$filename does not exist."; done < sample_ids.txt
```
```
## merge gvcf file
glnexus_cli --config gatk -m 230 --threads 42 gvcfs/*.g.vcf.gz > wheat_gwas.bcf
bcftools view --threads 42 wheat_gwas.bcf | bgzip -@ 15 -c > wheat_gwas.vcf.gz
bcftools index wheat_gwas.vcf.gz --threads 42
```
## filtration and imputation
```
## get header and first n lines
bcftools view -h wheat_gwas.vcf.gz
bcftools view -H wheat_gwas.vcf.gz|hean -n 10

## separate snp and indel
bcftools view -v snps wheat_gwas.vcf.gz -Oz -o wheat_gwas_snps.vcf.gz --threads 42
bcftools index wheat_gwas_snps.vcf.gz --threads 42
bcftools view -v indels wheat_gwas.vcf.gz -Oz -o wheat_gwas_indels.vcf.gz --threads 42
bcftools index wheat_gwas_indels.vcf.gz --threads 42

## count total, multi allelic sites, and remove multiallelic sites
bcftools view -H wheat_gwas.vcf.gz | wc -l
bcftools view -m2 -M2 -H input.vcf | wc -l
vcftools --gzvcf wheat_gwas_snps.vcf.gz --max-alleles 2 --recode --out wheat_gwas_snps_biallelic.vcf.gz
bgzip -@ 42 -c wheat_gwas_snps_biallelic.vcf.recode.vcf > wheat_gwas_snps_biallelic.vcf.gz
bcftools index --threads 42 wheat_gwas_snps_biallelic.vcf.gz

## filtration
bcftools filter -i 'QUAL > 20 && FMT/DP >= 10 && FMT/GQ >= 20' -S . wheat_gwas_snps_biallelic.vcf.gz -Oz --threads 42 -o bcftools_filtered_wheat_gwas_snps_biallelic.vcf.gz
bcftools index --threads 42 bcftools_filtered_wheat_gwas_snps_biallelic.vcf.gz
vcftools --gzvcf bcftools_filtered_wheat_gwas_snps_biallelic.vcf.gz --maf 0.01 --max-missing 0.01 --recode --out bcftools_vcftools_filtered_wheat_gwas_snps_biallelic

vcftools --gzvcf wheat_gwas_snps_biallelic.vcf.gz --minQ 30 --minDP 10 --minGQ 20 --maf 0.01 --max-missing 0.05 --recode --out filtered_wheat_gwas_snps_biallelic
## extract sample names and count
bcftools query -l wheat_gwas.vcf.gz
## remove one sample
bcftools view -s ^002FSD01 --threads 40 -o updated_filtered_wheat_gwas_snps_biallelic.recode.vcf.gz -Oz filtered_wheat_gwas_snps_biallelic.recode.vcf.gz
vcftools --gzvcf updated_filtered_wheat_gwas_snps_biallelic.recode.vcf.gz --minQ 30 --minDP 10 --minGQ 20 --maf 0.01 --max-missing 0.05 --recode --out final_snp
bgzip final_snp.recode.vcf
bcftools index --threads 42 final_snp.recode.vcf.gz

## remove crossing lines
dffinal['VCFsampleID'].dropna().to_csv('nonCrossing.txt', index=False, header=False)
bcftools view --threads 30 -S nonCrossing.txt -Oz -o nonCrossing_snp final_snp.recode.vcf.gz

```
## GWAS with gapit
```
/Users/yongjia/Desktop/workstation/Students/qituo_wheat/GWAS
## vcf to hapmap
run_pipeline.pl -Xmx2G -importGuess nonCrossing_final.recode.vcf -export nonCrossing_final -exportType HapmapDiploid

## kinship
plink --vcf nonCrossing_final.recode.vcf --make-bed --allow-extra-chr --double-id --vcf-half-call missing --out plink

## run gapit
conda activate gapit
sudo jupyter notebook --allow-root &
```
