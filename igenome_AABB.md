## RNAseq and exome for variant calling against wheatCS
https://expert.cheekyscientist.com/how-to-do-variant-calling-from-rnaseq-ngs-data/
## indexing
```bash
## extract AABB
samtools faidx Triticum_aestivum.IWGSC.dna_rm.toplevel.fa
samtools faidx Triticum_aestivum.IWGSC.dna_rm.toplevel.fa 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B Un > wheatCS_AABB.fasta
samtools faidx wheatCS_AABB.fasta

## picard dictionary
PJAR=/data/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar
java -jar $PJAR CreateSequenceDictionary R=Triticum_aestivum.IWGSC.dna_rm.toplevel.fa O=Triticum_aestivum.IWGSC.dna_rm.toplevel.fa.dict
java -jar $PJAR CreateSequenceDictionary R=wheatCS_AABB.fasta O=wheatCS_AABB.fasta.dict
```
## STAR mapping
```bash
STAR --runThreadN 30 --genomeDir /data/igenome/wheatCS/genome_index \
        --readFilesIn trimmed.TR.R1.fq.gz trimmed.TR.R2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix igenomeRNA \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
```
## prepare bam
```bash
## mark duplicates
gatk MarkDuplicatesSpark -I $BAM -O "dedup_"$BAM ## mark duplicates

## add RG tags
samtools addreplacerg -O BAM -@ 15 -o $OUT_BAM -r '@RG\tID:CRR289962\tSM:CRR289962\tPL:ILLUMINA' IN_BAM
## index
samtools index -c -@ 15 $BAM ## index bam, use -c for long chromosome
```
## all in one
```bash
#!/bin/bash --login

#SBATCH --job-name=gatk
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=120G
#SBATCH --export=NONE

module load gatk4/4.2.5.0--hdfd78af_0
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
BAM="igenomeRNA.bam"
REF="wheatCS_AABB.fasta"
srun --export=all -n 1 -c 30  gatk MarkDuplicates -I $BAM -O "dedup_"$BAM -M marked_dup_metrics.txt
srun --export=all -n 1 -c 30 samtools addreplacerg -O BAM -@ 30 -o "RG_dedup_"$BAM -r '@RG\tID:HT621\tSM:HT621\tPL:ILLUMINA' "dedup_"$BAM
srun --export=all -n 1 -c 30 samtools index -c -@ 30 "RG_dedup_"$BAM
srun --export=all -n 1 -c 30 gatk SplitNCigarReads -R $REF -I "RG_dedup_"$BAM -O "split_RG_dedup_"$BAM ## note: the dictionary file should be Triticum_aestivum.IWGSC.dna_rm.toplevel.dict
srun --export=all -n 1 -c 30 samtools index -c -@ 30 "split_RG_dedup_"$BAM
srun --export=all -n 1 -c 30 gatk HaplotypeCaller -R $REF -I "split_RG_dedup_"$BAM -ERC GVCF -O "split_RG_dedup_"$BAM".g.vcf"
```
## merge gvcf and genotyping using glnexus
```bash
#!/bin/bash --login

#SBATCH --job-name=gln
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=980G
#SBATCH --export=NONE
module load bcftools/1.15--haf5b3da_0 ## bgzip

srun --export=all -n 1 -c 128 glnexus_cli --config gatk -m 980 --threads 128 glnexus/*.g.vcf.gz > igenome.bcf
srun --export=all -n 1 -c 128 bcftools view --threads 128 igenome.bcf | bgzip -@ 24 -c > igenome.vcf.gz
srun --export=all -n 1 -c 128 bcftools index --threads 128 -c igenome.vcf.gz
```
## filter snps for single copy gene regions
```bash
## sort and index
bcftools sort igenome.vcf.gz -Oz -o sorted_igenome.vcf.gz -T .
bcftools index -c --threads 20 sorted_igenome.vcf.gz

## intersect single copy genes
bedtools intersect -a sorted_igenome.vcf.gz -b updated_wheatCS_AB_single.bed > updated_wheatCS_AB_single.bed.vcf
bgzip updated_wheatCS_AB_single.bed.vcf
bcftools index -c updated_wheatCS_AB_single.bed.vcf.gz

## filter snps
bcftools view -v snps updated_wheatCS_AB_single.bed.vcf.gz > snp.vcf

## add ref genotype to vcf
awk '{print $0 "\t" "0/0"}' snp.vcf > snp_addRef.vcf

## add header
zcat sorted_igenome.vcf.gz | grep "^#" > header.vcf
vi header.vcf ## add wheatCS_AABB
cat header.vcf snp_addRef.vcf > tmp && mv tmp snp_addRef.vcf

## vcf to fasta and phylogeny
python /data/tools/vcf2phylip/vcf2phylip.py -i snp_addRef.vcf -f

python /data/tools/vcf2phylip/vcf2phylip.py -i snp_addRef.vcf -f -r --output-prefix resolve-IUPAC-min30 -m 30
iqtree -s resolve-IUPAC-min30.min30.fasta -B 1000 -T 30

python /data/tools/vcf2phylip/vcf2phylip.py -i snp_addRef.vcf -f -r --output-prefix resolve-IUPAC-wheatCS -m 30 -o wheatCS_AABB ##iqtree may use first taxon as outgroup
iqtree -s resolve-IUPAC-wheatCS.min30.fasta -B 1000 -T 30

```
