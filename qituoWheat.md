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
```
cat qituo_gwas/sample_ids.txt |while read R;do (sed "s/SAMPLE/$R/g" bwa.conf > qituo_gwas/$R".conf");done
```
```
while IFS= read -r filename; do [[ ! -e $filename"_sort.bam" ]] && echo "$filename does not exist."; done < sample_ids.txt
```
