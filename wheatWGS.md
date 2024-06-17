## index genome
export REF=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta
PJAR=/scratch/pawsey0399/yjia/tools/miniconda3/envs/nf-env/share/picard-2.18.29-0/picard.jar

srun --export=all -n 1 -c 64 samtools faidx $REF
srun --export=all -n 1 -c 64 java -jar $PJAR CreateSequenceDictionary \
   R=201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz \
   O=201216_Fielder_pseudomolecules_V1+unanchored_contigs.dict
## add sample information to bam if not added during alignment
srun --export=all -n 1 -c 128 samtools addreplacerg -O BAM -@ 128 -o updated_CRR289962_sort.bam  -r '@RG\tID:CRR289962\tSM:CRR289962\tPL:ILLUMINA' CRR289962_sort.bam
srun --export=all -n 1 -c 128 samtools index -@ 128 updated_CRR289962_sort.bam

## wheat WGS variant calling

```bash
# bwa.conf
#!/bin/bash --login

#SBATCH --job-name=bwa
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

#conda activate nf-env
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
module load gatk4/4.2.5.0--hdfd78af_0
GENOME=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta.gz
export REF=/scratch/pawsey0399/yjia/wheat/WGS/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta
srun --export=all -n 1 -c 64  bwa mem -t 64 -R '@RG\tID:SAMPLE\tSM:SAMPLE\tPL:ILLUMINA' $GENOME SAMPLE_f1.fastq.gz SAMPLE_r2.fastq.gz | samtools view -@ 28 -Sb - | samtools sort -@ 28 -o SAMPLE_sort.bam
srun --export=all -n 1 -c 64 samtools index -@ 64 -c SAMPLE_sort.bam
#srun --export=all -n 1 -c 10 gatk HaplotypeCaller -R $REF -I SAMPLE_sort.bam -ERC GVCF -O SAMPLE_2D_sort.bam.g.vcf -L chr2D
#srun --export=all -n 1 -c 10 bgzip -@ 10 SAMPLE_2D_sort.bam.g.vcf ## compress gvcf file
```
```bash
cat third_sample_ids|while read R;do (sed "s/SAMPLE/$R/g" bwa.conf > ./third_fastq/$R".conf");done
```
## install glnexus
https://github.com/dnanexus-rnd/GLnexus
```bash
# Clone repo
git clone https://github.com/dnanexus-rnd/GLnexus.git
cd GLnexus
git checkout vX.Y.Z  # optional, check out desired revision
```
```bash
## glnexus.conf
#!/bin/bash --login

#SBATCH --job-name=gln
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=6:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=980G
#SBATCH --export=NONE
module load bcftools/1.15--haf5b3da_0
srun --export=all -n 1 -c 64 glnexus_cli --config gatk --bed TraesFLD2D01G513900_updown1Mb.bed -m 980 --threads 64 2Dvcf/*.g.vcf.gz > more_TraesFLD2D01G513900_updown1Mb.bcf
srun --export=all -n 1 -c 64 bcftools view --threads 64 more_TraesFLD2D01G513900_updown1Mb.bcf | bgzip -@ 4 -c > more_TraesFLD2D01G513900_updown1Mb.vcf.gz
```
## snp annotation
### build snpeff database
```bash
cp cds.fa genes.gff protein.fa sequences.fa sequences.fa.fai /data/tools/snpEff/data/fielder
echo "fielder.genome : fielder" >> snpEff.config
java -jar snpEff.jar build -gff3 -v fielder
java -Xmx8g -jar snpEff.jar fielder /data/wheat/fielder/TraesFLD2D01G513900.vcf.gz >  /data/wheat/fielder/TraesFLD2D01G513900.annotated.vcf
bgzip TraesFLD2D01G513900.annotated.vcf
tabix -C TraesFLD2D01G513900.annotated.vcf.gz
bcftools annotate -x ^FORMAT/GT TraesFLD2D01G513900.annotated.vcf.gz > TraesFLD2D01G513900.annotated.simple.vcf
```
