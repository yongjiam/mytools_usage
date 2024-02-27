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
