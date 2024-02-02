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
REF="Triticum_aestivum.IWGSC.dna_rm.toplevel.fa"
srun --export=all -n 1 -c 30  gatk MarkDuplicates -I $BAM -O "dedup_"$BAM -M marked_dup_metrics.txt
srun --export=all -n 1 -c 30 samtools addreplacerg -O BAM -@ 30 -o "RG_dedup_"$BAM -r '@RG\tID:HT621\tSM:HT621\tPL:ILLUMINA' "dedup_"$BAM
srun --export=all -n 1 -c 30 samtools index -c -@ 30 "RG_dedup_"$BAM
srun --export=all -n 1 -c 30 gatk SplitNCigarReads -R $REF -I "RG_dedup_"$BAM -O "split_RG_dedup_"$BAM ## note: the dictionary file should be Triticum_aestivum.IWGSC.dna_rm.toplevel.dict
srun --export=all -n 1 -c 30 samtools index -c -@ 30 "split_RG_dedup_"$BAM
srun --export=all -n 1 -c 30 gatk HaplotypeCaller -R $REF -I "split_RG_dedup_"$BAM -ERC GVCF -O "split_RG_dedup_"$BAM".g.vcf"
```
