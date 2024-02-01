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

