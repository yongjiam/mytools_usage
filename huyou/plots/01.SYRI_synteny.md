## install
```
mamba create --name SYRI syri
mamba install plotsr
```
## alignment
```
minimap2 -ax asm5 --eqx -t 30 ragtag_hifihic_hap1.fasta ragtag_hifihic_hap2.fasta > ragtag_hifihic_hap1_hap2.sam
samtools view -b -@ 30 ragtag_hifihic_hap1_hap2.sam  > ragtag_hifihic_hap1_hap2.bam
```
## plot
```
syri -c ragtag_hifihic_hap1_hap2.bam -r ragtag_hifihic_hap1.fasta -q ragtag_hifihic_hap2.fasta -k -F B
plotsr --sr syri.out --genomes genomes.txt -o syri.pdf
```
```
## genomes.txt
#file	name	tags
ragtag_hifihic_hap1.fasta	hap1	lw:1.5
ragtag_hifihic_hap2.fasta	hap2	lw:1.5
```
