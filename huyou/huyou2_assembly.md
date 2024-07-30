## huyou2 assembly, data source
https://www.nature.com/articles/s41597-024-03437-3
https://www.ebi.ac.uk/ena/browser/view/PRJNA1091318

## data
SRR28430796	Huyou_pacbio_HiFi
SRR28430795	Huyou_illumina
SRR28430802	Huyou_ripe_fruit_RNA
SRR28430798	Huyou_leaf_RNA
SRR28430799	Huyou_root_RNA
SRR28430801	Huyou_flower_RNA
SRR28430797	Huyou_HiC
SRR28430800	Huyou_stem_RNA
SRR28430803	Huyou_unripe_fruit_RNA

## hifiasm assembly
```
#!/bin/bash --login

#SBATCH --job-name=huyou
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=14:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=980G
#SBATCH --export=NONE

module load samtools/1.15--h3843a85_0
module load singularity/4.1.0-slurm
#srun --export=all -n 1 -c 128  samtools fastq -@ 128 ./ccs/m64257e_211030_130656.ccs.bam > hifi_ccs.fastq
IMAGE=/scratch/pawsey0399/yjia/huyou/containers/hifiasm_latest.sif
srun --export=all -n 1 -c 32 singularity exec --bind ${PWD}:${PWD} $IMAGE hifiasm -o huyou2.asm -t 32 \
	--h1 ./SRR28430797_1.fastq.gz \
	--h2 ./SRR28430797_2.fastq.gz \
	SRR28430796_subreads.fastq.gz
```
## purge duplicates
```
## gfa to fasta
awk 'BEGIN {FS="\t"} $1=="S" {print ">"$2"\n"$3}' huyou2.asm.hic.hap1.p_ctg.gfa > huyou2.asm.hic.hap1.p_ctg.fasta
awk 'BEGIN {FS="\t"} $1=="S" {print ">"$2"\n"$3}' huyou2.asm.hic.hap2.p_ctg.gfa > huyou2.asm.hic.hap2.p_ctg.fasta

## purge dups
#!/bin/bash --login

#SBATCH --job-name=purge
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

conda activate assembly

## run purge_dups on draft asm
hap_asm=huyou2.asm.hic.hap1.p_ctg.fasta
HIFIseq=SRR28430796_subreads.fastq.gz
NAME=huyou2
HISTPY=/scratch/pawsey0399/yjia/tools/purge_dups/scripts/hist_plot.py

##step1
srun --export=all -n 1 -c 128   minimap2 -xasm20 -t 128 $hap_asm $HIFIseq | gzip -c - > $NAME.paf.gz
srun --export=all -n 1 -c 128 pbcstat $NAME.paf.gz
srun --export=all -n 1 -c 128 calcuts PB.stat > cutoffs 2>calcults.log
srun --export=all -n 1 -c 128 python3 $HISTPY -c cutoffs PB.stat PB.cov.png ## run this to determine the low, mid, high cutoffs

####srun --export=all -n 1 -c 128 calcuts -l 5 -m 25 -u 204 PB.stat > cutoffs_manual 2>calcults.log
srun --export=all -n 1 -c 128 split_fa $hap_asm > $hap_asm.split
srun --export=all -n 1 -c 128 minimap2 -xasm5 -DP $hap_asm.split $hap_asm.split | gzip -c - > $hap_asm.split.self.paf.gz

##step2
srun --export=all -n 1 -c 128 purge_dups -2 -T cutoffs -c PB.base.cov $hap_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
#####srun --export=all -n 1 -c 128 purge_dups -2 -T cutoffs_manual -c PB.base.cov $hap_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

##step3
srun --export=all -n 1 -c 128 get_seqs -e dups.bed $hap_asm
```
