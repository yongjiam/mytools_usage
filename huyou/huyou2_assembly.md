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
## hic pipeline
```
## nextflow_hap2.sh
export NXF_HOME=$PWD
nextflow run WarrenLab/hic-scaffolding-nf \
    --contigs ./hap2_purged.fa \
    --r1Reads ./rawdata/SRR28430797_1.fastq.gz \
    --r2Reads ./rawdata/SRR28430797_2.fastq.gz \
    --juicer-tools-jar /data/huyou/juicer2/juicer_tools_1.22.01.jar \
    --extra-yahs-args "-e GATC"

## nextflow.config
process {
    memory = '100 GB'
    time = '1d'

    withName: 'CHROMAP_ALIGN' {
        cpus = 32
        publishDir = [ path: 'out_hap2/chromap', mode: 'copy' ]
    }
    withName: 'YAHS_SCAFFOLD' { publishDir = [ path: 'out_hap2/scaffolds', mode: 'copy' ] }
    withName: 'JUICER_PRE' { publishDir = [ path: 'out_hap2/juicebox_input', mode: 'copy' ] }
    withName: 'PRINT_VERSIONS' { publishDir = [ path: 'out_hap2/', mode: 'copy' ] }
    withName: 'ASSEMBLY_STATS' { publishDir = [ path: 'out_hap2/scaffolds', mode: 'copy' ] }
}

profiles {
    lewis {
        process {
            executor = 'slurm'
            queue = 'BioCompute'
            clusterOptions = '--account=warrenlab'
            conda = '/storage/hpc/group/warrenlab/users/esrbhb/mambaforge/envs/chromap-yahs'
        }

        conda.enabled = true

        params {
            juicerToolsJar = '/storage/htc/warrenlab/users/esrbhb/software/juicer_tools_1.11.09_jcuda.0.8.jar'
        }
    }

    conda {
        process.conda = "$baseDir/conda.yml"
        conda.enabled = true
    }
}

manifest {
    defaultBranch = 'main'
    homePage = 'https://github.com/WarrenLab/hic-scaffolding-nf'
    author = 'Edward S. Rice'
    version = '0.0.1'
}
```
## assign to chromosomes using ragtag against SWO
```
#!/bin/bash --login

#SBATCH --job-name=rag2
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=6:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

conda activate base
REF=/scratch/pawsey0399/yjia/huyou/pangenomes/huazhong_downloads/SWO.v3.0.genome.fa

srun --export=all -n 1 -c 64   ragtag.py scaffold $REF hap2_out_scaffolds_final.fa -t 64 -o ./ragtag_hap2 &> log2.txt
```
## gene model prediction using gemoma
```
#!/bin/bash --login

#SBATCH --job-name=hap1
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

module load star/2.7.10a--h9ee0642_0
module load samtools/1.15--h3843a85_0

srun --export=all -n 1 -c 64  STAR --runThreadN 64 --runMode genomeGenerate \
       --genomeDir ./hap1_updated_star_index  \
       --genomeSAindexNbases 13 \
       --genomeFastaFiles hap1_ragtag.fasta

## maping
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430798_1.fastq hifiasm_purge/trimmed.SRR28430798_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430798 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430799_1.fastq hifiasm_purge/trimmed.SRR28430799_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430799 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430800_1.fastq hifiasm_purge/trimmed.SRR28430800_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430800 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430801_1.fastq hifiasm_purge/trimmed.SRR28430801_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430801 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430802_1.fastq hifiasm_purge/trimmed.SRR28430802_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430802 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
srun --export=all -n 1 -c 64 STAR --runThreadN 64 --genomeDir ./hap1_updated_star_index \
        --readFilesIn hifiasm_purge/trimmed.SRR28430803_1.fastq hifiasm_purge/trimmed.SRR28430803_2.fastq \
        --outFileNamePrefix hap1updated_RNAmap_SRR28430803 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard

## merge bams
srun --export=all -n 1 -c 64 samtools merge -@ 64 hap1merged.bam hap1updated_RNAmap*Aligned.sortedByCoord.out.bam


GEMOMAP=/scratch/pawsey0399/yjia/tools/gemoma18/GeMoMa-1.8.jar
srun --export=all -n 1 -c 64 java -Xmx110G -jar $GEMOMAP CLI GeMoMaPipeline threads=64 tblastn=False \
	AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=H1Y \
	p=true \
	pc=true \
	o=true \
	t=hap1_ragtag.fasta \
	outdir=gemoma_hap1/ \
	s=own i=SWO a=SWO.v3.0.gene.model.gff3 g=SWO.v3.0.genome.fa \
	s=own i=ZGYCC a=ZGYCC.v2.0.gene.model.gff3 g=ZGYCC.v2.0.genome.fa \
	s=own i=JZ a=JZ.v1.0.gene.model.gff3 g=JZ.v1.0.genome.fa \
	r=MAPPED ERE.m=hap1merged.bam
```
## swap chromosomes between SD1 and SD2, based on kpart results
```
## change chromosome ID
#/scratch/pawsey0399/yjia/huyou2/swap_chroms/
sed 's/>chr/>SD1chr/;s/>scaff/>SD1scaff/' SD1_ragtag.fasta > SD1_changeID.fasta
sed 's/>chr/>SD2chr/;s/>scaff/>SD2scaff/' SD2_ragtag.fasta > SD2_changeID.fasta

## chr2, chr3, chr8
sed -i 's/SD1chr2/SD2chr2/; s/SD1chr3/SD2chr3/; s/SD1chr8/SD2chr8/' SD1_changeID.fasta
sed -i 's/SD2chr2/SD1chr2/; s/SD2chr3/SD1chr3/; s/SD2chr8/SD1chr8/' SD2_changeID.fasta


```
