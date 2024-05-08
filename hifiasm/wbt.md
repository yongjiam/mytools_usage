## genome assembly WBT wild barley
### 1.Data
```
-rw-r--r-- 1 yjia pawsey0399 26G Jan 14  2023 WBT_m64268e_230112_082652.hifi_reads.bam
-rw-r--r-- 1 yjia pawsey0399 26G Jan 18  2023 WBT_m64292e_230115_234637.hifi_reads.bam
-rw-r--r-- 1 yjia pawsey0399 25G Jan 19  2023 WBT_m64292e_230117_084848.hifi_reads.bam
-rw-r--r-- 1 yjia pawsey0399 26G Jan 20  2023 WBT_m64292e_230118_175426.hifi_reads.bam

-rw-r--r-- 1 yjia pawsey0399 123G May  3 11:01 Sample_8_WBT_S8_R1_001.trimmed.fq.gz
-rw-r--r-- 1 yjia pawsey0399 128G May  3 11:49 Sample_8_WBT_S8_R2_001.trimmed.fq.gz
```
### 2.Preprocessing
#### trim hic reads
>fastp -i ${R1} -I ${R2} -o ${R1_out} -O ${R2_out} --thread ${SLURM_NTASKS_PER_NODE} --trim_front1 5 --trim_front2 5 --cut_mean_quality 30 --qualified_quality_phred 20 --detect_adapter_for_pe --cut_front --cut_tail --length_required 30 --correction --cut_window_size 5 --overrepresentation_analysis --trim_poly_g -j ${Report_prefix}.fastp.json -h ${Report_prefix}.fastp.html
#### merge hifi fastq from bams
```
cat X201SC22091003-Z01-F004_02.tar.partaa X201SC22091003-Z01-F004_02.tar.partab > merged.tar
tar -xvf merged.tar
conda activate base
srun --export=all -n 1 -c 10 bam2fastq --split-barcodes -o out \
WBT_m64268e_230112_082652.hifi_reads.bam WBT_m64292e_230115_234637.hifi_reads.bam WBT_m64292e_230117_084848.hifi_reads.bam WBT_m64292e_230118_175426.hifi_reads.bam
```
### 3.hifiasm assembling
```
module load samtools/1.15--h3843a85_0
module load singularity/3.11.4-slurm
#srun --export=all -n 1 -c 128  samtools fastq -@ 128 ./ccs/m64257e_211030_130656.ccs.bam > hifi_ccs.fastq
srun --export=all -n 1 -c 64 singularity exec /scratch/pawsey0399/yjia/huyou/containers/hifiasm_latest.sif hifiasm -o wbt.asm -t 64 -l0 \
	--primary \
	--h1 /scratch/pawsey0399/yjia/WBT/Sample_8_WBT_S8_R1_001.trimmed.fq.gz \
	--h2 /scratch/pawsey0399/yjia/WBT/Sample_8_WBT_S8_R2_001.trimmed.fq.gz \
	/scratch/pawsey0399/yjia/WBT/hifi_reads/out.fastq.gz
```
### 4.purge duplication
```
## or install by even in setonix
git clone https://github.com/dfguan/purge_dups.git
cd purge_dups/src && make

## run purge_dups on draft asm
hap_asm=huyou_k19.asm.hic.hap2.p_ctg.fasta
HIFI=hifi_ccs.fastq
HISTPY=/scratch/pawsey0399/yjia/tools/purge_dups/scripts/hist_plot.py
##step1
srun --export=all -n 1 -c 128   minimap2 -xasm20 -t 128 $hap_asm $HIFI | gzip -c - > $HIFI.paf.gz
srun --export=all -n 1 -c 128 pbcstat $HIFI.paf.gz
srun --export=all -n 1 -c 128 calcuts PB.stat > cutoffs 2>calcults.log
####srun --export=all -n 1 -c 128 python3 $HISTPY -c cutoffs PB.stat PB.cov.png ## run this to determine the low, mid, high cutoffs
####srun --export=all -n 1 -c 128 calcuts -l 5 -m 25 -u 204 PB.stat > cutoffs_manual 2>calcults.log
srun --export=all -n 1 -c 128 split_fa $hap_asm > $hap_asm.split
srun --export=all -n 1 -c 128 minimap2 -xasm5 -DP $hap_asm.split $hap_asm.split | gzip -c - > $hap_asm.split.self.paf.gz
##step2
srun --export=all -n 1 -c 128 purge_dups -2 -T cutoffs -c PB.base.cov $hap_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
#####srun --export=all -n 1 -c 128 purge_dups -2 -T cutoffs_manual -c PB.base.cov $hap_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
##step3
srun --export=all -n 1 -c 128 get_seqs -e dups.bed $hap_asm
##step4 Merge hap.fa and $hap_asm and redo the above steps to get a decent haplotig set
```
### 5.hic scaffolding 
#### build singularity image
```
singularity build yjhicpipe.sif docker://yongjia111/yjhicpipe:latest
```
#### nextflow.sh
```
. /opt/conda/etc/profile.d/conda.sh
conda activate hic-scaffolding-nf
export NXF_HOME=/scratch/pawsey0399/yjia/WBT/hifionly/
nextflow run WarrenLab/hic-scaffolding-nf \
    -c /scratch/pawsey0399/yjia/WBT/hifionly/nextflow.config \
    --contigs /scratch/pawsey0399/yjia/WBT/hifionly/wbt_hifionly.asm.p_ctg.fasta \
    --r1Reads /scratch/pawsey0399/yjia/WBT/Sample_8_WBT_S8_R1_001.trimmed.fq.gz \
    --r2Reads /scratch/pawsey0399/yjia/WBT/Sample_8_WBT_S8_R2_001.trimmed.fq.gz \
    --juicer-tools-jar /scratch/pawsey0399/yjia/WBT/juicer_tools_1.22.01.jar \
    --extra-yahs-args "-e GATC"
```
#### nextflow.config
```
process {
    memory = '490 GB'
    time = '1d'

    withName: 'CHROMAP_ALIGN' {
        cpus = 64
        publishDir = [ path: 'out_hifionly/chromap', mode: 'copy' ]
    }
    withName: 'YAHS_SCAFFOLD' { publishDir = [ path: 'out_hifionly/scaffolds', mode: 'copy' ] }
    withName: 'JUICER_PRE' { publishDir = [ path: 'out_hifionly/juicebox_input', mode: 'copy' ] }
    withName: 'PRINT_VERSIONS' { publishDir = [ path: 'out_hifionly/', mode: 'copy' ] }
    withName: 'ASSEMBLY_STATS' { publishDir = [ path: 'out_hifionly/scaffolds', mode: 'copy' ] }
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
#### run hic pipeline
```
##hicpipe.conf
module load singularity/3.11.4-slurm
IMAGE=/scratch/pawsey0399/yjia/WBT/yjhicpipe.sif
srun --export=all -n 1 -c 64 singularity exec $IMAGE bash nextflow.sh
```
