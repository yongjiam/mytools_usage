## TEtools docker images
https://github.com/Dfam-consortium/TETools

## denovoRE.conf
```bash 
#!/bin/bash --login

#SBATCH --job-name=denovoRE
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

module load singularity/3.11.4-slurm
MAKER=/scratch/pawsey0399/yjia/huyou/containers/yongmaker.sif
GENOME=/scratch/pawsey0399/yjia/huyou/pangenomes/company/huyou.hap1.genome.fa
#srun --export=all -n 1 -c 64 singularity exec $MAKER bash denovoRE.sh

srun --export=all -n 1 -c 64 singularity build tetools.sif docker://dfam/tetools:latest ## build sif image
srun --export=all -n 1 -c 128 singularity run docker://dfam/tetools:latest bash denovoRE.sh ## run repeatmodeler
```

## ## denovoRE.sh
```bash
#. /opt/miniconda3/etc/profile.d/conda.sh
#conda activate maker
BuildDatabase -name HAP1 -engine ncbi /scratch/pawsey0399/yjia/huyou/pangenomes/company/huyou.hap1.genome.fa
RepeatModeler -threads 64 -engine ncbi -database HAP1 2>&1 | tee repeatmodeler.log
```
