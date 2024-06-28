## use gemoma to prediction genes in 76 barley genomes

### whole genome prediction based on morexV3
```bash
## gemoma_all.conf
#!/bin/bash --login

#SBATCH --job-name=VARIETY
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

#GEMOMAP=/scratch/pawsey0399/yjia/tools/gemoma-1.9-0/GeMoMa-1.9.jar
GEMOMAP="/scratch/pawsey0399/yjia/tools/gemoma18/GeMoMa-1.8.jar"
srun --export=all -n 1 -c 128   java -Xmx230G -jar $GEMOMAP CLI GeMoMaPipeline threads=128 tblastn=False \
	AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=VARIETY01G \
	p=true \
	pc=true \
	o=true \
	t=genome.fasta.gz \
	outdir=gemoma_output2_allMorexV3 \
	s=own i=morex a=/scratch/pawsey0399/yjia/shunlin/morexV3/all.gff3.gz g=/scratch/pawsey0399/yjia/shunlin/morexV3/genome.fasta.gz
```
```
## copy gemoma_all.conf to all genome folders
ls --color=never -d */|while read R;do cp gemoma_all.conf $R;done
## add variety name
for i in $(ls --color=never -d */ | grep -v GeMoMa);do VAR=$(basename $i); sed -i "s/VARIETY/$VAR/" $i"gemoma_all.conf";done
## submit jobs
for i in $(ls --color=never -d */ | grep -v GeMoMa);do cd $i; sbatch gemoma_all.conf;cd -;done
```
