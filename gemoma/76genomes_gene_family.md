## use gemoma to predict target genes only
```
#!/bin/bash
## gemoma.sh
GEMOMAP="/scratch/pawsey0399/yjia/tools/gemoma18/GeMoMa-1.8.jar"
for i in $(ls --color=never -d */ | grep -v GeMoMa)
do
        java -jar $GEMOMAP CLI GeMoMaPipeline threads=128 AnnotationFinalizer.r=NO tblastn=False \
        t=$i"genome.fasta.gz" outdir=$i"gemoma_output_PLATZ"  GAF.f="iAA>=0.75 and pAA>=0.75" Extractor.f=False selected=PLATZ_morexV3_gene_id.txt \
        a=../morexV3/all.gff3.gz g=../morexV3/genome.fasta.gz
done
```
```
#### gemoma.conf
#!/bin/bash --login

#SBATCH --job-name=gemoma
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

#conda activate base
srun --export=all -n 1 -c 128  bash gemoma.sh &> log.txt
```
```
## sum_gemoma_output.sh
echo "gemoma_output_folder_name:"
read OUTPUT
echo "query gene list file:"
read GENES
## add variety ID to protein fasta and merge
for i in $(ls --color=never -d */ | grep -v GeMoMa);do VAR=$(basename $i); sed "s/>/>$VAR\_/" $i$OUTPUT"/predicted_proteins.fasta";done > "merged_"$OUTPUT".fasta"

## add variety ID to gff and merge
for i in $(ls --color=never -d */ | grep -v GeMoMa);do VAR=$(basename $i)"_"; (awk '!/^#/{print $0}' $i$OUTPUT"/final_annotation.gff" | sed "s/^/$VAR/");done > "merged_"$OUTPUT".gff"

# another way
for i in $(ls --color=never -d */ | grep -v GeMoMa)
do
	VAR=$(basename $i)
	TYPE=$(grep $VAR variety_name_and_type_match.csv|cut -d ',' -f4)
	cat $GENES | while read GENE;do COPY=$(grep $VAR "merged_"$OUTPUT".fasta"| grep $GENE | wc -l); echo -e $GENE '\t' $COPY '\t' $VAR '\t' $TYPE;done
done > "merged_"$OUTPUT".count"
```
```
## merge_output_fasta.sh
echo "gemoma_output_folder_name:"
read OUTPUT
echo "query_gene_list_file:"
read GENES
## add variety ID to protein fasta and merge
for i in $(ls --color=never -d */);do VAR=$(basename $i); sed "s/>/>$VAR\_/" $i$OUTPUT"/predicted_proteins.fasta";done > "merged_"$OUTPUT".fasta"

## read gene id and count copy number and add type info
cat $GENES | while read R1;
do
	# use uniq -c to count VAR number
	grep $R1 "merged_"$OUTPUT".fasta" | cut -d "_" -f1 | sed 's/>//' | sort | uniq -c | while read C V; do echo -e $R1 '\t' $C '\t' $V '\t' $(grep $V variety_name_and_type_match.csv|cut -d ',' -f4);done > $R1".count"
	# sort by type and var name
	sort -k4,4 -k2,2 -o $R1".count" $R1".count"
done
```
### plot the output of gene family search
```
## CNV_distribution_AMY6H.ipynb
```
