## select target genome with chromosome scale data
(base) yongjia@ChengdaoServer:/media/yongjia/Elements/yongjia/cactus_test/citrus_genomes
```
ls *.fa|while read R;do assembly-stats $R;done > assembly-stats

grep -B2 N50 assembly-stats ## look for N50 and n values
```
```
## compared to previous orthofinder species,
## replaced Citrus_reticulata with newly published chromosomes
## https://academic.oup.com/hr/article/11/8/uhae177/7700677
## added Citrus_hindsii
## replaced Citrus_australasica with GCA_029618585.genome.fa


Citrus_reticulata	./Citrus_reticulata_Chachi.chr.fa
Citrus_australasica	./GCA_029618585.genome.fa
Citrus_clementina	./GCF.v1.0.genome.fa
Citrus_hindsii	./GJ.v2.0.genome.fa
H1Y	./hap1_updated.fasta
H2Y	./hap2_updated.fasta
Atalantia_buxfoliata	./HKC.v2.0.genome.fa
Citrus_grandis	./HWB.v1.0.genome.fa
SD1	./SD1_updated.fasta
SD2	./SD2_updated.fasta
Citrus_sinensis	./SWO.v3.0.genome.fa
Citrus_maxima	./XGF.v1.0.genome.fa
Poncirus_trifoliata	./ZK.v1.0.genome.fa
```
## ragtag to sort and assign genomes using SWO as reference
```
## install ragtag
mamba create -n ragtag bioconda::ragtag

### ragtag_list
Citrus_reticulata	./Citrus_reticulata_Chachi.chr.fa
Citrus_australasica	./GCA_029618585.genome.fa
Citrus_clementina	./GCF.v1.0.genome.fa
Citrus_hindsii	./GJ.v2.0.genome.fa
H1Y	./hap1_updated.fasta
H2Y	./hap2_updated.fasta
Atalantia_buxfoliata	./HKC.v2.0.genome.fa
Citrus_grandis	./HWB.v1.0.genome.fa
SD1	./SD1_updated.fasta
SD2	./SD2_updated.fasta
Citrus_maxima	./XGF.v1.0.genome.fa
Poncirus_trifoliata	./ZK.v1.0.genome.fa


## run ragtag
cat ragtag_list |while read R1 R2;do ragtag.py scaffold SWO.v3.0.genome.fa $R2 -t 60 -o $R1;done

## clean genomes
mkdir cleaned_citrus_genomes
cut -f1 ragtag_list|while read R;do sed 's/_RagTag//' $R"/ragtag.scaffold.fasta" > cleaned_citrus_genomes/$R".fa";done ## replace ragtag in IDs
mkdir cleaned_citrus_chromosome_only
ls *fa|while read R;do bioawk -c fastx '{if ($name ~ /chr[1-9]/) print ">"$name"\n"$seq}' $R > ../cleaned_citrus_chromosome_only/$R;done ## retain only chromosomes
```
## run cactus
```
## create seqfile for cactus
ls *.fa|while read R;do echo -e $(echo $R|cut -d '.' -f1)"\t./"$R;done > citrus_genomes.seqfile

## cactus.sh
cactus-pangenome ./js-pg ./citrus_genomes.seqfile --outDir citrus_pan --outName citrus_pan --reference Citrus_siniensis --noSplit --gbz clip full --gfa clip full --xg clip full --odgi --vcf --giraffe clip --haplo clip --vcfReference Citrus_siniensis --logFile citrus_pan.log  --coordinationDir /data/tmp --batchLogsDir ./batch-logs --consMemory 230Gi --indexMemory 230Gi --mgMemory 230Gi --mgCores 60 --mapCores 8 --consCores 60 --indexCores 60 --giraffe clip

##
mkdir tmp
singularity shell -B ${PWD}:/data /media/yongjia/Elements/yongjia/containers/cactus281.sif
Singularity> cd /data && bash cactus.sh
```
