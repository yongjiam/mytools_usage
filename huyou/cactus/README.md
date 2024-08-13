## run minigraph cactus for selected citrus genomes
https://github.com/ComparativeGenomicsToolkit/cactus/tree/master
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md
minigraph cactus for aligning the same species
Progressive Cactus for aligning different species

## install cactus 
https://github.com/ComparativeGenomicsToolkit/cactus/releases
```
singularity build cactus29.sif docker:quay.io/comparative-genomics-toolkit/cactus:v2.9.0
```

## test run cactus
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/README.md
### location: 
(base) yongjia@ChengdaoServer:/media/yongjia/Elements/yongjia/cactus_test$
### download 4-t2t-orangs-mc-2023v2.seqfile
```
mPonAbe1_pri	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/655/GCA_028885655.2/GCA_028885655.2.fa.gz
mPonAbe1_alt	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/685/GCA_028885685.2/GCA_028885685.2.fa.gz
mPonPyg2.1	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/625/GCA_028885625.2/GCA_028885625.2.fa.gz
mPonPyg2.2	https://hgdownload.soe.ucsc.edu/hubs/GCA/028/885/525/GCA_028885525.2/GCA_028885525.2.fa.gz
```
### command cactus.sh
```
cactus-pangenome ./js-pg ./4-t2t-orangs-mc-2023v2.seqfile --outDir 4-t2t-orangs-mc-2023v2 --outName 4-t2t-orangs-mc-2023v2 --reference mPonAbe1_pri mPonAbe1_alt \
--noSplit --gbz clip full --gfa clip full --xg clip full --odgi --vcf --giraffe clip --haplo clip --vcfReference mPonAbe1_pri mPonAbe1_alt --logFile 4-t2t-orangs-mc-2023v2.log  \
--coordinationDir /data/tmp --batchLogsDir ./batch-logs --consMemory 230Gi --indexMemory 230Gi --mgMemory 230Gi --mgCores 60 --mapCores 8 --consCores 60 --indexCores 60 --giraffe clip
```
### run using singularity in tmux
```
singularity shell -B ${PWD}:/data /media/yongjia/Elements/yongjia/containers/cactus29.sif
bash cactus.sh &> log.txt
```
