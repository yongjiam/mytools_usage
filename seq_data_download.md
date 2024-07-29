## guide to download bulk sequencing data
### Option1: use NCBI sra-tools
```bash
## install sra-tools
conda install bioconda::sra-tools

## get a list of SRR ID
#go to NCBI https://www.ncbi.nlm.nih.gov/
#search project ID, such as PRJEB26543
#click the "Send to" button, select "File" -> "Accession list"

## download SRA
cat SRR_ids|while read R;do prefetch -O /data/igenome/emmer_exome $R;done > log.txt

## download fastq
cat SRR_ids|while read R;do fastq-dump --split-files $R && gzip $R"_1.fastq" && gzip $R"_2.fastq"

## extract fastq from sra file
cat SRR_ids|while read R;do echo "fastq-dump --gzip --skip-technical --readids --dumpbase --split-3 --clip  ."/$R/$R".sra";done > fastq-dump_commands
```
### Option2: use sra explorer 
```bash
#Go to at https://sra-explorer.info/
#Search Project ID, such as PRJNA1029507
#save the list of fastq link addresses
#download using wget/curl
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA005878/CRR290051/CRR290051_f1.fastq.gz
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com:443/data/gsapub/CRA005878/CRR290054/CRR290054_f1.fastq.gz
...
## or in batch in setonix
srun --export=all -n 1 -c 10  xargs -P 10 -I {} sh -c "{}" < wget_commands
## or in batch in linux
cat curls | xargs -n 1 -P 4 ./download.sh

####download.sh
#!/bin/bash
wget "$1"
```

### Option3: use ENA website to get fastq links
```bash
# go to https://www.ebi.ac.uk/ena/browser/home
# search for project id
# save fastq links
# then refer to option2
```


