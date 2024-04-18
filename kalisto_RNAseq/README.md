# pipeline of RNAseq analysis using kalisto

## prepare sample_pair file
```
ls *.gz | paste - - > sample_pair
```

## trim fastq with fastp
```
while read R1 R2;
do
	fastp \
-i ./$R1 \
-I ./$R2 \
-o ./trimmed.$R1 \
-O ./trimmed.$R2 \
--qualified_quality_phred 20 \
--cut_front \
--cut_tail \
--length_required 30 \
--detect_adapter_for_pe \
--thread 32
done < sample_pairs
```
## quantification using kalisto
```
##build index
kallisto index -i transcripts.idx transcripts.fa.gz

## quantification
ls 'trimmed*fq.gz' | paste - - | while read R1 R2;
do
	D=$(echo $R1 | cut -d '.' -f2)
	kallisto \
	quant \
	-i transcripts.idx \
	-t 32 \
	-b 100 \
	-o $D \
	$R1 \
	$R2
done
```
## merge multiple sample into gene expression matrix
```
## merge whole abundance.tsv for different SRR samples
for i in $(ls --color=never -d */);do SRR=$(echo $i | cut -d "_" -f1); SAM=$(grep -m1 $SRR file_report.txt | awk '{print $NF}'; \
(sed "s/tpm/$SAM/" $i"abundance.tsv" > $i$SRR".tpm");done

## merge single gene abundance.tsv for different SRR samples
for i in $(ls --color=never -d */);do SRR=$(basename $i); SAM=$(grep -m1 $SRR ../filereport_read_run_PRJNA889532_tsv.txt | awk '{print
$NF}');echo $SAM $SRR $(grep -m1 HORVU.MOREX.r3.4HG0409010.1 $i"abundance.tsv");done > tong_PRJNA889532_tpm.tsv
```
