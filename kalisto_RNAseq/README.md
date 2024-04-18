# pipeline of RNAseq analysis using kalisto

## prepare sample_pair file
```
ls *.gz | paste - - > sample_pair
```
## prepare transcript.fa
```
#!/bin/bash
## change sequence ID in fasta
while read line; do
    if [[ $line == ">"* ]]; then
	    ID=">"$(echo "$line" | grep -o 'NLLv2_murdoch[^;]*')
	    echo $ID
    else
        echo "$line"
    fi
done < transcripts.fa
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
### notes on multi-mapped reads
https://www.sciencedirect.com/science/article/pii/S2001037020303032 \
https://www.nature.com/articles/nbt.3519 \
User
i am using kalisto tool to process RNA sequencing data. How does this tool handle reads that mapped to multiple loci or transcripts? And how does this compare to other RNA sequencing tools such star and hisat and tophat2
ChatGPT
Kalisto is a quantification tool specifically designed for RNA-Seq data analysis. It uses an EM algorithm to probabilistically assign reads to transcripts and estimate abundances. When dealing with reads that map to multiple loci or transcripts, Kalisto employs an expectation-maximization (EM) algorithm to probabilistically assign those reads to their most likely originating transcripts. This probabilistic approach accounts for multi-mapping reads by distributing them across potential transcript targets based on their likelihood of originating from each transcript.

Comparing Kalisto to other popular RNA-Seq alignment and quantification tools like STAR, HISAT, and Tophat2:

STAR: STAR is an aligner that performs both alignment and quantification. It employs a two-pass mapping strategy to improve alignment accuracy, allowing for the detection of novel splice junctions. STAR also provides a variety of output formats, including read counts and transcripts per million (TPM). However, unlike Kalisto, STAR does not directly estimate transcript abundances but rather provides aligned reads that can be used for downstream quantification.
HISAT: HISAT (Hierarchical Indexing for Spliced Alignment of Transcripts) is another aligner commonly used for RNA-Seq data analysis. Similar to STAR, HISAT utilizes a two-step mapping process to improve alignment accuracy, particularly for reads spanning splice junctions. It generates alignment files that can be used for downstream analysis, including quantification. Like STAR, HISAT does not directly estimate transcript abundances.
Tophat2: Tophat2 is an older aligner and quantifier for RNA-Seq data that aligns reads to the genome using the Bowtie2 aligner and then assembles transcripts and estimates abundances. However, it has been largely superseded by HISAT, which offers improved performance and accuracy.
Compared to these aligners, Kalisto provides a distinct advantage in its ability to directly estimate transcript abundances from RNA-Seq data through its probabilistic model. This can result in more accurate quantification, particularly in cases where reads map to multiple loci or transcripts, as Kalisto accounts for this ambiguity in its estimation process. Additionally, Kalisto is generally faster than aligners like STAR and HISAT, as it does not perform alignment but focuses solely on quantification. However, the choice of tool depends on the specific requirements of the analysis and the preferences of the user.




User
in the expectation-maximization (EM) algorithm used by kalisto, if a read could map to two transcripts with identical sequence, does it assign this read to both of the two transcripts ?
ChatGPT
In the expectation-maximization (EM) algorithm used by Kalisto, if a read could map to two transcripts with identical sequences, it does not assign the read to both transcripts equally. Instead, Kalisto assigns the read probabilistically to the transcripts based on their relative abundances in the sample.

Here's a simplified explanation of how this process works:

Initial Assignment: In the expectation step (E-step) of the EM algorithm, Kalisto initially assigns the read to both transcripts with equal probability.
Maximization: In the maximization step (M-step), Kalisto updates the abundance estimates of the transcripts based on the initial assignments of reads. This step involves recalculating the relative abundances of transcripts in the sample based on the assigned reads.
Updated Assignment: In the next iteration of the E-step, Kalisto uses the updated abundance estimates to probabilistically assign reads to transcripts. This assignment takes into account the updated abundance estimates, ensuring that reads are more likely to be assigned to transcripts with higher estimated abundances.
Convergence: The EM algorithm iterates between the E-step and M-step until convergence is reached, meaning that the assignment of reads and the abundance estimates of transcripts stabilize.
During this iterative process, even if two transcripts have identical sequences, Kalisto will assign reads to them based on their estimated abundances. Transcripts with higher estimated abundances will have a higher probability of receiving reads, reflecting their presumed higher expression levels in the sample. Therefore, the assignment of reads to transcripts is not necessarily equal in cases where transcripts have identical sequences but differ in abundance estimates.
