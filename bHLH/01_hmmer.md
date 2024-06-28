## 01. hmmer search for gene family in reference genome
```
## workdir: /Users/yongjia/Desktop/Programs/HMMER
## workdir: /Volumes/Elements5T/barley_annotation_V3/gene_annotation
hmmfetch Pfam-A.hmm PLATZ > ./hmm_files/PLATZ.hmm
hmmpress ./hmm_files/PLATZ.hmm
hmmscan --domtblout hmmscan_PLATZ_output.txt  -E 1e-5 /Users/yongjia/Desktop/Programs/HMMER/hmm_files/PLATZ.hmm  Hv_Morex.pgsb.Jul2020.aa.fa
```
## hmmscan output parser
```
## to target gene cds and protein sequence, and also domain sequences
https://github.com/yongjiam/jupyter_notebooks/blob/main/hmmscan_output_parser.ipynb
```
