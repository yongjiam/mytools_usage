## filter for black white lines
paste white_black.ids white_black.ids  > tmp && mv tmp white_black.ids
plink --bfile 2071_all_bed --maf 0.01 --geno 0.5 --keep white_black.ids --make-bed --out filtered_samples

awk '{print $1"_"$1, $2}' white_black.phenotype.txt > tmp && mv tmp white_black.phenotype.txt ## somehow the hapmap ids were doubled

## plink binary to vcf
plink --bfile filtered_samples --recode vcf --out filtered_samples
## vcf to hapmap
run_pipeline.pl -Xmx2G -importGuess filtered_samples.vcf -export filtered_samples -exportType HapmapDiploid
cut -f2-3 white_black.phenotype.txt > tmp && mv tmp white_black.phenotype.txt ## phenotype only need one ID column versus 2 in plink

## phenotype data
IID_IID white-black
YL0413_YL0413 1
YL0414_YL0414 1
YL0417_YL0417 1

## genotype data

## pangenome synteny
```
## data download
https://www.zenodo.org/record/7367881
https://www.nature.com/articles/s41588-023-01423-w#data-availability

## decompress tar.gz gffs and fasta files
ls --color=never *.tar.gz|while read R;do tar -xzvf $R;done

## extract protein files, sample_id file contain windows symbols and not recoganized by bash
gffread -w C10.trans.fa -x C10.cds.fa -y C10.pep.fa -g ./C10.fa C10.gff -F
gffread -w C12.trans.fa -x C12.cds.fa -y C12.pep.fa -g ./C12.fa C12.gff -F
gffread -w C13.trans.fa -x C13.cds.fa -y C13.pep.fa -g ./C13.fa C13.gff -F

## remove alternative transcript from pep.fa
cat sample_list |while read R;do (bioawk -c fastx '{if ($name ~ /T01/) print ">"$name"\n"$seq}' $R".pep.fa" > $R".pep.T01.fa");done
cat fasta_file_list |while read R;do (bioawk -c fastx '{gsub(/\./, "", $seq); print ">"$name"\n"$seq}' $R > tmp);mv tmp $R;done

## some sequence containing "." and remove them
bioawk -c fastx '{gsub(/\./, "", $seq); print ">"$name"\n"$seq}' pan110.pep_T01.fa > clean_pan110.pep_T01.fa
cat fasta_file_list |while read R;do (bioawk -c fastx '{gsub(/\./, "", $seq); print ">"$name"\n"$seq}' $R > tmp);mv tmp $R;done

## diamond blastp
srun --export=all -n 1 -c 64  diamond makedb --in pan110.pep_T01.fa --threads 64 -d pan110
srun --export=all -n 1 -c 64  diamond blastp -d pan110 -q clean_pan110.pep_T01.fa --threads 64 --evalue 1e-10 --outfmt 6 -o pan110.blast

## pangenes
miniprot -t64 -d C10.mpi C10.fa ## index
miniprot --outs=0.95 --no-cs -Iut64 L10.fa cd-hit/db95.fasta >L10.paf ## miniprot blast

conda activate base
srun --export=all -n 1 -c 64 bash miniprot_run.sh &> log_run.txt
srun --export=all -n 1 -c 64   pangene -p0 *.paf > graph95.gfa
srun --export=all -n 1 -c 64   k8 /software/projects/pawsey0399/yjia/tools/pangene/pangene.js gfa2matrix graph95.gfa > graph95_gene_presence_absence.Rtab
srun --export=all -n 1 -c 64   k8 /software/projects/pawsey0399/yjia/tools/pangene/pangene.js call graph95.gfa > graph95_bubble.txt
```
