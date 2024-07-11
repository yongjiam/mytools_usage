## hic nextflow_join.sh on nimbus
```
export NXF_HOME=/data/huyou/juicer2
nextflow run WarrenLab/hic-scaffolding-nf \
    --contigs /data/huyou/hifiasm_purged_hap1hap2merge/purged_join.fa \
    --r1Reads /data/huyou/raw_data/HIC/changshanhuyou-1/changshanhuyou-1_R1.fastq.gz \
    --r2Reads /data/huyou/raw_data/HIC/changshanhuyou-1/changshanhuyou-1_R2.fastq.gz \
    --juicer-tools-jar /data/huyou/juicer2/juicer_tools_1.22.01.jar \
    --extra-yahs-args "-e GATC"
```
## juicebox on mac
```
## run minimap2  without review
conda activate ragtag
ragtag.py scaffold $REF join_scaffolds_review.fa -t 30 -o ./ragtag_hifihic_join

## plot ragtag.scaffold.asm.paf from the ragtag output to guide juicebox
out_JBAT.review.assembly ## for juicer_post

```
## juicer_post.sh on nimbus
```
conda activate hic-scaffolding-nf
juicer post -o out_JBAT_review juicebox_input/out_JBAT.review.assembly juicebox_input/out_JBAT.liftover.agp /data/huyou/hifiasm_purged_hap1hap2merge/purged_join.fa

## produce out_JBAT_review.FINAL.fa / join_ragtag.sh
```
## merge unanchored scaffolds
```
## merge all unanchored scaffolds
sed -i 's/chrUn_RagTag/scaffold_RagTag/' join_ragtag.fasta
```
```
#!/bin/bash

# Variables
input_fasta="join_ragtag.fasta"    # Replace with your input FASTA file path
#input_fasta="tmp.fasta"
output_fasta="join.fasta"  # Replace with your desired output file path
#output_fasta="processed_tmp.fasta"
gap_size=100                           # Number of 'N's to insert between scaffold sequences

# Create a gap sequence of 100 'N's
gap=$(printf 'N%.0s' $(seq 1 $gap_size))

# Temporary files for scaffold and non-scaffold sequences
temp_scaffold_fasta="temp_scaffold.fasta"
temp_non_scaffold_fasta="temp_non_scaffold.fasta"

# Extract sequences with "scaffold" in their ID
bioawk -c fastx '$name ~ /scaffold/ {print ">"$name; print $seq}' $input_fasta > $temp_scaffold_fasta

# Extract sequences without "scaffold" in their ID
bioawk -c fastx '$name !~ /scaffold/ {print ">"$name; print $seq}' $input_fasta > $temp_non_scaffold_fasta

# Append non-scaffold sequences to the output file
cat $temp_non_scaffold_fasta > $output_fasta

# Merge scaffold sequences with 100 'N's between them
{
    echo ">chrUn_RagTag"
    bioawk -c fastx '{
        if (NR > 1) printf("%s", "'$gap'");
        printf("%s", $seq);
    }' $temp_scaffold_fasta
    echo
} >> $output_fasta

# Cleanup temporary files
rm $temp_scaffold_fasta $temp_non_scaffold_fasta

echo "Processed sequences saved to $output_fasta"
```
