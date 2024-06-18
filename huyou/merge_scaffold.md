## merge small scaffolds
```
#!/bin/bash

# Variables
input_fasta="ragtag_hifihic_join/ragtag.scaffold.fasta"    # Replace with your input FASTA file path
#input_fasta="tmp.fasta"
output_fasta="processed_ragtag_hifihic_join.fasta"  # Replace with your desired output file path
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
    echo ">ragtag_chrUn"
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
