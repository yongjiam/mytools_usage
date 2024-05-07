## tuning hifiasm parameters to get better assembly
```bash
By default, hifiasm use -l3 for purge duplication, which use -s 0.55 similarity cuttoff
You can use lower -s for more aggressive purging or -l0 for no purging

# Assemble inbred/homozygous genomes (-l0 disables duplication purging)
hifiasm -o CHM13.asm -t32 -l0 CHM13-HiFi.fa.gz 2> CHM13.asm.log

# Assemble heterozygous genomes with built-in duplication purging (-s 0.55)
hifiasm -o HG002.asm -t32 HG002-file1.fq.gz HG002-file2.fq.gz
```
