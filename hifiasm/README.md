## tuning hifiasm parameters to get better assembly
https://hifiasm.readthedocs.io/en/latest/faq.html#faq

### which output assembly to use
If parental data is available, "*dip.hap*.p_ctg.gfa" produced in trio-binning mode should be always preferred. \

Otherwise if Hi-C data is available, "*hic.hap*.p_ctg.gfa" produced in Hi-C mode is the best choice. Both trio-binning mode and Hi-C mode generate fully-phased assemblies. \

If you only have HiFi reads, hifiasm in default outputs "*bp.hap*.p_ctg.gfa". The primary/alternate assemblies can be also produced by using --primary. All these HiFi-only assemblies are not fully-phased. See blog here for more details.

```bash
By default, hifiasm use -l3 for purge duplication, which use -s 0.55 similarity cuttoff
You can use lower -s for more aggressive purging or -l0 for no purging

# Assemble inbred/homozygous genomes (-l0 disables duplication purging)
hifiasm -o CHM13.asm -t32 -l0 CHM13-HiFi.fa.gz 2> CHM13.asm.log

# Assemble heterozygous genomes with built-in duplication purging (-s 0.55)
hifiasm -o HG002.asm -t32 HG002-file1.fq.gz HG002-file2.fq.gz
```
