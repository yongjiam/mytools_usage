## tuning hifiasm parameters to get better assembly
https://hifiasm.readthedocs.io/en/latest/faq.html#faq

### which output assembly to use
```bash
srun --export=all -n 1 -c 64 singularity exec --bind ${PWD}:${PWD} hifiasm_latest.sif hifiasm -o huyou.asm -t 64 \
	--h1 ./HIC/changshanhuyou-1_R1.fq.gz \
	--h2 ./HIC/changshanhuyou-1_R2.fq.gz \
	hifi_ccs.fastq
```
If parental data is available, \*dip.hap\*.p_ctg.gfa produced in trio-binning mode should be always preferred. \

Otherwise if Hi-C data is available, \*hic.hap\*.p_ctg.gfa produced in Hi-C mode is the best choice. Both trio-binning mode and Hi-C mode generate fully-phased assemblies. \

If you only have HiFi reads, hifiasm in default outputs \*bp.hap\*.p_ctg.gfa. The primary/alternate assemblies can be also produced by using --primary. All these HiFi-only assemblies are not fully-phased. See blog here for more details.

### purging duplication level with -l and -s
```bash
By default, hifiasm use -l3 for purge duplication, which use -s 0.55 similarity cuttoff
You can use lower -s for more aggressive purging or -l0 for no purging

# Assemble inbred/homozygous genomes (-l0 disables duplication purging)
hifiasm -o CHM13.asm -t32 -l0 CHM13-HiFi.fa.gz 2> CHM13.asm.log

# Assemble heterozygous genomes with built-in duplication purging (-s 0.55)
hifiasm -o HG002.asm -t32 HG002-file1.fq.gz HG002-file2.fq.gz
```
### hifiasm misidentifies coverage threshold for homozygous reads
```bash
## check hifiasm log file
## use -k to manually set homozygous read coverage, heterozygous coverage 2x homozygous coverage
srun --export=all -n 1 -c 16 singularity exec --bind ${PWD}:${PWD} ./containers/hifiasm_latest.sif hifiasm -o huyou_k19.asm -t 16 \
	--h1 ./HIC/changshanhuyou-1_R1.fq.gz \
	--h2 ./HIC/changshanhuyou-1_R2.fq.gz \
	-k 19 \
	hifi_ccs.fastq
```
