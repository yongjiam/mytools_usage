## tuning hifiasm parameters to get better assembly
https://hifiasm.readthedocs.io/en/latest/faq.html#faq

### which output assembly to use
```bash
srun --export=all -n 1 -c 64 singularity exec --bind ${PWD}:${PWD} hifiasm_latest.sif hifiasm -o huyou.asm -t 64 \
	--h1 ./HIC/changshanhuyou-1_R1.fq.gz \
	--h2 ./HIC/changshanhuyou-1_R2.fq.gz \
	hifi_ccs.fastq
```
>If parental data is available, \*dip.hap\*.p_ctg.gfa produced in trio-binning mode should be always preferred.

>Otherwise if Hi-C data is available, \*hic.hap\*.p_ctg.gfa produced in Hi-C mode is the best choice. Both trio-binning mode and Hi-C mode generate fully-phased assemblies.

>If you only have HiFi reads, hifiasm in default outputs \*bp.hap\*.p_ctg.gfa. The primary/alternate assemblies can be also produced by using --primary. All these HiFi-only assemblies are not fully-phased. See blog here for more details.

### purging duplication level with -l and -s
```bash
By default, hifiasm use -l3 for purge duplication, which use -s 0.55 similarity cuttoff
You can use lower -s for more aggressive purging or -l0 for no purging

# Assemble inbred/homozygous genomes (-l0 disables duplication purging)
hifiasm -o CHM13.asm -t32 -l0 CHM13-HiFi.fa.gz 2> CHM13.asm.log

# Assemble heterozygous genomes with built-in duplication purging (-s 0.55)
hifiasm -o HG002.asm -t32 HG002-file1.fq.gz HG002-file2.fq.gz
```
### hifiasm misidentifies coverage threshold for homozygous reads, use -k or --hom-cov to fix this
```bash
## check hifiasm log file, hifiasm use -k 51 by default
## use -k to specify the k-mer length which will correct homozygous read coverage, heterozygous coverage 2x homozygous coverage
srun --export=all -n 1 -c 16 singularity exec --bind ${PWD}:${PWD} ./containers/hifiasm_latest.sif hifiasm -o huyou_k19.asm -t 16 \
	--h1 ./HIC/changshanhuyou-1_R1.fq.gz \
	--h2 ./HIC/changshanhuyou-1_R2.fq.gz \
	-k 19 \
	hifi_ccs.fastq
```
### hifiasm log file interpretation
https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#loginter

Hifiasm prints several information for quick debugging, including:

>k-mer plot: showing how many k-mers appear a certain number of times. For homozygous samples, there should be one peak around read coverage. For heterozygous samples, there should two peaks, where the smaller peak is around the heterozygous read coverage and the larger peak is around the homozygous read coverage. For example, issue10 indicates the heterozygous read coverage and the homozygous read coverage are 28 and 57, respectively. Issue49 is another good example. Weird k-mer plot like issue93 is often caused by insufficient coverage or presence of contaminants.

>homozygous coverage: coverage threshold for homozygous reads. Hifiasm prints it as: [M::purge_dups] homozygous read coverage threshold: X. If it is not around homozygous coverage, the final assembly might be either too large or too small. To fix this issue, please set --hom-cov to homozygous coverage.

>number of het/hom bases: how many bases in unitig graph are heterozygous and homozygous during Hi-C phased assembly. Hifiasm prints it as: [M::stat] # heterozygous bases: X; # homozygous bases: Y. Given a heterozygous sample, if there are much more homozygous bases than heterozygous bases, hifiasm fails to identify correct coverage threshold for homozygous reads. In this case, please set --hom-cov to homozygous coverage.

### hifiasm -h
```
Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>
Options:
  Input/Output:
    -o STR       prefix of output files [hifiasm.asm]
    -t INT       number of threads [1]
    -h           show help information
    --version    show version number
  Overlap/Error correction:
    -k INT       k-mer length (must be <64) [51]
    -w INT       minimizer window size [51]
    -f INT       number of bits for bloom filter; 0 to disable [37]
    -D FLOAT     drop k-mers occurring >FLOAT*coverage times [5.0]
    -N INT       consider up to max(-D*coverage,-N) overlaps for each oriented read [100]
    -r INT       round of correction [3]
    -z INT       length of adapters that should be removed [0]
    --max-kocc   INT
                 employ k-mers occurring <INT times to rescue repetitive overlaps [2000]
    --hg-size    INT(k, m or g)
                 estimated haploid genome size used for inferring read coverage [auto]
  Assembly:
    -a INT       round of assembly cleaning [4]
    -m INT       pop bubbles of <INT in size in contig graphs [10000000]
    -p INT       pop bubbles of <INT in size in unitig graphs [0]
    -n INT       remove tip unitigs composed of <=INT reads [3]
    -x FLOAT     max overlap drop ratio [0.8]
    -y FLOAT     min overlap drop ratio [0.2]
    -i           ignore saved read correction and overlaps
    -u           post-join step for contigs which may improve N50; 0 to disable; 1 to enable
                 [1] and [1] in default for the UL+HiFi assembly and the HiFi assembly, respectively
    --hom-cov    INT
                 homozygous read coverage [auto]
    --lowQ       INT
                 output contig regions with >=INT% inconsistency in BED format; 0 to disable [70]
    --b-cov      INT
                 break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [0]
    --h-cov      INT
                 break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [-1]
    --m-rate     FLOAT
                 break contigs at positions with <=FLOAT*coverage exact overlaps;
                 only work with '--b-cov' or '--h-cov'[0.75]
    --primary    output a primary assembly and an alternate assembly
    --ctg-n      INT
                 remove tip contigs composed of <=INT reads [3]
  Trio-partition:
    -1 FILE      hap1/paternal k-mer dump generated by "yak count" []
    -2 FILE      hap2/maternal k-mer dump generated by "yak count" []
    -3 FILE      list of hap1/paternal read names []
    -4 FILE      list of hap2/maternal read names []
    -c INT       lower bound of the binned k-mer's frequency [2]
    -d INT       upper bound of the binned k-mer's frequency [5]
    --t-occ      INT
                 forcedly remove unitigs with >INT unexpected haplotype-specific reads;
                 ignore graph topology; [60]
    --trio-dual  utilize homology information to correct trio phasing errors
  Purge-dups:
    -l INT       purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    -s FLOAT     similarity threshold for duplicate haplotigs in read-level [0.75 for -l1/-l2, 0.55 for -l3]
    -O INT       min number of overlapped reads for duplicate haplotigs [1]
    --purge-max  INT
                 coverage upper bound of Purge-dups [auto]
    --n-hap      INT
                 number of haplotypes [2]
  Hi-C-partition:
    --h1 FILEs   file names of Hi-C R1  [r1_1.fq,r1_2.fq,...]
    --h2 FILEs   file names of Hi-C R2  [r2_1.fq,r2_2.fq,...]
    --seed INT   RNG seed [11]
    --s-base     FLOAT
                 similarity threshold for homology detection in base-level;
                 -1 to disable [0.5]; -s for read-level (see <Purge-dups>)
    --n-weight   INT
                 rounds of reweighting Hi-C links [3]
    --n-perturb  INT
                 rounds of perturbation [10000]
    --f-perturb  FLOAT
                 fraction to flip for perturbation [0.1]
    --l-msjoin   INT
                 detect misjoined unitigs of >=INT in size; 0 to disable [500000]
  Ultra-Long-integration:
    --ul FILEs   file names of Ultra-Long reads [r1.fq,r2.fq,...]
    --ul-rate    FLOAT
                 error rate of Ultra-Long reads [0.2]
    --ul-tip     INT
                 remove tip unitigs composed of <=INT reads for the UL assembly [6]
    --path-max   FLOAT
                 max path drop ratio [0.6]; higher number may make the assembly cleaner
                 but may lead to more misassemblies
    --path-min   FLOAT
                 min path drop ratio [0.2]; higher number may make the assembly cleaner
                 but may lead to more misassemblies
    --ul-cut     INT
                 filter out <INT UL reads during the UL assembly [0]
  Dual-Scaffolding:
    --dual-scaf  output scaffolding
    --scaf-gap   INT
                 max gap size for scaffolding [3000000]
  Telomere-identification:
    --telo-m     STR
                 telomere motif at 5'-end; CCCTAA for human [NULL]
    --telo-p     INT
                 non-telomeric penalty [1]
    --telo-d     INT
                 max drop [2000]
    --telo-s     INT
                 min score for telomere reads [500]
Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz
See `https://hifiasm.readthedocs.io/en/latest/' or `man ./hifiasm.1' for complete documentation.
```
