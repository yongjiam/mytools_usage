### data
1. illunina
   -rw-r--r-- 1 yjia pawsey0399  21G Mar  2 11:37 1_R2.fq.gz\
   -rw-r--r-- 1 yjia pawsey0399  19G Mar  2 11:37 1_R1.fq.gz\
   -rw-r--r-- 1 yjia pawsey0399 9.6K Mar  2 11:37 1.quality.png\
   -rw-r--r-- 1 yjia pawsey0399 6.4K Mar  2 11:37 1.quality.pdf\
   -rw-r--r-- 1 yjia pawsey0399  26K Mar  2 11:37 1.base.png\
   -rw-r--r-- 1 yjia pawsey0399  21K Mar  2 11:37 1.base.pdf\
   
3. hifi
   -rw-r--r-- 1 yjia pawsey0399 1.5K Mar  2 11:37 m64257e_211030_130656.ccs.subreadset.xml\
   -rw-r--r-- 1 yjia pawsey0399  17M Mar  2 11:37 m64257e_211030_130656.ccs.bam.pbi\
   -rw-r--r-- 1 yjia pawsey0399   64 Mar  2 11:37 m64257e_211030_130656.ccs.bam.md5\
   -rw-r--r-- 1 yjia pawsey0399  21G Mar  2 11:37 m64257e_211030_130656.ccs.bam\
3. HIC
   -rw-r--r-- 1 yjia pawsey0399 7.7G Mar  2 11:37 changshanhuyou-1_R2.fq.gz\
   -rw-r--r-- 1 yjia pawsey0399 7.3G Mar  2 11:37 changshanhuyou-1_R1.fq.gz\
   -rw-r--r-- 1 yjia pawsey0399  10K Mar  2 11:37 changshanhuyou-1.quality.png\
   -rw-r--r-- 1 yjia pawsey0399 6.5K Mar  2 11:37 changshanhuyou-1.quality.pdf\
   -rw-r--r-- 1 yjia pawsey0399  36K Mar  2 11:37 changshanhuyou-1.base.png\
   -rw-r--r-- 1 yjia pawsey0399  22K Mar  2 11:37 changshanhuyou-1.base.pdf\
4. RNAseq
   -rwx------ 1 ubuntu ubuntu 2.2G Mar 17 05:00 R1.fq.gz*
   -rwx------ 1 ubuntu ubuntu 2.2G Mar 17 05:01 R2.fq.gz*
## genome assembly
```bash
#### huyou.conf
#!/bin/bash --login

#SBATCH --job-name=huyou
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --mem=980G
#SBATCH --export=NONE

module load samtools/1.15--h3843a85_0
module load singularity/3.11.4-slurm
#srun --export=all -n 1 -c 128  samtools fastq -@ 128 ./ccs/m64257e_211030_130656.ccs.bam > hifi_ccs.fastq
srun --export=all -n 1 -c 64 singularity exec --bind ${PWD}:${PWD} hifiasm_latest.sif hifiasm -o huyou.asm -t 64 \
	--h1 ./HIC/changshanhuyou-1_R1.fq.gz \
	--h2 ./HIC/changshanhuyou-1_R1.fq.gz \
	hifi_ccs.fastq
```
## survey using genomescope
#### count chromosome and contig number for 21 genome data download from http://citrus.hzau.edu.cn/download.php
```bash
ls *.gff3|while read R;do echo $(echo $R|cut -d '.' -f1); echo $(awk '!/^#/{print $1}' $R|sort|uniq|wc -l);done |paste - - > chromosome_count

## chromosome_count
CMJ	10	Citrus_grandis_majiayou
GJ	10	Fortunella_hindsii
HWB	10	Citrus_grandis_wanbaiyou
HZYT	10	Citrus_maxima_huazhouyou
SWO	10	Citrus_sinensis
ZK	10	Poncirus_trifoliata
GCF	94	Citrus_clementina
AEG	138	Aegle_marmelos
MSYG	160	Citrus_mangshanensis
CGI	170	Citropsis_gilletiana
ZGYCC	205	Citrus_ichangensis
HKC	221	Atalantia_buxfoliata
JLX	253	Murraya_paniculata
LW	270	luvunga_scandens
AZM	331	Citrus_australasica
HP	388	Clausena_lansium
XGF	493	Citrus_maxima_majia
HH	501	Citrus_hongheensis
SYT	1001	Luvunga_scandens
JZ	4141	Citrus_reticulata
RL	4465	Citrus_medica
```
#### count chromosome length from gff3 file
```bash
awk '$0 !~ /^#/ {chromosome[$1]=$5} END {for (chr in chromosome) print "Chromosome", chr, ": Total length =", chromosome[chr]}' ./pangenomes/Citrus_changshan-huyou.gene.gff

```
## orthofinder for 21 citrus species
```bash
## get primary sequence from fasta file
for i in *.protein.fa
do
	VAR=$(echo $i|cut -d "." -f1-2)
	bioawk -c fastx '{if ($name ~ /\.1$/) print ">"$name "\n" $seq}' $i > "primary_"$VAR".protein.fasta"
done

## run orthofinder /data/huyou/orthofinder/
orthofinder -f selected
```
## genome model prediction using gemoma (company annotation does not match genome.fa)
reference genomes from phytozome (citrus database data throw errors in gemoma)
```bash
##install gemoma 1.9
conda install -c bioconda gemoma

## run gemoma.sh in setonix
## genome level prediction
GeMoMa GeMoMaPipeline threads=64 tblastn=False \
	AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=HY \
	p=false \
	o=true \
	t=./huyou.hap1.genome.fa \
	outdir=hap1/ \
	s=own i=Ccl a=./phytozome/Cclementina_182_v1.0.gene.gff3.gz g=./phytozome/Cclementina_182_v1.fa.gz \
	s=own i=Csi a=./phytozome/Csinensis_154_v1.1.gene.gff3.gz g=./phytozome/Csinensis_154_v1.fa.gz \
	s=own i=Ptr a=./phytozome/Ptrifoliata_565_v1.3.1.gene.gff3.gz g=./phytozome/Ptrifoliata_565_v1.3.fa.gz
```
## genome model prediction using liftoff (using company annotation)
reference genomes from phytozome (citrus database data throw errors in gemoma)
```bash
##install gemoma 1.9
conda install -c bioconda gemoma

## run gemoma.sh in setonix
## genome level prediction
GeMoMa GeMoMaPipeline threads=64 tblastn=False \
	AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=HY \
	p=false \
	o=true \
	t=./huyou.hap1.genome.fa \
	outdir=hap1/ \
	s=own i=Ccl a=./phytozome/Cclementina_182_v1.0.gene.gff3.gz g=./phytozome/Cclementina_182_v1.fa.gz \
	s=own i=Csi a=./phytozome/Csinensis_154_v1.1.gene.gff3.gz g=./phytozome/Csinensis_154_v1.fa.gz \
	s=own i=Ptr a=./phytozome/Ptrifoliata_565_v1.3.1.gene.gff3.gz g=./phytozome/Ptrifoliata_565_v1.3.fa.gz
```

## synteny dotplot
1. gene-based
   mcscan
   
   jcvi
   https://github.com/tanghaibao/jcvi
   ```bash
   ## install
   pip install jcvi

   ## gff2bed
   ls *.gff*|while read R;do VAR=$(echo $R|cut -d '.' -f1);python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only $R -o $VAR".bed";done
   ## extract cds
   ls *.gff*|while read R;do V=$(echo $R|cut -d '.' -f1);gffread -x $V".cds" -g $V*.genome.fa $R -F;done
   ls *.cds|while read R;do python -m jcvi.formats.fasta format $R "formated_"$R;done

   ##pairwise synteny
   python -m jcvi.compara.catalog ortholog SWO HWB --no_strip_names ##produce anchor file and dotplot by default
   python -m jcvi.compara.catalog ortholog SWO HWB --cscore=.99 --no_strip_names ## 1:1 orthologous region only
   #or
   python -m jcvi.graphics.dotplot SWO.HWB.anchors

   ##multiple species
   python -m jcvi.compara.catalog ortholog huyou SWO --cscore=.99 --no_strip_names
   python -m jcvi.compara.synteny screen --minspan=30 --simple huyou.SWO.anchors huyou.SWO.anchors.new
   
   python -m jcvi.compara.catalog ortholog huyou HWB --cscore=.99 --no_strip_names
   python -m jcvi.compara.synteny screen --minspan=30 --simple huyou.HWB.anchors huyou.HWB.anchors.new

   ## create layout
	   # y, xstart, xend, rotation, color, label, va,  bed
	 .7,     .1,    .8,      15,      , SWO, top, SWO.bed
	 .5,     .1,    .8,       0,      , huyou, top, huyou.bed
	 .3,     .1,    .8,     -15,      , HWB, bottom, HWB.bed
	# edges
	e, 0, 1, huyou.SWO.anchors.simple
	e, 1, 2, huyou.HWB.anchors.simple

   ## modify huyou chr id
	s/chr1/ch5/
	s/chr2/ch3/
	s/chr3/ch2/
	s/chr4/ch8/
	s/chr5/ch9/
	s/chr6/ch7/
	s/chr7/ch4/
	s/chr8/ch1/
	s/chr9/ch6/
   sed -i -f sed_huyou_chr huyou.bed

   ## create seqids
   cut -f1 SWO.bed |sort|uniq|tr '\n' ',' |sed 's/,$/\n/' > seqids
   cut -f1 huyou.bed |sort |uniq|grep chr|tr '\n' ','|sed 's/,$/\n/' >> seqids
   cut -f1 HWB.bed|sort|uniq|tr '\n' ',' |sed 's/,$/\n/' >> seqids

   ## create plot
   python -m jcvi.graphics.karyotype seqids layout   
   ```
2. genome-based
   mummer
   https://www.nature.com/articles/s41588-022-01015-0 \
   ```bash
   git clone https://github.com/mummer4/mummer
   cd mummer && autoreconf -fi
   ./configure --prefix=/data/tools/mybin
   ```
   minimap2
   https://github.com/lh3/minimap2
   ```bash
   ## installation
   git clone https://github.com/lh3/minimap2
   cd minimap2 && make

   ## align two reference genome
   minimap2 -cx asm5 hap1.genome.fa hap2.genome.fa > minimap_aln.paf  # intra-species asm-to-asm alignment

   ## plot the results in R
   ```R
   # install packages
   install.packages("pafr")
   install.packages("patchwork")
   install.packages("tidyverse")
   
   # Load necessary libraries
   library(pafr)      # For handling PAF (Pairwise mApping Format) files
   library(patchwork) # For creating composite plots
   library(tidyverse) # For data manipulation and visualization
   # Read the PAF file into a data frame
   df <- read_paf("minimap_aln.paf")
   # Display the first few rows of the data frame
   df %>% as.data.frame() %>% head()

   ####### Create a dotplot visualization from the PAF data
   dotplot(df, order_by='provided',
           ordering= list(c("H2_ch1","H2_ch2","H2_ch3","H2_ch4","H2_ch5","H2_ch6","H2_ch7","H2_ch8","H2_ch9"),
           c("H1_ch1","H1_ch2","H1_ch3","H1_ch4","H1_ch5","H1_ch6","H1_ch7","H1_ch8","H1_ch9")),
           label_seqs = TRUE)
   
   ####### plot each chromosome separately in subplots
   install.packages("gridExtra")
   library(gridExtra)
   # Create an empty list to store the ggplot objects
   plots_list <- list()
   # Loop from 1 to 9
   for (i in 1:9) {
     # Create ggplot for each iteration
       CHR <- paste0("ch", i)
       CHR1 <- paste0("H1_ch", i)
       CHR2 <- paste0("H2_ch", i)
       plot <- dotplot(df, order_by='provided',
               ordering= list(c(CHR2),
               c(CHR1)),
               label_seqs = FALSE) +
       scale_x_continuous(breaks = c(0, 10000000, 20000000, 30000000, 40000000, 45000000), labels = c("0", "10", "20", "30", "40", "45"))+
       scale_y_continuous(breaks = c(0, 10000000, 20000000, 30000000, 40000000, 45000000), labels = c("0", "10", "20", "30", "40", "45"))+
       labs(x = CHR) +
       labs(y = "")
       
     # Name the ggplot object
     plot_name <- paste0("P", i)
     
     # Store the ggplot object in the list
     plots_list[[plot_name]] <- plot
   }
   
   # Arrange plots in a 3x3 grid
   final_plot <- grid.arrange(grobs = plots_list, nrow = 3, ncol = 3)
   ggsave("final_plot.pdf", final_plot, width = 6, height = 6)

   ######### Create a dotplot visualization for a single chromosome from the PAF data
   dotplot(df, order_by='provided',
           ordering= list(c("H2_ch2"),
           c("H1_ch2")),
           label_seqs = TRUE)
   # Create a synteny plot for a specific chromosome pair
   # (Adjust 'q_chrom' and 't_chrom' to match the desired chromosomes)
   plot2 <- plot_synteny(df, q_chrom = "H2_ch2", t_chrom = "H1_ch2", centre = TRUE)
   ggsave("synteny_chr2_plot.pdf", plot2, width = 18, height = 6)

   #############Create a coverage plot, filling based on the query sequences
   plot_coverage(df, fill = "qname")
   
   ############
   ```
   
