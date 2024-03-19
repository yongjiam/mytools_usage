## correct the chromosome name for hap1 and hap2
```bash
## sort by chromosome lengths
seqkit sort -lr hap1.fasta > tmp && mv tmp hap1.fasta
seqkit sort -lr hap2.fasta > tmp && mv tmp hap2.fasta

## long to short, named as chr1-9
### sed_id_hap1
s/h1tg000012l/H1_ch1/
s/h1tg000002l/H1_ch2/
s/h1tg000001l/H1_ch3/
s/h1tg000004l/H1_ch4/
s/h1tg000017l/H1_ch5/
s/h1tg000010l/H1_ch6/
s/h1tg000032l/H1_ch7/
s/h1tg000027l/H1_ch8/
s/h1tg000029l/H1_ch9/
### sed_id_hap2
s/h2tg000011l/H2_ch1/
s/h2tg000005l/H2_ch2/
s/h2tg000001l/H2_ch3/
s/h2tg000017l/H2_ch5/
s/h2tg000003l/H2_ch4/
s/h2tg000009l/H2_ch6/
s/h2tg000031l/H2_ch7/
s/h2tg000026l/H2_ch8/
s/h2tg000027l/H2_ch9/

## rename chromosome using other citrus species as references
##### Citrus sinensis v3.0
#####

## use jcvi to match chromosome IDs
## run liftoff
/data/huyou/liftoff
liftoff -g huyou.gff -o hap1.liftoff.gff3 \
	-u hap1.unmapped \
	-copies \ ## look for extra copies
	-sc 0.95 \ ## copy identify threshold
	-p 15 \ ##paralel processes
  huyou.hap1.genome.fa huyou.genome.fa ## target and reference

## hap1.bed
python -m jcvi.formats.gff bed --type=mRNA --primary_only hap1.liftoff.gff3 -o hap1.bed
sed -i -f sed_id_hap1 hap1.bed ## change chromosome id
sed -i 's/Ccha/H1_Ccha/' hap1.bed ## change gene id

## hap1.cds
gffread -x hap1.cds -g huyou.hap1.genome.fa hap1.liftoff.gff3 -F
sed -i 's/Ccha/H1_Ccha/' hap1.cds

## pairwise jcvi with SWO
##### dotplot
python -m jcvi.compara.catalog ortholog SWO hap1 --cscore=.99 --no_strip_names ## (gene id starts with Cs_ont_*.2+ not found in bed)
python -m jcvi.graphics.dotplot SWO.hap1.anchors

## produce synteny
python -m jcvi.compara.synteny screen --minspan=30 --simple SWO.hap1.anchors SWO.hap1.anchors.new

## prepare seqids and layout file
#### seqids
	chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chrUn
	H1_ch1,H1_ch2,H1_ch3,H1_ch4,H1_ch5,H1_ch6,H1_ch7,H1_ch8,H1_ch9
#### layout
	# y, xstart, xend, rotation, color, label, va,  bed
	 .7,     .1,    .8,      15,      , SWO, top, SWO.bed
	 .5,     .1,    .8,       0,      , hap1, top, hap1.bed
	# edges
	e, 0, 1, SWO.hap1.anchors.simple
## produce synteny plot
python -m jcvi.compara.synteny screen --minspan=30 --simple SWO.hap1.anchors SWO.hap1.anchors.new
python -m jcvi.graphics.karyotype seqids layout
mv karyotype.pdf karyotype_SWO.hap1.pdf
```
