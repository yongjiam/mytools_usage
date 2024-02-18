## deepvariant for PACBIO HIFI variant call
https://github.com/google/deepvariant
### singularity image
```bash
singularity build deepvariant.sif docker://google/deepvariant:1.6.0
```
### read alignment minimap2

```bash
 ## install minimap2 https://github.com/lh3/minimap2
conda install bioconda::minimap2

## or
git clone https://github.com/lh3/minimap2
cd minimap2 && make

## long read alignment
# long sequences against a reference genome
./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam
# create an index first and then map
./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa
./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam
# use presets (no test data)
./minimap2 -ax map-pb ref.fa pacbio.fq.gz > aln.sam       # PacBio CLR genomic reads
./minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads
./minimap2 -ax map-hifi ref.fa pacbio-ccs.fq.gz > aln.sam # PacBio HiFi/CCS genomic reads (v2.19 or later)
./minimap2 -ax asm20 ref.fa pacbio-ccs.fq.gz > aln.sam    # PacBio HiFi/CCS genomic reads (v2.18 or earlier)
./minimap2 -ax sr ref.fa read1.fa read2.fa > aln.sam      # short genomic paired-end reads
./minimap2 -ax splice ref.fa rna-reads.fa > aln.sam       # spliced long reads (strand unknown)
./minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # noisy Nanopore Direct RNA-seq
./minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # Final PacBio Iso-seq or traditional cDNA
./minimap2 -ax splice --junc-bed anno.bed12 ref.fa query.fa > aln.sam  # prioritize on annotated junctions
./minimap2 -cx asm5 asm1.fa asm2.fa > aln.paf             # intra-species asm-to-asm alignment
./minimap2 -x ava-pb reads.fa reads.fa > overlaps.paf     # PacBio read overlap
./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf    # Nanopore read overlap

# man page for detailed command line options
man ./minimap2.1

```
### run variant call
```bash
##deepvariant.sh
mkdir tmp
singularity exec --bind ${PWD}:${PWD},tmp:/tmp deepvariant.sif /opt/deepvariant/bin/run_deepvariant \
	    --model_type PACBIO \
	    --ref Clipper.V1.fasta \
	    --reads sorted.Clipper-Buloke.markdup.bam \
	    --output_vcf deepvariant_output.vcf.gz \ ## required
	    --output_gvcf deepvariant_output.gvcf.gz \
	    --num_shards 128
```
