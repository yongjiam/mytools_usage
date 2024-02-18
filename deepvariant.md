## deepvariant for PACBIO HIFI variant call
### singularity image
```bash
singularity build deepvariant.sif docker://google/deepvariant:1.6.0
```
### run variant call
```bash
##deepvariant.sh
singularity exec --bind ${PWD}:${PWD} deepvariant.sif /opt/deepvariant/bin/run_deepvariant \
	    --model_type PACBIO \
	    --ref Clipper.V1.fasta \
	    --reads sorted.Clipper-Buloke.markdup.bam \
	    --output_vcf deepvariant_output.vcf.gz \ ## required
	    --output_gvcf deepvariant_output.gvcf.gz \
	    --num_shards 128
```
