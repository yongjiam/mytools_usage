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
## from 
```
## juicer_post.sh on nimbus
```
juicer post -o out_JBAT_review juicebox_input/out_JBAT.review.assembly juicebox_input/out_JBAT.liftover.agp /data/huyou/hifiasm_purged_hap1hap2merge/purged_join.fa
```
