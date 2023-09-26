## resources
https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/#:~:text=You%20can%20use%20FastQC%20to,FastQC%20on%20your%20fastq%20files.


## fastqc on fastq.gz
find /scratch/pawsey0399/yjia/skylar/WGS/ -name "*fq.gz"|grep trimmed| while read R;
do
	fastqc $R -o trimmed_fastqc_output -t 128  ## need to mkdir output dir in advance
done

## multiqc on fastqc outputs
conda activate bio
srun --export=all -n 1 -c 64 multiqc -i skylar_WGS_qc -n WGS_qc fastqc_output &> multiqc_log.txt
