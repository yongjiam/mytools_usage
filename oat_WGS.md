## 

## modify chromosome position for part2.vcf.gz
```bash
#!/bin/bash --login

#SBATCH --job-name=bcftools
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE

#conda activate base
module load bcftools/1.15--haf5b3da_0
srun --export=all -n 1 -c 128 bcftools concat QUAL200MQ50DP6_more_HaplotypeCaller*.gz --threads 128 -Oz > updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz
srun --export=all -n 1 -c 128 bcftools index updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz --threads 128
srun --export=all -n 1 -c 128 plink --vcf updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz --make-bed --allow-extra-chr --double-id --out updated_QUAL200MQ50DP6_more_merged_snps
srun --export=all -n 1 -c 128 plink --bfile updated_QUAL200MQ50DP6_more_merged_snps --maf 0.05 --geno 0.2 --mind 0.2 --make-bed --out filtered_updated_QUAL200MQ50DP6_more_merged_snps --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile filtered_updated_QUAL200MQ50DP6_more_merged_snps --recode --out changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 sed -i 's/_part1//g;s/_part2//g' changeChr.map
srun --export=all -n 1 -c 128 plink --file changeChr --make-bed --out changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 plink2 --bfile changeChr --set-all-var-ids @:# --make-bed --out named_changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile named_changeChr --indep-pairwise 50 10 0.3 --out LD_pruned --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile named_changeChr --extract LD_pruned.prune.in --out LD_pruned --make-bed --allow-extra-chr
```

## for tree 
```bash
module load bcftools/1.15--haf5b3da_0
#srun --export=all -n 1 -c 128 bcftools concat QUAL200MQ50DP6_more_HaplotypeCaller*.gz --threads 128 -Oz > updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz
#srun --export=all -n 1 -c 128 bcftools index updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz --threads 128
#srun --export=all -n 1 -c 128 plink --vcf updated_QUAL200MQ50DP6_more_merged.snps.vcf.gz --make-bed --allow-extra-chr --double-id --out updated_QUAL200MQ50DP6_more_merged_snps
srun --export=all -n 1 -c 128 plink --bfile updated_QUAL200MQ50DP6_more_merged_snps --maf 0.05 --geno 0.01 --make-bed --out tree_filtered_updated_QUAL200MQ50DP6_more_merged_snps --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile tree_filtered_updated_QUAL200MQ50DP6_more_merged_snps --recode --out tree_changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 sed -i 's/_part1//g;s/_part2//g' tree_changeChr.map
srun --export=all -n 1 -c 128 plink --file tree_changeChr --make-bed --out tree_changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 plink2 --bfile tree_changeChr --set-all-var-ids @:# --make-bed --out tree_named_changeChr --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile tree_named_changeChr --indep-pairwise 50 10 0.01 --out tree_LD_pruned --allow-extra-chr
srun --export=all -n 1 -c 128 plink --bfile tree_named_changeChr --extract tree_LD_pruned.prune.in --out tree_LD_pruned --make-bed --allow-extra-chr
srun --export=all -n 1 -c 128 plink2 --bfile tree_LD_pruned --recode vcf id-paste=iid --out tree_LD_pruned --allow-extra-chr
```
