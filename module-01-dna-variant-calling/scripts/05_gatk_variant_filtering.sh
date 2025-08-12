#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=gatk_var_fil
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=gatk_var_fil.o%j
#SBATCH --error=gatk_var_fil.e%j

# Input
INVCF=../results/sample1_gatk_raw.vcf.gz
OUTVCF=../results/sample1_gatk_filtered.vcf.gz

# Apply hard filters to SNPs
gatk VariantFiltration \
  -V $INVCF \
  -O $OUTVCF \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "FS_gt_60" --filter-expression "FS > 60.0" \
  --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0"

