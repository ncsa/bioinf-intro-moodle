#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=bcf_var_fil
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=bcf_var_fil.o%j
#SBATCH --error=bcf_var_fil.e%j

set -e

INVCF=../results/sample1_bcftools_raw.vcf.gz
OUTVCF=../results/sample1_bcftools_filtered.vcf.gz

bcftools filter \
  -e 'QUAL<30 || DP<10' \
  -s LOWQUAL \
  -Oz -o $OUTVCF \
  $INVCF

bcftools index $OUTVCF

