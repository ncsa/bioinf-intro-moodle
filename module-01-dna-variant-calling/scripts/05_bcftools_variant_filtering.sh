#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=bcf_var_fil
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=bcf_var_fil.o%j
#SBATCH --error=bcf_var_fil.e%j

set -e

bcftools filter \
  -e 'QUAL<30 || DP<10' \
  -s LOWQUAL \
  -Oz -o ../results/sample1_bcftools_filtered.vcf.gz \
  ../results/sample1_bcftools.vcf.gz

bcftools index ../results/sample1_bcftools_filtered.vcf.gz

