#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=bcf_metrics
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=bcf_metrics.o%j
#SBATCH --error=bcf_metrics.e%j

set -e

bcftools stats ../results/sample1_bcftools_filtered.vcf.gz > ../results/sample1_bcftools_filtered.stats

plot-vcfstats -p ../results/stats_plots ../results/sample1_bcftools_filtered.stats

