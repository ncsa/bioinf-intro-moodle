#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=bcf_metrics
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=bcf_metrics.o%j
#SBATCH --error=bcf_metrics.e%j

set -e

# Input files
VCF=../results/sample1_bcftools_filtered.vcf.gz
BAM=../results/sample1_bqsr.bam
REF=../data/hg38_chr20_expanded.fa
STATS=../results/sample1_bcftools_filtered.stats
PLOT_DIR=../results/stats_plots

mkdir -p "$PLOT_DIR"

# 1. Basic sanity check
echo "Variant count:"
bcftools view -H "$VCF" | wc -l

echo "QUAL distribution:"
bcftools query -f '%QUAL\n' "$VCF" | sort -n | uniq -c | tail

# 2. Mean coverage over region
samtools depth -r 20:5000000-15000000 "$BAM" | \
awk '{sum+=$3} END {print "Mean depth:", sum/NR}'

# 3. VCF stats + plotting
bcftools stats "$VCF" > "$STATS"
plot-vcfstats -p "$PLOT_DIR" "$STATS"


