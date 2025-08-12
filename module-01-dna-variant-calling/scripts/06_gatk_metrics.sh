#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=gatk_metrics
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=gatk_metrics.o%j
#SBATCH --error=gatk_metrics.e%j

# Input
REF=../data/hg38_chr20.fa
INBAM=../results/sample1_bqsr.bam
VCF=../results/sample1_gatk_filtered.vcf.gz
SAMPLE=sample1

# Output
ALIGN_METRICS=../results/${SAMPLE}_alignment_metrics.txt
VARIANT_METRICS=../results/${SAMPLE}_variant_metrics.txt

# 1. Basic sanity check
echo "Variant count:"
bcftools view -H "$VCF" | wc -l

echo "QUAL distribution:"
bcftools query -f '%QUAL\n' "$VCF" | sort -n | uniq -c | tail

# 2. Mean coverage over region
samtools depth -r 20:5000000-15000000 "$BAM" | \
awk '{sum+=$3} END {print "Mean depth:", sum/NR}'

# 3. Collect alignment summary
gatk CollectAlignmentSummaryMetrics \
  R=$REF \
  I=$INBAM \
  O=$ALIGN_METRICS

# 4. Collect variant calling summary
gatk CollectVariantCallingMetrics \
  INPUT=$VCF \
  OUTPUT=$VARIANT_METRICS \
  DBSNP=../data/dbsnp_chr20.vcf.gz \
  SEQUENCE_DICTIONARY=../data/hg38_chr20.dict

