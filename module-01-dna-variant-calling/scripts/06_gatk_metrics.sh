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
VCF=../results/sample1_filtered.vcf.gz
SAMPLE=sample1

# Output
ALIGN_METRICS=../results/${SAMPLE}_alignment_metrics.txt
VARIANT_METRICS=../results/${SAMPLE}_variant_metrics.txt

# Collect alignment summary
gatk CollectAlignmentSummaryMetrics \
  R=$REF \
  I=$INBAM \
  O=$ALIGN_METRICS

# Collect variant calling summary
gatk CollectVariantCallingMetrics \
  INPUT=$VCF \
  OUTPUT=$VARIANT_METRICS \
  DBSNP=../data/dbsnp_chr20.vcf.gz \
  SEQUENCE_DICTIONARY=../data/hg38_chr20.dict

