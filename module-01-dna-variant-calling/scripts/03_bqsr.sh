#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=bqsr
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=bqsr.o%j
#SBATCH --error=bqsr.e%j

# Input
REF=../data/hg38_chr20_expanded.fa
INBAM=../results/sample1_marked.bam
DBSNP=../data/dbsnp_chr20_overlap.vcf.gz          # Known variant sites for BQSR
RECAL=../results/sample1_recal.table
POST_RECAL=../results/sample1_post_recal.table
OUTBAM=../results/sample1_bqsr.bam

# Step 1: Generate recalibration table
gatk BaseRecalibrator \
  -R $REF \
  -I $INBAM \
  --known-sites $DBSNP \
  -O $RECAL

# Step 2: Apply the recalibration
gatk ApplyBQSR \
  -R $REF \
  -I $INBAM \
  --bqsr-recal-file $RECAL \
  -O $OUTBAM

# Optional: Generate post-recalibration table for comparison
gatk BaseRecalibrator \
  -R $REF \
  -I $OUTBAM \
  --known-sites $DBSNP \
  -O $POST_RECAL

