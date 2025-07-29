#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=hpcaller
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=hpcaller.o%j
#SBATCH --error=hpcaller.e%j

# Input
REF=../data/hg38_chr20_expanded.fa
INBAM=../results/sample1_bqsr.bam
OUTVCF=../results/sample1_raw.vcf.gz

# Run HaplotypeCaller
gatk HaplotypeCaller \
  -R $REF \
  -I $INBAM \
  -O $OUTVCF \
  -ERC NONE \
  -stand-call-conf 2.0 \
  --minimum-mapping-quality 0

# You can set -ERC to GVCF if you plan joint genotyping later

