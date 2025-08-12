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
OUTVCF=../results/sample1_bcftools_raw.vcf.gz

# Call variants
bcftools mpileup -f $REF $INBAM -Ou \
  | bcftools call -mv -Oz -o $OUTVCF

# Index and inspect
bcftools index $OUTVCF
bcftools view $OUTVCF | grep -v "^#" | wc -l


