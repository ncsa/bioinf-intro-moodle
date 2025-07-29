#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=markdup
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=markdup.o%j
#SBATCH --error=markdup.e%j


# Input
INBAM=../results/sample1_sorted.bam
OUTBAM=../results/sample1_marked.bam
METRICS=../results/sample1_markdup_metrics.txt

# Run MarkDuplicates
gatk MarkDuplicates \
  -I $INBAM \
  -O $OUTBAM \
  -M $METRICS

# Index the output BAM
samtools index $OUTBAM

