#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=align
#SBATCH --account=      # the account to charge
#SBATCH --partition=    # which queue are we using
#SBATCH --output=align.o%j
#SBATCH --error=align.e%j 

# Input
REF=../data/hg38_chr20_expanded.fa
FASTQ1=../data/sample_R1_trimmed.fastq.gz
FASTQ2=../data/sample_R2_trimmed.fastq.gz
SAMPLE=sample1

# Output
OUTBAM=../results/${SAMPLE}_sorted.bam

# Index the reference if needed
bwa index $REF

# Align and sort
bwa mem -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" $REF $FASTQ1 $FASTQ2 | \
  samtools sort -o $OUTBAM

samtools index $OUTBAM

