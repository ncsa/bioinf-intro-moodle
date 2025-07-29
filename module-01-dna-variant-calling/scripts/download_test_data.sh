#!/bin/bash
set -e
mkdir -p ../data && cd ../data

# Download full GRCh38 reference from Ensembl (primary assembly)
#wget -O GRCh38_full.fa.gz ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gunzip GRCh38_full.fa.gz

# Extract chr20 only (adjust name if needed)
samtools faidx GRCh38_full.fa 20 > hg38_chr20.fa

# Index and dictionary
samtools faidx hg38_chr20.fa
gatk CreateSequenceDictionary -R hg38_chr20.fa

# Download dbSNP VCF (chr20 only)
wget -O dbsnp_chr20.vcf.gz https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
tabix -p vcf dbsnp_chr20.vcf.gz || echo "You may need to subset chr20 manually for testing."

# Match headings
# treat dbsnp files to match the chr20 headings
bcftools view -r NC_000020.11 dbsnp_chr20.vcf.gz -Oz -o dbsnp_chr20_subset.vcf.gz

# Create new header with corrected contig
bcftools view -h dbsnp_chr20_subset.vcf.gz | sed 's/NC_000020.11/20/g' > new_header.txt

# Reheader VCF
bcftools reheader -h new_header.txt -o dbsnp_chr20_fixed.vcf.gz dbsnp_chr20_subset.vcf.gz
tabix -p vcf dbsnp_chr20_fixed.vcf.gz


# Download sample FASTQ from 1000 Genomes (HG00258, small read subset)
wget -O sample_R1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -O sample_R2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# Trim the input data to keep samples small
gzcat sample_R1.fastq.gz | head -n 2000000 | gzip > sample_R1_trimmed.fastq.gz
gzcat sample_R2.fastq.gz | head -n 2000000 | gzip > sample_R2_trimmed.fastq.gz


echo "All test data downloaded to ../data/"
