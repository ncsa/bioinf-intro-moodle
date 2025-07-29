# Module 1 - DNA Variant Calling on Small Human Genome Region

This repository demonstrates a small-scale GATK single-sample variant calling pipeline using a ~1MB region of human chromosome 20.

## 📁 Folder Structure

After running `scripts/download_test_data.sh`, you will find a `data/` directory created. 

You will also need to manually create an empty `results` directory under `module-01-dna-variant-calling` with the command `mkdir results` 

module-01-dna-variant-calling/

├── data/ # Input FASTQs, reference genome, known variants

├── results/ # BAMs, VCFs, metrics

├── scripts/ # Step-by-step shell scripts

└── README.md


## 🔧 Scripts Overview

| Script                   | Description                             |
|--------------------------|-----------------------------------------|
| `01_align.sh`            | Align FASTQs to reference and sort      |
| `02_markdup.sh`          | Mark PCR duplicates and index BAM       |
| `03_bqsr.sh`             | Recalibrate base quality scores         |
| `04_haplotypecaller.sh`  | Call raw variants with HaplotypeCaller  |
| `05_variant_filtering.sh`| Apply hard filters to raw VCF           |
| `06_metrics.sh`          | Collect QC metrics                      |

## 🧪 Data Requirements

- Small human genome region (e.g., chr20:10Mb–11Mb)
- Paired-end FASTQ files (~100k reads each)
- Reference FASTA (`hg38_chr20.fa`) + `.fai` and `.dict`
- Known sites VCF (e.g., dbSNP limited to chr20)

## 💻 How to Submit Jobs  

For each script, use a text editor such as vim or nano to put in the account name and queue.

```bash
cd scripts/
sbatch 01_align.sh
sbatch 02_markdup.sh
sbatch 03_bqsr.sh
sbatch 04_haplotypecaller.sh
sbatch 05_variant_filtering.sh
sbatch 06_metrics.sh
```

## bcftools vs GATK

our example data was small. therefore, 05_variant_filtering and 06_metrics.sh won't discover anything. Therefore, we deployed bcftools instead of GATK for these steps
