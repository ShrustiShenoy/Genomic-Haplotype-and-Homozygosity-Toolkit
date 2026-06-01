# Genomic Haplotype and Homozygosity Toolkit

A bioinformatics toolkit for investigating founder effects, shared haplotypes, genome-wide homozygosity, and runs of homozygosity (ROH) from single-sample VCF datasets.

The toolkit contains two complementary analysis modules:

1. **Extended Haplotype Analyzer**
2. **Genome-Wide Homozygosity Estimator**

These tools are designed for rare disease studies, founder variant investigations, population genetics analyses, and genomic ancestry assessment using hard-filtered VCF files.

---

## Features

### Extended Haplotype Analyzer

* Variant-of-interest (VOI) centered haplotype analysis
* Automatic chromosome naming detection
* Haplotype carrier identification
* Identity-by-State (IBS) similarity estimation
* Founder haplotype assessment
* Publication-ready visualizations
* Compatible with WGS and low-coverage genome datasets

Outputs include:

* Haplotype block grid
* IBS similarity heatmap
* Alternate allele frequency plots
* VOI summary panels
* TSV genotype reports

---

### Genome-Wide Homozygosity Estimator

* Genome-wide homozygosity estimation
* Heterozygosity assessment
* Runs of Homozygosity (ROH) detection
* F-statistic estimation
* FROH calculation
* Per-chromosome homozygosity profiling
* Publication-quality multi-panel reports

Outputs include:

* Homozygosity statistics
* ROH segments
* Chromosome-level homozygosity heatmaps
* Genome ideograms
* Summary tables

---

## Scientific Background

### Haplotype Analysis

Individuals carrying the same pathogenic variant may share surrounding genomic markers inherited from a common ancestor. By comparing phased and unphased haplotypes surrounding a Variant of Interest (VOI), the software evaluates potential founder effects and shared ancestry.

The toolkit calculates pairwise Identity-by-State (IBS) similarity among haplotypes:

* IBS ≥ 0.95 → likely shared founder haplotype
* IBS 0.70–0.95 → shared ancestry
* IBS < 0.70 → independent haplotypes

### Homozygosity Analysis

Genome-wide homozygosity is estimated from callable SNP sites.

The software identifies:

* Homozygous reference calls
* Homozygous alternate calls
* Heterozygous calls
* No-call sites

Runs of Homozygosity (ROH) are detected using a sliding-window approach and summarized as:

FROH = Total ROH Length / Total Autosomal Genome Length

FROH is commonly used as a genomic estimator of autozygosity and parental relatedness.

---

## Requirements

Python 3.9+

Required packages:

```bash
pip install cyvcf2 pandas numpy matplotlib scipy configparser
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/USERNAME/Genomic-Haplotype-and-Homozygosity-Toolkit.git

cd Genomic-Haplotype-and-Homozygosity-Toolkit
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Running Haplotype Analysis

Configure:

```ini
[PATHS]
vcf_dir=/path/to/vcfs
results_dir=./results

[VCF_FILES]
sample1=sample1
sample2=sample2
```

Run:

```bash
python haplotype_analysis_50kb.py config.ini
```

---

## Running Homozygosity Analysis

Configure:

```ini
[SAMPLES]
sample1=/path/sample1.vcf.gz
sample2=/path/sample2.vcf.gz
```

Run:

```bash
python homozygosity_visualiser2.py homozygosity_config.in
```

---

## Input Requirements

Each sample must have:

```text
sample.hard-filtered.vcf.gz
sample.hard-filtered.vcf.gz.tbi
```

Generated using bgzip and tabix indexing.

---

## Applications

* Founder variant investigations
* Rare disease studies
* Population genetics
* Carrier screening research
* Consanguinity assessment
* Autozygosity studies
* Clinical genomics research

---

## Disclaimer

This software is intended for research use only and should not be used as the sole basis for clinical decision-making.

---

## Citation

If you use this software in academic work, please cite the repository and include the software version used in your analyses.

