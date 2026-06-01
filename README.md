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

Genome-wide homozygosity is estimated from callable SNP sites across autosomes (chr1–22), rather than from every genomic base.

The software identifies:

* Homozygous reference calls
* Homozygous alternate calls
* Heterozygous calls
* No-call sites

Important notes on genomic coverage:

* Analysis is restricted to autosomal chromosomes (chr1–22)
* Only high-quality SNPs are included (filtered by genotype quality and depth)
* Indels and multi-allelic variants are excluded by default
* Non-variant regions of the genome are not analyzed
* Sex chromosomes (chrX, chrY) and mitochondrial DNA are not included unless explicitly configured

Thus, the estimates reflect genome-wide patterns based on variant sites, rather than complete base-by-base genome coverage.

Runs of Homozygosity (ROH) are detected using a sliding-window approach and summarized as:

FROH = Total ROH Length / Total Autosomal Genome Length

FROH is commonly used as a genomic estimator of autozygosity and parental relatedness.

---

