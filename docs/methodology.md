# Methodology

## Extended Haplotype Analyzer

The software extracts all variants within a user-defined genomic interval surrounding a Variant of Interest (VOI).

For each sample:

1. Variants are parsed from indexed VCF files.
2. Genotypes are converted into haplotype matrices.
3. Carrier haplotypes are identified.
4. Pairwise Identity-by-State (IBS) scores are computed.
5. Shared haplotypes are visualized.

IBS similarity is calculated as:

IBS = Matching Alleles / Comparable Alleles

excluding missing positions.

---

## Homozygosity Estimator

The software scans autosomal chromosomes and classifies each callable SNP as:

* Homozygous reference
* Homozygous alternate
* Heterozygous
* No-call

Genome-wide homozygosity is computed as:

Homozygosity (%) =
(HomRef + HomAlt) / Callable Sites × 100

ROH detection uses a sliding-window strategy based on:

* Minimum SNP count
* Minimum physical length
* Minimum homozygosity fraction

FROH is estimated as:

FROH = Total ROH Length / Total Autosomal Length

using GRCh38 autosomal chromosome sizes.
