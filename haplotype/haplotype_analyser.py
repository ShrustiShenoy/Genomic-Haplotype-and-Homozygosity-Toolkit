#!/usr/bin/env python3
"""
================================================================================
  GENOME-WIDE HOMOZYGOSITY ESTIMATOR
  
  Estimates percentage of homozygosity for WGS/LCG samples from hard-filtered
  VCF files. Produces:
    ✓ Global homozygosity statistics (per sample)
    ✓ Runs of Homozygosity (ROH) detection via sliding window
    ✓ Per-chromosome homozygosity breakdown
    ✓ Zygosity class breakdown (HomRef / HomAlt / Het / No-call)
    ✓ F-statistic (inbreeding coefficient estimate)
    ✓ Publication-quality multi-panel figure (white background, Nature style)
    ✓ TSV summary table

  Compatible with:
    - Hard-filtered single-sample VCFs (.vcf.gz + .tbi)
    - WGS and LCG (low-coverage genome) data
    - Samples with or without phasing

  Dependencies:
    pip install cyvcf2 matplotlib pandas numpy scipy

  Usage:
    python3 homozygosity_estimator.py          [uses defaults below]
    python3 homozygosity_estimator.py config.ini
================================================================================
"""

# ── Standard Library ──────────────────────────────────────────────────────────
import os
import sys
import warnings
import argparse
import configparser
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings("ignore")

# ── Third-party ───────────────────────────────────────────────────────────────
try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    from matplotlib.gridspec import GridSpec
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.patheffects as pe
    from cyvcf2 import VCF
except ImportError as e:
    sys.exit(
        f"\n[ERROR] Missing dependency: {e}\n"
        "Install with:\n"
        "  pip install cyvcf2 matplotlib pandas numpy\n"
    )

# ══════════════════════════════════════════════════════════════════════════════
#  MATPLOTLIB GLOBAL STYLE  — publication-ready, white background
# ══════════════════════════════════════════════════════════════════════════════
matplotlib.rcParams.update({
    "font.family":        "DejaVu Sans",
    "font.size":          10,
    "axes.titlesize":     12,
    "axes.labelsize":     10,
    "axes.linewidth":     0.8,
    "axes.edgecolor":     "#333333",
    "axes.facecolor":     "white",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.major.size":   3.5,
    "ytick.major.size":   3.5,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "grid.color":         "#e5e5e5",
    "grid.linewidth":     0.6,
    "legend.frameon":     True,
    "legend.framealpha":  0.9,
    "legend.edgecolor":   "#cccccc",
    "figure.facecolor":   "white",
    "savefig.facecolor":  "white",
    "pdf.fonttype":       42,   # editable text in PDF/Illustrator
    "ps.fonttype":        42,
})

# ══════════════════════════════════════════════════════════════════════════════
#  PUBLICATION COLOUR PALETTE
#  Inspired by Nature/Cell figure standards:
#    - Colourblind-safe (Okabe-Ito palette base)
#    - High contrast on white
#    - Distinct hues per zygosity class
# ══════════════════════════════════════════════════════════════════════════════
PUB_COLORS = {
    # Zygosity classes
    "hom_ref":  "#4878CF",   # steel blue
    "hom_alt":  "#D65F5F",   # muted red
    "het":      "#6ACC65",   # muted green
    "nocall":   "#B8B8B8",   # neutral grey

    # ROH
    "roh":      "#E8A838",   # amber (warm, distinct from blue/red/green)

    # Per-sample accent colours (up to 8 samples)
    # Okabe-Ito colourblind-safe palette
    "samples": [
        "#4878CF",   # blue
        "#D65F5F",   # red
        "#6ACC65",   # green
        "#B47CC7",   # purple
        "#C4AD66",   # tan/gold
        "#77BEDB",   # light blue
        "#D6A76A",   # orange-tan
        "#E87070",   # salmon
    ],

    # Chromosome heatmap — white → orange → dark red
    "heatmap_cmap": ["#FFFFFF", "#FDE8C8", "#F4A460", "#C1440E", "#7B1600"],
}

# ══════════════════════════════════════════════════════════════════════════════
#  DEFAULT PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
DEFAULTS = {
    "samples": {
        "1326262332":    "/mnt/disk2/Shrusti/haplotype_check/vcfs/1326262332.hard-filtered.vcf.gz",
        "1326003255_LPG": "/mnt/disk2/Shrusti/haplotype_check/vcfs/1326003255_LPG.hard-filtered.vcf.gz",
    },
    "results_dir":      "./homozygosity_results",
    "chromosomes":      [f"chr{i}" for i in range(1, 23)],
    "numeric_chroms":   False,
    "min_gq":           20,
    "min_dp":           8,
    "skip_multiallelic": True,
    "skip_indels":      True,
    "roh_window_snps":  50,
    "roh_min_snps":     50,
    "roh_min_hom_frac": 0.95,
    "roh_min_length_bp": 500_000,
    "output_dpi":       300,
}

# ══════════════════════════════════════════════════════════════════════════════
#  hg38 autosome sizes
# ══════════════════════════════════════════════════════════════════════════════
AUTOSOME_SIZE_HG38 = {
    "chr1":  248_956_422, "chr2":  242_193_529, "chr3":  198_295_559,
    "chr4":  190_214_555, "chr5":  181_538_259, "chr6":  170_805_979,
    "chr7":  159_345_973, "chr8":  145_138_636, "chr9":  138_394_717,
    "chr10": 133_797_422, "chr11": 135_086_622, "chr12": 133_275_309,
    "chr13": 114_364_328, "chr14": 107_043_718, "chr15": 101_991_189,
    "chr16":  90_338_345, "chr17":  83_257_441, "chr18":  80_373_285,
    "chr19":  58_617_616, "chr20":  64_444_167, "chr21":  46_709_983,
    "chr22":  50_818_468,
}
for i in range(1, 23):
    AUTOSOME_SIZE_HG38[str(i)] = AUTOSOME_SIZE_HG38[f"chr{i}"]

TOTAL_AUTOSOME_HG38 = sum(AUTOSOME_SIZE_HG38[f"chr{i}"] for i in range(1, 23))


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1 — VCF PARSING
# ══════════════════════════════════════════════════════════════════════════════
def parse_vcf_for_homozygosity(vcf_path, chromosomes,
                                min_gq=20, min_dp=8,
                                skip_multiallelic=True,
                                skip_indels=True):
    per_chrom = {}
    all_sites = []

    for chrom in chromosomes:
        per_chrom[chrom] = dict(
            hom_ref=0, hom_alt=0, het=0, nocall=0,
            positions=[], zygosity=[]
        )

    vcf = VCF(vcf_path)
    chrom_set = set(chromosomes)

    for variant in vcf:
        chrom = variant.CHROM
        if chrom not in chrom_set:
            continue
        if skip_multiallelic and len(variant.ALT) != 1:
            continue
        if skip_indels:
            if (len(variant.REF) != 1) or (len(variant.ALT[0]) != 1):
                continue

        gt = variant.genotypes[0]
        a1, a2 = gt[0], gt[1]

        gq_arr = variant.format("GQ")
        dp_arr = variant.format("DP")
        gq = int(gq_arr[0][0]) if gq_arr is not None else 0
        dp = int(dp_arr[0][0]) if dp_arr is not None else 0

        if a1 < 0 or a2 < 0 or gq < min_gq or dp < min_dp:
            zyg = 0
            per_chrom[chrom]["nocall"] += 1
        elif a1 == 0 and a2 == 0:
            zyg = 1
            per_chrom[chrom]["hom_ref"] += 1
        elif a1 == a2:
            zyg = 2
            per_chrom[chrom]["hom_alt"] += 1
        else:
            zyg = 3
            per_chrom[chrom]["het"] += 1

        per_chrom[chrom]["positions"].append(variant.POS)
        per_chrom[chrom]["zygosity"].append(zyg)
        all_sites.append((chrom, variant.POS, zyg))

    vcf.close()

    global_counts = dict(hom_ref=0, hom_alt=0, het=0, nocall=0)
    for chrom, counts in per_chrom.items():
        for key in global_counts:
            global_counts[key] += counts[key]

    return {"per_chrom": per_chrom, "global": global_counts, "all_sites": all_sites}


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2 — STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
def compute_homozygosity_stats(global_counts):
    hom_ref = global_counts["hom_ref"]
    hom_alt = global_counts["hom_alt"]
    het     = global_counts["het"]
    nocall  = global_counts["nocall"]

    total_callable = hom_ref + hom_alt + het
    total_sites    = total_callable + nocall

    if total_callable == 0:
        return {k: float("nan") for k in
                ["hom_pct", "het_pct", "hom_ref_pct", "hom_alt_pct",
                 "f_stat", "het_hom_ratio", "total_callable", "total_sites",
                 "nocall_pct", "hom_ref", "hom_alt", "het", "nocall"]}

    hom_pct      = 100.0 * (hom_ref + hom_alt) / total_callable
    het_pct      = 100.0 * het / total_callable
    hom_ref_pct  = 100.0 * hom_ref  / total_callable
    hom_alt_pct  = 100.0 * hom_alt  / total_callable
    nocall_pct   = 100.0 * nocall / total_sites if total_sites > 0 else 0.0
    het_hom_ratio = het / (hom_ref + hom_alt) if (hom_ref + hom_alt) > 0 else float("nan")

    total_alleles    = 2 * total_callable
    alt_allele_count = (1 * het) + (2 * hom_alt)
    p = alt_allele_count / total_alleles if total_alleles > 0 else 0.0
    q = 1.0 - p

    expected_het = 2 * p * q
    observed_het = het / total_callable
    f_stat = (expected_het - observed_het) / expected_het if expected_het > 0 else float("nan")

    return {
        "hom_ref": hom_ref, "hom_alt": hom_alt,
        "het": het,         "nocall": nocall,
        "total_callable": total_callable,
        "total_sites":    total_sites,
        "hom_pct":        round(hom_pct, 4),
        "het_pct":        round(het_pct, 4),
        "hom_ref_pct":    round(hom_ref_pct, 4),
        "hom_alt_pct":    round(hom_alt_pct, 4),
        "nocall_pct":     round(nocall_pct, 4),
        "het_hom_ratio":  round(het_hom_ratio, 4),
        "alt_af":         round(p, 6),
        "exp_het":        round(expected_het, 6),
        "obs_het":        round(observed_het, 6),
        "f_stat":         round(f_stat, 6),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3 — ROH DETECTION
# ══════════════════════════════════════════════════════════════════════════════
def detect_roh(per_chrom, window_snps=50, min_snps=50,
               min_hom_frac=0.95, min_length_bp=500_000):
    roh_list = []

    for chrom, data in per_chrom.items():
        positions = np.array(data["positions"])
        zygosity  = np.array(data["zygosity"])

        callable_mask = zygosity > 0
        pos_c = positions[callable_mask]
        zyg_c = zygosity[callable_mask]

        n = len(pos_c)
        if n < window_snps:
            continue

        is_hom = (zyg_c == 1) | (zyg_c == 2)

        cum_hom           = np.concatenate([[0], np.cumsum(is_hom)])
        window_hom_counts = cum_hom[window_snps:] - cum_hom[:-window_snps]
        window_hom_frac   = window_hom_counts / window_snps
        in_roh_window     = window_hom_frac >= min_hom_frac

        snp_in_roh = np.zeros(n, dtype=bool)
        for wi in range(len(in_roh_window)):
            if in_roh_window[wi]:
                snp_in_roh[wi: wi + window_snps] = True

        padded      = np.concatenate([[False], snp_in_roh, [False]])
        transitions = np.diff(padded.astype(int))
        starts_idx  = np.where(transitions ==  1)[0]
        ends_idx    = np.where(transitions == -1)[0]

        for s_idx, e_idx in zip(starts_idx, ends_idx):
            n_snps    = e_idx - s_idx
            start_bp  = int(pos_c[s_idx])
            end_bp    = int(pos_c[e_idx - 1])
            length_bp = end_bp - start_bp
            hom_frac  = float(is_hom[s_idx:e_idx].mean())

            if n_snps >= min_snps and length_bp >= min_length_bp:
                roh_list.append({
                    "chrom": chrom, "start": start_bp, "end": end_bp,
                    "length_bp": length_bp, "n_snps": n_snps,
                    "hom_frac": round(hom_frac, 4),
                })

    return roh_list


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4 — FROH
# ══════════════════════════════════════════════════════════════════════════════
def compute_froh(roh_list):
    if not roh_list:
        return {"froh": 0.0, "total_roh_bp": 0, "n_roh": 0,
                "mean_roh_length_mb": 0.0, "largest_roh_mb": 0.0}

    total_roh_bp = sum(r["length_bp"] for r in roh_list)
    lengths_mb   = [r["length_bp"] / 1e6 for r in roh_list]

    return {
        "froh":               round(total_roh_bp / TOTAL_AUTOSOME_HG38, 6),
        "total_roh_bp":       total_roh_bp,
        "n_roh":              len(roh_list),
        "mean_roh_length_mb": round(float(np.mean(lengths_mb)), 3),
        "largest_roh_mb":     round(float(np.max(lengths_mb)), 3),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 5 — PUBLICATION-READY FIGURE
#
#  Layout (white background, Nature/Cell journal style):
#
#   ┌──────────────────────┬──────────────────────────────────────┐
#   │  Panel A             │  Panel B                             │
#   │  Zygosity stacked    │  Homozygosity % lollipop + F/FROH   │
#   │  bar (proportional)  │  dot plot                            │
#   ├──────────────────────┴──────────────────────────────────────┤
#   │  Panel C  Per-chromosome homozygosity heatmap (full width)  │
#   ├─────────────────────────────────────────────────────────────┤
#   │  Panel D  ROH genome ideogram (full width)                  │
#   └─────────────────────────────────────────────────────────────┘
# ══════════════════════════════════════════════════════════════════════════════

def _panel_label(ax, label, x=-0.12, y=1.06):
    """Add bold panel letter (A, B, C…) — standard Nature figure style."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="top", ha="left",
            color="#111111")


def plot_homozygosity_report(results, roh_by_sample, output_path, dpi=300):
    """
    Four-panel publication figure.

    Panel A — Proportional stacked bar: zygosity class composition per sample
    Panel B — Lollipop chart: homozygosity % with F-statistic annotation
    Panel C — Heatmap: per-chromosome homozygosity % (samples × chromosomes)
    Panel D — Genome ideogram: ROH segments along autosomes
    """
    sample_ids    = list(results.keys())
    n_samples     = len(sample_ids)
    sample_colors = PUB_COLORS["samples"][:n_samples]
    chromosomes   = [f"chr{i}" for i in range(1, 23)]

    # ── Figure geometry ───────────────────────────────────────────────────────
    fig = plt.figure(figsize=(18, 20), facecolor="white")
    gs  = GridSpec(
        4, 2, figure=fig,
        height_ratios=[1.1, 1.1, 1.0, 1.2],
        hspace=0.55, wspace=0.38,
        left=0.08, right=0.96, top=0.94, bottom=0.04,
    )

    ax_a = fig.add_subplot(gs[0, 0])   # Zygosity stacked bar
    ax_b = fig.add_subplot(gs[0, 1])   # Lollipop hom%
    ax_c = fig.add_subplot(gs[1, :])   # Per-chrom heatmap
    ax_d = fig.add_subplot(gs[2, :])   # ROH genome ideogram
    ax_e = fig.add_subplot(gs[3, :])   # Summary stats table

    # ══════════════════════════════════════════════════════════════════════
    #  PANEL A — Proportional stacked bar (zygosity composition)
    # ══════════════════════════════════════════════════════════════════════
    zyg_order  = ["hom_ref", "hom_alt", "het", "nocall"]
    zyg_labels = ["Hom REF (0/0)", "Hom ALT (1/1)", "Heterozygous (0/1)", "No-call"]
    zyg_colors = [PUB_COLORS[k] for k in zyg_order]

    x_pos  = np.arange(n_samples)
    bar_w  = 0.5
    bottoms = np.zeros(n_samples)

    for zk, zl, zc in zip(zyg_order, zyg_labels, zyg_colors):
        # Proportional: each segment = % of total_sites
        vals = np.array([
            100.0 * results[s]["global"][zk] /
            max(results[s]["stats"]["total_sites"], 1)
            for s in sample_ids
        ])
        ax_a.bar(x_pos, vals, width=bar_w, bottom=bottoms,
                 color=zc, label=zl, edgecolor="white", linewidth=0.6)
        bottoms += vals

    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(sample_ids, rotation=25, ha="right", fontsize=9)
    ax_a.set_ylabel("Proportion of all sites (%)", fontsize=9)
    ax_a.set_ylim(0, 105)
    ax_a.set_title("Zygosity Class Composition", fontsize=11, fontweight="bold", pad=8)
    ax_a.yaxis.set_major_formatter(mticker.FormatStrFormatter("%g%%"))
    ax_a.grid(axis="y", alpha=0.5)
    ax_a.legend(fontsize=7.5, loc="lower right",
                ncol=1, handlelength=1.2, handletextpad=0.5)
    _panel_label(ax_a, "A")

    # ══════════════════════════════════════════════════════════════════════
    #  PANEL B — Lollipop chart: Homozygosity % per sample
    #            with F-statistic and FROH annotated
    # ══════════════════════════════════════════════════════════════════════
    hom_pcts = [results[s]["stats"]["hom_pct"] for s in sample_ids]
    het_pcts = [results[s]["stats"]["het_pct"] for s in sample_ids]

    y_pos = np.arange(n_samples)

    # Hom% lollipop (left side)
    ax_b.hlines(y_pos, 0, hom_pcts, color="#cccccc", linewidth=1.5, zorder=1)
    ax_b.scatter(hom_pcts, y_pos, color=PUB_COLORS["hom_alt"],
                 s=80, zorder=3, label="Homozygosity %", edgecolors="white", linewidths=0.8)

    # Het% markers (hollow)
    ax_b.scatter(het_pcts, y_pos, color=PUB_COLORS["het"],
                 s=80, zorder=3, label="Heterozygosity %",
                 marker="D", edgecolors=PUB_COLORS["het"], linewidths=1.5,
                 facecolors="white")

    # Annotate F-statistic and FROH
    for yi, sid in enumerate(sample_ids):
        f    = results[sid]["stats"]["f_stat"]
        froh = roh_by_sample[sid]["froh_stats"]["froh"]
        n_roh = roh_by_sample[sid]["froh_stats"]["n_roh"]
        ax_b.text(
            101, yi,
            f"F={f:+.4f}   F_ROH={froh:.4f}   n_ROH={n_roh}",
            va="center", ha="left", fontsize=7.5, color="#555555",
            fontfamily="monospace",
        )

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(sample_ids, fontsize=9)
    ax_b.set_xlabel("Percentage of callable SNPs (%)", fontsize=9)
    ax_b.set_xlim(0, 100)
    ax_b.set_title("Homozygosity & Heterozygosity\nwith F-statistic & F_ROH",
                   fontsize=11, fontweight="bold", pad=8)
    ax_b.legend(fontsize=8, loc="lower right")
    ax_b.grid(axis="x", alpha=0.5)
    ax_b.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g%%"))

    # Reference line at 50% (HWE expectation for diploid)
    ax_b.axvline(50, color="#bbbbbb", lw=0.8, linestyle="--", alpha=0.8,
                 label="_nolegend_")
    ax_b.text(50.5, -0.55, "50%\n(HWE)", fontsize=6.5, color="#aaaaaa", va="top")

    _panel_label(ax_b, "B")

    # ══════════════════════════════════════════════════════════════════════
    #  PANEL C — Per-chromosome homozygosity heatmap
    #            Rows = samples, Columns = chromosomes
    # ══════════════════════════════════════════════════════════════════════
    hom_matrix = np.full((n_samples, len(chromosomes)), np.nan)

    for si, sid in enumerate(sample_ids):
        for ci, chrom in enumerate(chromosomes):
            cd = results[sid]["per_chrom"].get(chrom, {})
            hr  = cd.get("hom_ref", 0)
            ha  = cd.get("hom_alt", 0)
            ht  = cd.get("het", 0)
            tot = hr + ha + ht
            if tot > 0:
                hom_matrix[si, ci] = 100.0 * (hr + ha) / tot

    # Custom white→orange→dark-red colormap
    cmap_hom = LinearSegmentedColormap.from_list(
        "hom_heat", PUB_COLORS["heatmap_cmap"], N=256
    )

    im = ax_c.imshow(
        hom_matrix, cmap=cmap_hom,
        vmin=0, vmax=100,
        aspect="auto", interpolation="nearest",
    )

    # X-axis: chromosome numbers
    ax_c.set_xticks(range(len(chromosomes)))
    ax_c.set_xticklabels(
        [c.replace("chr", "") for c in chromosomes],
        fontsize=8
    )
    ax_c.set_xlabel("Chromosome", fontsize=9)

    # Y-axis: sample labels
    ax_c.set_yticks(range(n_samples))
    ax_c.set_yticklabels(sample_ids, fontsize=9)

    # Cell value annotations
    for si in range(n_samples):
        for ci in range(len(chromosomes)):
            val = hom_matrix[si, ci]
            if not np.isnan(val):
                # White text on dark cells, dark text on light cells
                txt_color = "white" if val > 70 else "#333333"
                ax_c.text(ci, si, f"{val:.0f}",
                          ha="center", va="center",
                          fontsize=6, color=txt_color, fontweight="bold")

    # Colorbar
    cbar = fig.colorbar(im, ax=ax_c, orientation="vertical",
                        fraction=0.012, pad=0.01, shrink=0.85)
    cbar.set_label("Homozygosity (%)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax_c.set_title("Per-Chromosome Homozygosity (%)",
                   fontsize=11, fontweight="bold", pad=8)
    ax_c.tick_params(top=False, bottom=True)
    _panel_label(ax_c, "C")

    # ══════════════════════════════════════════════════════════════════════
    #  PANEL D — ROH Genome Ideogram
    #            Horizontal bars per sample; ROH segments in amber
    #            Chromosomes as alternating light/white bands
    # ══════════════════════════════════════════════════════════════════════
    chrom_sizes = {c: AUTOSOME_SIZE_HG38.get(c, 200_000_000) for c in chromosomes}

    # Build cumulative offsets (Mb)
    cum_offset = {}
    offset_mb  = 0.0
    for chrom in chromosomes:
        cum_offset[chrom] = offset_mb
        offset_mb += chrom_sizes[chrom] / 1e6
    total_genome_mb = offset_mb

    # Alternating chromosome bands (very subtle grey/white)
    for ci, chrom in enumerate(chromosomes):
        band_start = cum_offset[chrom]
        band_end   = band_start + chrom_sizes[chrom] / 1e6
        band_color = "#F7F7F7" if ci % 2 == 0 else "#FFFFFF"
        ax_d.axvspan(band_start, band_end,
                     facecolor=band_color, alpha=1.0, zorder=0,
                     linewidth=0)
        # Chromosome label
        mid = (band_start + band_end) / 2
        ax_d.text(mid, n_samples + 0.05,
                  chrom.replace("chr", ""),
                  ha="center", va="bottom", fontsize=6.5, color="#666666")

    # Chromosome separator lines
    for chrom in chromosomes:
        ax_d.axvline(cum_offset[chrom], color="#dddddd", lw=0.5, zorder=1)

    # Track height and spacing
    track_h = 0.55
    for si, sid in enumerate(sample_ids):
        y_ctr = si

        # Grey baseline track
        ax_d.hlines(y_ctr, 0, total_genome_mb,
                    colors="#dddddd", lw=0.8, zorder=1)

        # ROH segments
        for roh in roh_by_sample[sid]["roh_list"]:
            x0 = cum_offset[roh["chrom"]] + roh["start"] / 1e6
            w  = roh["length_bp"] / 1e6
            rect = mpatches.FancyArrow(
                x0, y_ctr, w, 0,
                width=track_h,
                head_width=track_h,
                head_length=0,
                length_includes_head=True,
                facecolor=PUB_COLORS["roh"],
                edgecolor="none",
                zorder=2, alpha=0.85,
            )
            # Simpler: use barh
            ax_d.barh(
                y_ctr, w, left=x0, height=track_h,
                color=PUB_COLORS["roh"], edgecolor="none",
                alpha=0.85, zorder=2,
            )

    ax_d.set_yticks(range(n_samples))
    ax_d.set_yticklabels(sample_ids, fontsize=9)
    ax_d.set_xlim(0, total_genome_mb)
    ax_d.set_ylim(-0.6, n_samples)
    ax_d.set_xlabel("Genomic position (Mb, GRCh38 autosomes)", fontsize=9)
    ax_d.grid(axis="x", alpha=0.3, linewidth=0.5)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)
    ax_d.spines["left"].set_visible(False)
    ax_d.tick_params(left=False)

    roh_patch = mpatches.Patch(facecolor=PUB_COLORS["roh"], alpha=0.85,
                                label=f"ROH segment (≥{500_000//1000} kb)")
    ax_d.legend(handles=[roh_patch], fontsize=8, loc="upper right")
    ax_d.set_title("Runs of Homozygosity (ROH) — Autosome Ideogram",
                   fontsize=11, fontweight="bold", pad=8)
    _panel_label(ax_d, "D")

    # ══════════════════════════════════════════════════════════════════════
    #  PANEL E — Summary statistics table
    # ══════════════════════════════════════════════════════════════════════
    ax_e.axis("off")

    # Build table data
    col_labels = [
        "Sample", "Callable SNPs",
        "Hom %", "Het %", "No-call %",
        "Het/Hom", "F-stat", "n_ROH", "Total ROH (Mb)", "F_ROH",
        "Interpretation",
    ]

    table_data = []
    for sid in sample_ids:
        st   = results[sid]["stats"]
        froh = roh_by_sample[sid]["froh_stats"]

        froh_v = froh["froh"]
        if froh_v >= 0.25:
            interp = "High (≥1st cousin)"
        elif froh_v >= 0.10:
            interp = "Elevated (2nd cousin)"
        elif froh_v >= 0.05:
            interp = "Mild (distant)"
        elif froh_v >= 0.01:
            interp = "Low (background)"
        else:
            interp = "Normal (outbred)"

        table_data.append([
            sid,
            f"{st['total_callable']:,}",
            f"{st['hom_pct']:.2f}%",
            f"{st['het_pct']:.2f}%",
            f"{st['nocall_pct']:.2f}%",
            f"{st['het_hom_ratio']:.3f}",
            f"{st['f_stat']:+.4f}",
            str(froh["n_roh"]),
            f"{froh['total_roh_bp']/1e6:.2f}",
            f"{froh_v:.4f}",
            interp,
        ])

    tbl = ax_e.table(
        cellText    = table_data,
        colLabels   = col_labels,
        loc         = "center",
        cellLoc     = "center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1.0, 1.8)

    # Style header row
    for ci in range(len(col_labels)):
        cell = tbl[(0, ci)]
        cell.set_facecolor("#2C5F8A")
        cell.set_text_props(color="white", fontweight="bold", fontsize=8.5)
        cell.set_edgecolor("white")

    # Style data rows — alternating
    row_bg = ["#F0F4FA", "#FFFFFF"]
    for ri in range(1, len(table_data) + 1):
        for ci in range(len(col_labels)):
            cell = tbl[(ri, ci)]
            cell.set_facecolor(row_bg[(ri - 1) % 2])
            cell.set_edgecolor("#DDDDDD")

    ax_e.set_title("Summary Statistics Table",
                   fontsize=11, fontweight="bold", pad=12)
    _panel_label(ax_e, "E", x=-0.01)

    # ══════════════════════════════════════════════════════════════════════
    #  Figure title and save
    # ══════════════════════════════════════════════════════════════════════
    fig.suptitle(
        "Genome-Wide Homozygosity and Runs of Homozygosity (ROH) Analysis\n"
        f"n={n_samples} sample{'s' if n_samples > 1 else ''}  |  "
        "GRCh38 autosomes (chr1–22)  |  SNPs only  |  GQ≥20, DP≥8",
        fontsize=13, fontweight="bold", y=0.97, color="#111111",
    )

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  ✓ Publication figure saved: {output_path}")


# ══════════════════════════════════════════════════════════════════════════════
#  STEP 6 — TSV / TEXT OUTPUT
# ══════════════════════════════════════════════════════════════════════════════
def write_summary_tsv(results, roh_by_sample, output_path):
    rows = []
    for sid, data in results.items():
        st   = data["stats"]
        froh = roh_by_sample[sid]["froh_stats"]
        froh_v = froh["froh"]

        if froh_v >= 0.25:
            interp = "HIGH FROH (≥0.25): consistent with 1st-cousin or closer"
        elif froh_v >= 0.10:
            interp = "ELEVATED FROH (0.10-0.25): consistent with 2nd-cousin"
        elif froh_v >= 0.05:
            interp = "MILDLY ELEVATED FROH (0.05-0.10): possible distant relatedness"
        elif froh_v >= 0.01:
            interp = "LOW FROH (0.01-0.05): background relatedness"
        else:
            interp = "NORMAL FROH (<0.01): consistent with outbred individual"

        rows.append({
            "sample_id":      sid,
            "total_callable": st["total_callable"],
            "hom_ref":        st["hom_ref"],
            "hom_alt":        st["hom_alt"],
            "het":            st["het"],
            "nocall":         st["nocall"],
            "hom_pct":        st["hom_pct"],
            "het_pct":        st["het_pct"],
            "hom_ref_pct":    st["hom_ref_pct"],
            "hom_alt_pct":    st["hom_alt_pct"],
            "nocall_pct":     st["nocall_pct"],
            "het_hom_ratio":  st["het_hom_ratio"],
            "alt_af":         st["alt_af"],
            "obs_het":        st["obs_het"],
            "exp_het":        st["exp_het"],
            "f_stat":         st["f_stat"],
            "n_roh":          froh["n_roh"],
            "total_roh_mb":   round(froh["total_roh_bp"] / 1e6, 3),
            "mean_roh_mb":    froh["mean_roh_length_mb"],
            "largest_roh_mb": froh["largest_roh_mb"],
            "froh":           froh["froh"],
            "interpretation": interp,
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"  ✓ Summary TSV: {output_path}")
    return df


def write_roh_tsv(roh_by_sample, output_path):
    rows = []
    for sid, data in roh_by_sample.items():
        for roh in data["roh_list"]:
            rows.append({
                "sample_id": sid,
                "chrom":     roh["chrom"],
                "start":     roh["start"],
                "end":       roh["end"],
                "length_bp": roh["length_bp"],
                "length_mb": round(roh["length_bp"] / 1e6, 3),
                "n_snps":    roh["n_snps"],
                "hom_frac":  roh["hom_frac"],
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["sample_id", "length_bp"], ascending=[True, False])
    df.to_csv(output_path, sep="\t", index=False)
    print(f"  ✓ ROH segments TSV: {output_path}")


def write_text_report(results, roh_by_sample, output_path):
    with open(output_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("  GENOME-WIDE HOMOZYGOSITY ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        for sid, data in results.items():
            st   = data["stats"]
            froh = roh_by_sample[sid]["froh_stats"]

            f.write(f"Sample: {sid}\n")
            f.write("-" * 60 + "\n")
            f.write(f"  Total callable SNPs  : {st['total_callable']:>12,}\n")
            f.write(f"  Homozygous REF (0/0) : {st['hom_ref']:>12,}  ({st['hom_ref_pct']:6.2f}%)\n")
            f.write(f"  Homozygous ALT (1/1) : {st['hom_alt']:>12,}  ({st['hom_alt_pct']:6.2f}%)\n")
            f.write(f"  Heterozygous  (0/1)  : {st['het']:>12,}  ({st['het_pct']:6.2f}%)\n")
            f.write(f"  No-call / filtered   : {st['nocall']:>12,}  ({st['nocall_pct']:6.2f}%)\n\n")
            f.write(f"  ── Homozygosity Metrics ──\n")
            f.write(f"  Homozygosity %        : {st['hom_pct']:8.4f}%\n")
            f.write(f"  Heterozygosity %      : {st['het_pct']:8.4f}%\n")
            f.write(f"  Het / Hom ratio       : {st['het_hom_ratio']:8.4f}\n")
            f.write(f"  Observed Het rate     : {st['obs_het']:8.6f}\n")
            f.write(f"  Expected Het (HWE)    : {st['exp_het']:8.6f}\n")
            f.write(f"  F-statistic (AF-based): {st['f_stat']:+8.6f}\n\n")
            f.write(f"  ── Runs of Homozygosity (ROH) ──\n")
            f.write(f"  N ROH segments (≥500kb): {froh['n_roh']:>8}\n")
            f.write(f"  Total ROH length        : {froh['total_roh_bp']/1e6:>8.2f} Mb\n")
            f.write(f"  Mean ROH length         : {froh['mean_roh_length_mb']:>8.3f} Mb\n")
            f.write(f"  Largest ROH             : {froh['largest_roh_mb']:>8.3f} Mb\n")
            f.write(f"  FROH (genome fraction)  : {froh['froh']:>8.6f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("Output files:\n")
        f.write("  homozygosity_summary.tsv   — per-sample statistics\n")
        f.write("  roh_segments.tsv           — all ROH segments\n")
        f.write("  homozygosity_report.png    — 5-panel figure\n")
        f.write("  analysis_report.txt        — this file\n")
        f.write("=" * 80 + "\n")

    print(f"  ✓ Text report: {output_path}")


# ══════════════════════════════════════════════════════════════════════════════
#  CONFIG LOADER
# ══════════════════════════════════════════════════════════════════════════════
def load_config(config_file):
    cfg = DEFAULTS.copy()

    if not os.path.exists(config_file):
        print(f"  [INFO] No config.ini at '{config_file}' — using built-in defaults.")
        return cfg

    ini = configparser.ConfigParser()
    ini.read(config_file)

    if "PATHS"   in ini: cfg["results_dir"] = ini["PATHS"].get("results_dir", cfg["results_dir"])
    if "SAMPLES" in ini: cfg["samples"]     = {k: v for k, v in ini["SAMPLES"].items()}
    if "FILTERS" in ini:
        cfg["min_gq"]            = int(ini["FILTERS"].get("min_gq",   cfg["min_gq"]))
        cfg["min_dp"]            = int(ini["FILTERS"].get("min_dp",   cfg["min_dp"]))
        cfg["skip_multiallelic"] = ini["FILTERS"].getboolean("skip_multiallelic", cfg["skip_multiallelic"])
        cfg["skip_indels"]       = ini["FILTERS"].getboolean("skip_indels",       cfg["skip_indels"])
    if "ROH" in ini:
        cfg["roh_window_snps"]   = int(ini["ROH"].get("roh_window_snps",   cfg["roh_window_snps"]))
        cfg["roh_min_snps"]      = int(ini["ROH"].get("roh_min_snps",      cfg["roh_min_snps"]))
        cfg["roh_min_hom_frac"]  = float(ini["ROH"].get("roh_min_hom_frac",  cfg["roh_min_hom_frac"]))
        cfg["roh_min_length_bp"] = int(ini["ROH"].get("roh_min_length_bp", cfg["roh_min_length_bp"]))
    if "VISUALIZATION" in ini:
        cfg["output_dpi"] = int(ini["VISUALIZATION"].get("output_dpi", cfg["output_dpi"]))

    return cfg


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="Genome-Wide Homozygosity Estimator — publication-ready output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("config", nargs="?", default="config.ini",
                        help="Path to config.ini (default: config.ini)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    if cfg.get("numeric_chroms"):
        cfg["chromosomes"] = [str(i) for i in range(1, 23)]

    results_dir = Path(cfg["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("  GENOME-WIDE HOMOZYGOSITY ESTIMATOR  (publication mode)")
    print("=" * 80)
    print(f"\n  Results directory  : {results_dir}")
    print(f"  Samples            : {list(cfg['samples'].keys())}")
    print(f"  Chromosomes        : chr1–22 (autosomes)")
    print(f"  Quality filters    : GQ≥{cfg['min_gq']}, DP≥{cfg['min_dp']}")
    print(f"  ROH parameters     : window={cfg['roh_window_snps']} SNPs, "
          f"min_frac={cfg['roh_min_hom_frac']}, "
          f"min_len={cfg['roh_min_length_bp']//1000} kb\n")

    # ── Check VCFs ─────────────────────────────────────────────────────────
    print("[STEP 1] Checking VCF files...")
    missing = []
    for sid, vcf_path in cfg["samples"].items():
        if not os.path.exists(vcf_path):
            missing.append(vcf_path)
        elif not os.path.exists(vcf_path + ".tbi"):
            print(f"  [WARNING] No tabix index: tabix -p vcf {vcf_path}")
        else:
            print(f"  [OK] {sid}  ({os.path.getsize(vcf_path)/1e6:.0f} MB)")

    if missing:
        print("\n[ERROR] Missing VCF files:")
        for m in missing: print(f"  {m}")
        sys.exit(1)

    # ── Parse VCFs ─────────────────────────────────────────────────────────
    print("\n[STEP 2] Parsing VCFs...")
    all_results = {}

    for sid, vcf_path in cfg["samples"].items():
        print(f"  → {sid}")
        raw   = parse_vcf_for_homozygosity(
            vcf_path, cfg["chromosomes"],
            cfg["min_gq"], cfg["min_dp"],
            cfg["skip_multiallelic"], cfg["skip_indels"],
        )
        stats = compute_homozygosity_stats(raw["global"])
        all_results[sid] = {
            "per_chrom": raw["per_chrom"],
            "global":    raw["global"],
            "stats":     stats,
        }
        print(f"    Callable: {stats['total_callable']:,} SNPs  |  "
              f"Hom: {stats['hom_pct']:.2f}%  |  "
              f"Het: {stats['het_pct']:.2f}%  |  "
              f"F: {stats['f_stat']:+.4f}")

    # ── ROH ────────────────────────────────────────────────────────────────
    print("\n[STEP 3] Detecting ROH...")
    roh_by_sample = {}
    for sid in all_results:
        roh_list   = detect_roh(all_results[sid]["per_chrom"],
                                cfg["roh_window_snps"], cfg["roh_min_snps"],
                                cfg["roh_min_hom_frac"], cfg["roh_min_length_bp"])
        froh_stats = compute_froh(roh_list)
        roh_by_sample[sid] = {"roh_list": roh_list, "froh_stats": froh_stats}
        print(f"  {sid}: {froh_stats['n_roh']} ROH  |  "
              f"{froh_stats['total_roh_bp']/1e6:.2f} Mb  |  "
              f"FROH={froh_stats['froh']:.4f}")

    # ── Write outputs ───────────────────────────────────────────────────────
    print("\n[STEP 4] Writing outputs...")
    write_summary_tsv(all_results, roh_by_sample,
                      str(results_dir / "homozygosity_summary.tsv"))
    write_roh_tsv(roh_by_sample,
                  str(results_dir / "roh_segments.tsv"))
    write_text_report(all_results, roh_by_sample,
                      str(results_dir / "analysis_report.txt"))

    # ── Figure ─────────────────────────────────────────────────────────────
    print("\n[STEP 5] Generating publication figure...")
    plot_homozygosity_report(
        results       = all_results,
        roh_by_sample = roh_by_sample,
        output_path   = str(results_dir / "homozygosity_report.png"),
        dpi           = cfg["output_dpi"],
    )

    print("\n" + "=" * 80)
    print("  COMPLETE")
    print(f"  Output: {results_dir}/")
    print("    ├── homozygosity_summary.tsv")
    print("    ├── roh_segments.tsv")
    print("    ├── analysis_report.txt")
    print("    └── homozygosity_report.png   ← 5-panel publication figure")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
