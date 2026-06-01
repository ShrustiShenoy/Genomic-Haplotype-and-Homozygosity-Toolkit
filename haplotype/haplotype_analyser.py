#!/usr/bin/env python3
"""
================================================================================
  EXTENDED HAPLOTYPE ANALYZER
  
  Supports:
    ✓ Single VOI (Variant of Interest) analysis
    ✓ Extended haplotype phasing and comparison
    ✓ WGS (whole genome) and LCG (low-coverage genome) samples
    ✓ IBS (Identity-by-State) similarity calculation
    ✓ Configurable via config.ini
  
  Dependencies:
    pip install cyvcf2 matplotlib pandas numpy configparser
  
  Usage:
    python3 haplotype_analyzer.py config.ini
    
  Output:
    <results_dir>/
    ├── haplotype_grid.png           [Panel A]
    ├── ibs_heatmap.png              [Panel B]
    ├── alt_allele_frequency.png     [Panel C]
    ├── voi_summary.png              [Panel D]
    ├── haplotype_results.tsv        [Full genotype table]
    ├── voi_calls.tsv                [VOI-only calls]
    └── analysis_summary.txt         [Text report]

================================================================================
"""

# ─── Standard Library ────────────────────────────────────────────────────────
import os
import sys
import warnings
import configparser
from pathlib import Path

warnings.filterwarnings("ignore")

# ─── Third-party Dependencies ───────────────────────────────────────────────
try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from cyvcf2 import VCF
except ImportError as e:
    sys.exit(
        f"\n[ERROR] Missing required package: {e}\n"
        "Install dependencies with:\n"
        "  pip install cyvcf2 matplotlib pandas numpy configparser\n"
    )

# ═══════════════════════════════════════════════════════════════════════════════
#  CONFIGURATION LOADER
# ═══════════════════════════════════════════════════════════════════════════════

class ConfigManager:
    """Load and validate configuration from INI file"""
    
    def __init__(self, config_file):
        if not os.path.exists(config_file):
            raise FileNotFoundError(
                f"Config file not found: {config_file}\n"
                "Please create config.ini with all required settings.\n"
            )
        
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        self._validate()
    
    def _validate(self):
        """Validate required sections and keys"""
        required_sections = {
            'PATHS': ['vcf_dir', 'results_dir'],
            'VCF_FILES': [],  # At least one sample
            'VARIANT_OF_INTEREST': [
                'chrom_prefix', 'chrom', 'position', 'ref_allele',
                'alt_allele', 'rsid', 'gene', 'cdna', 'protein'
            ],
            'ANALYSIS_PARAMS': ['flank_bp', 'include_multiallelic', 'include_indels'],
            'VISUALIZATION': ['output_dpi', 'color_scheme', 'separate_panels']
        }
        
        for section, keys in required_sections.items():
            if section not in self.config:
                raise ValueError(f"Missing section: [{section}]")
            
            if section != 'VCF_FILES':  # VCF_FILES items are dynamic
                for key in keys:
                    if key not in self.config[section]:
                        raise ValueError(f"Missing key: {key} in [{section}]")
        
        # At least one sample
        if len(self.config['VCF_FILES']) == 0:
            raise ValueError("No samples defined in [VCF_FILES] section")
    
    def get_samples(self):
        """Get list of (LIMS_ID, vcf_stem) tuples from config"""
        return [
            (lims_id, stem)
            for lims_id, stem in self.config['VCF_FILES'].items()
        ]
    
    def get_vcf_dir(self):
        return self.config['PATHS']['vcf_dir']
    
    def get_results_dir(self):
        return self.config['PATHS']['results_dir']
    
    def get_voi_params(self):
        """Return VOI parameters as dict"""
        chrom_prefix = self.config['VARIANT_OF_INTEREST'].get('chrom_prefix', '').strip()
        chrom = self.config['VARIANT_OF_INTEREST']['chrom'].strip()
        
        return {
            'chrom_prefix': chrom_prefix,
            'chrom': chrom,
            'pos': int(self.config['VARIANT_OF_INTEREST']['position']),
            'ref': self.config['VARIANT_OF_INTEREST']['ref_allele'].strip(),
            'alt': self.config['VARIANT_OF_INTEREST']['alt_allele'].strip(),
            'rsid': self.config['VARIANT_OF_INTEREST']['rsid'].strip(),
            'gene': self.config['VARIANT_OF_INTEREST']['gene'].strip(),
            'cdna': self.config['VARIANT_OF_INTEREST']['cdna'].strip(),
            'protein': self.config['VARIANT_OF_INTEREST']['protein'].strip(),
        }
    
    def get_flank(self):
        return int(self.config['ANALYSIS_PARAMS']['flank_bp'])
    
    def get_include_multiallelic(self):
        return self.config['ANALYSIS_PARAMS'].getboolean('include_multiallelic')
    
    def get_include_indels(self):
        return self.config['ANALYSIS_PARAMS'].getboolean('include_indels')
    
    def get_output_dpi(self):
        return int(self.config['VISUALIZATION']['output_dpi'])
    
    def get_color_scheme(self):
        scheme = self.config['VISUALIZATION'].get('color_scheme', 'dark').lower()
        return scheme if scheme in ['dark', 'light'] else 'dark'
    
    def get_separate_panels(self):
        return self.config['VISUALIZATION'].getboolean('separate_panels')

# ═══════════════════════════════════════════════════════════════════════════════
#  COLOR PALETTES
# ═══════════════════════════════════════════════════════════════════════════════

COLORS = {
    'dark': {
        'bg': "#0d0f1a",
        'panel_bg': "#14172a",
        'grid': "#252840",
        'text': "#dde1f9",
        'ref': "#3a7bd5",
        'alt': "#e74c3c",
        'voi': "#f5a623",
        'missing': "#22253a",
        'samples': ["#3a7bd5", "#2ecc71", "#9b59b6", "#e67e22", "#1abc9c", "#e91e63"],
        'row_bg': ["#181b2e", "#1c1f33"],
        'hdr_bg': "#1e2240",
    },
    'light': {
        'bg': "#f8f9fa",
        'panel_bg': "#ffffff",
        'grid': "#e0e0e0",
        'text': "#1a1a1a",
        'ref': "#2196F3",
        'alt': "#d32f2f",
        'voi': "#ff9800",
        'missing': "#cccccc",
        'samples': ["#1976D2", "#388E3C", "#7B1FA2", "#F57C00", "#00838F", "#C2185B"],
        'row_bg': ["#f5f5f5", "#fafafa"],
        'hdr_bg': "#e8eaf6",
    }
}

# ═══════════════════════════════════════════════════════════════════════════════
#  UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def safe_int(val):
    """Safely convert value to int, return None if fails"""
    try:
        return int(val)
    except (TypeError, ValueError):
        return None

def haplotype_similarity(h1: np.ndarray, h2: np.ndarray) -> float:
    """
    Calculate Identity-by-State (IBS) similarity
    
    Args:
        h1, h2: Haplotype arrays with values 0 (REF), 1 (ALT), -1 (missing)
    
    Returns:
        Fraction of positions where alleles match (excluding missing)
    """
    valid = (h1 >= 0) & (h2 >= 0)
    n = valid.sum()
    if n == 0:
        return float("nan")
    return float((h1[valid] == h2[valid]).sum()) / n

def style_ax(ax, colors):
    """Apply dark theme to matplotlib axis"""
    ax.set_facecolor(colors['panel_bg'])
    for sp in ax.spines.values():
        sp.set_edgecolor(colors['grid'])

def gt_string(a1, a2, phased):
    """Format genotype as string"""
    sep = "|" if phased else "/"
    def allele_code(a):
        return "." if a == -1 else str(a)
    return allele_code(a1) + sep + allele_code(a2)

# ═══════════════════════════════════════════════════════════════════════════════
#  VCF PARSING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_vcf_region(vcf_file: str, chrom: str, start: int, end: int,
                     skip_multiallelic=True, skip_indels=True) -> list[dict]:
    """
    Extract variants from specified region in single-sample VCF
    
    Args:
        vcf_file: Path to .vcf.gz file (must have .tbi index)
        chrom: Chromosome (with or without 'chr' prefix)
        start, end: Genomic coordinates (1-based, inclusive)
        skip_multiallelic: Skip sites with multiple ALT alleles
        skip_indels: Skip indels
    
    Returns:
        List of variant dicts with genotype and quality info
    """
    records = []
    region = f"{chrom}:{start}-{end}"
    
    vcf = VCF(vcf_file)
    
    try:
        iterator = vcf(region)
    except Exception as exc:
        print(f"  [WARN] Could not query region {region} in {vcf_file}: {exc}")
        vcf.close()
        return records
    
    for v in iterator:
        # Apply filters
        if skip_multiallelic and len(v.ALT) != 1:
            continue
        
        is_indel = len(v.REF) != 1 or len(v.ALT[0]) != 1
        if skip_indels and is_indel:
            continue
        
        # Extract genotype
        gt = v.genotypes[0]  # [allele1, allele2, phased]
        a1, a2 = gt[0], gt[1]
        phased = bool(gt[2])
        
        # Extract FORMAT fields
        dp_arr = v.format("DP")
        gq_arr = v.format("GQ")
        ad_arr = v.format("AD")
        
        dp = safe_int(dp_arr[0][0]) if dp_arr is not None else None
        gq = safe_int(gq_arr[0][0]) if gq_arr is not None else None
        
        ad_ref, ad_alt = None, None
        if ad_arr is not None:
            ad_ref = safe_int(ad_arr[0][0])
            ad_alt = safe_int(ad_arr[0][1])
        
        records.append({
            'chrom': v.CHROM,
            'pos': v.POS,
            'rsid': v.ID or ".",
            'ref': v.REF,
            'alt': v.ALT[0],
            'a1': a1,
            'a2': a2,
            'phased': phased,
            'dp': dp,
            'gq': gq,
            'ad_ref': ad_ref,
            'ad_alt': ad_alt,
            'filter': ";".join(v.FILTER) if v.FILTER else "PASS",
        })
    
    vcf.close()
    return records

def auto_detect_chrom_prefix(vcf_file: str) -> str:
    """
    Auto-detect if VCF uses 'chr'-prefixed chromosomes
    
    Returns:
        'chr' if chr-prefixed, '' if numeric
    """
    try:
        vcf = VCF(vcf_file)
        chroms = list(vcf.seqnames)[:5]  # Sample first 5
        vcf.close()
        
        has_chr_prefix = any(c.startswith('chr') for c in chroms)
        return 'chr' if has_chr_prefix else ''
    except Exception as e:
        print(f"[AUTO-DETECT] Could not detect prefix: {e}")
        return ''

# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class HaplotypeAnalyzer:
    """Main analysis pipeline"""
    
    def __init__(self, config_mgr: ConfigManager):
        self.cfg = config_mgr
        self.colors = COLORS[config_mgr.get_color_scheme()]
        
        # Setup directories
        self.vcf_dir = Path(config_mgr.get_vcf_dir())
        self.results_dir = Path(config_mgr.get_results_dir())
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Get VOI parameters
        voi_params = config_mgr.get_voi_params()
        self.gene = voi_params['gene']
        self.voi_chrom = voi_params['chrom']
        self.voi_pos = voi_params['pos']
        self.voi_ref = voi_params['ref']
        self.voi_alt = voi_params['alt']
        self.voi_rsid = voi_params['rsid']
        self.voi_cdna = voi_params['cdna']
        self.voi_protein = voi_params['protein']
        self.chrom_prefix = voi_params['chrom_prefix']
        
        # Analysis parameters
        self.flank = config_mgr.get_flank()
        self.region_start = self.voi_pos - self.flank
        self.region_end = self.voi_pos + self.flank
        
        # Construct full chromosome name
        if not self.chrom_prefix:
            # Auto-detect
            sample_stem = config_mgr.get_samples()[0][1]
            vcf_path = self.vcf_dir / (sample_stem + '.hard-filtered.vcf.gz')
            detected = auto_detect_chrom_prefix(str(vcf_path))
            self.chrom_prefix = detected
            print(f"[AUTO-DETECT] Chromosome prefix: '{self.chrom_prefix}'\n")
        
        self.voi_chrom_full = f"{self.chrom_prefix}{self.voi_chrom}"
        
        # Data storage
        self.all_records = {}
        self.all_positions = []
        self.pos_meta = {}
        self.hap_matrix = {}
        self.voi_strand = {}
        self.results_df = None
        self.sim_matrix = None
        self.strand_keys = None
        
        self.skip_multiallelic = not config_mgr.get_include_multiallelic()
        self.skip_indels = not config_mgr.get_include_indels()
    
    def run(self):
        """Execute full analysis pipeline"""
        print("="*80)
        print(f"  HAPLOTYPE ANALYSIS — {self.gene}  rs{self.voi_rsid}")
        print(f"  VOI: {self.voi_chrom_full}:{self.voi_pos}  {self.voi_ref}>{self.voi_alt}")
        print(f"  HGVS: {self.voi_cdna}  {self.voi_protein}")
        print("="*80)
        
        print("\n[STEP 1] Checking VCF files...")
        self._check_vcf_files()
        
        print("\n[STEP 2] Parsing VCF regions...")
        self._parse_vcfs()
        
        print("\n[STEP 3] Building haplotype matrix...")
        self._build_haplotype_matrix()
        
        print("\n[STEP 4] Identifying VOI carriers...")
        self._identify_carriers()
        
        print("\n[STEP 5] Computing IBS similarity...")
        self._compute_ibs()
        
        print("\n[STEP 6] Creating result tables...")
        self._create_results_tables()
        
        print("\n[STEP 7] Generating visualizations...")
        self._create_visualizations()
        
        print("\n[STEP 8] Writing summary report...")
        self._write_summary()
        
        print("\n" + "="*80)
        print("  ANALYSIS COMPLETE")
        print("="*80 + "\n")
    
    def _check_vcf_files(self):
        """Verify all VCF files exist and are indexed"""
        samples = self.cfg.get_samples()
        missing = []
        
        for lims_id, stem in samples:
            vcf_file = self.vcf_dir / (stem + '.hard-filtered.vcf.gz')
            tbi_file = vcf_file.with_suffix(vcf_file.suffix + '.tbi')
            
            if not vcf_file.exists():
                missing.append(str(vcf_file))
                print(f"  [MISSING] {vcf_file}")
            elif not tbi_file.exists():
                print(f"  [WARNING] No index: {tbi_file}")
                print(f"            Run: tabix -p vcf {vcf_file}")
            else:
                size_mb = vcf_file.stat().st_size / (1024*1024)
                file_type = "WGS" if size_mb > 500 else "LCG" if size_mb < 100 else "exome"
                print(f"  [OK] {vcf_file.name} ({size_mb:.1f} MB, {file_type})")
        
        if missing:
            raise FileNotFoundError(
                f"Missing VCF files:\n  " + "\n  ".join(missing)
            )
    
    def _parse_vcfs(self):
        """Parse VCF files and extract variants in VOI region"""
        samples = self.cfg.get_samples()
        
        for lims_id, stem in samples:
            vcf_file = self.vcf_dir / (stem + '.hard-filtered.vcf.gz')
            
            recs = parse_vcf_region(
                str(vcf_file),
                self.voi_chrom_full,
                self.region_start,
                self.region_end,
                skip_multiallelic=self.skip_multiallelic,
                skip_indels=self.skip_indels
            )
            
            self.all_records[lims_id] = recs
            
            # Status check
            voi = next((r for r in recs if r['pos'] == self.voi_pos), None)
            
            if voi is None:
                voi_status = "NOT FOUND"
            elif voi['a1'] == -1:
                voi_status = "NO-CALL"
            elif voi['a1'] == 0 and voi['a2'] == 0:
                voi_status = "HOM REF (WT)"
            elif voi['a1'] == 1 and voi['a2'] == 1:
                voi_status = "HOM ALT ***"
            else:
                voi_status = "HET"
            
            n_var = len(recs)
            print(f"  {lims_id:25s}  {n_var:4d} variants  |  VOI: {voi_status}")
    
    def _build_haplotype_matrix(self):
        """Build matrix of haplotypes across all samples"""
        # Union all positions
        all_pos = set()
        for recs in self.all_records.values():
            all_pos.update(r['pos'] for r in recs)
        
        self.all_positions = sorted(all_pos)
        n_pos = len(self.all_positions)
        
        # Metadata per position
        for recs in self.all_records.values():
            for r in recs:
                if r['pos'] not in self.pos_meta:
                    self.pos_meta[r['pos']] = {
                        'rsid': r['rsid'],
                        'ref': r['ref'],
                        'alt': r['alt'],
                        'is_voi': r['pos'] == self.voi_pos,
                    }
        
        # Build haplotype matrix
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        
        for lims_id in sample_ids:
            by_pos = {r['pos']: r for r in self.all_records[lims_id]}
            h1 = np.full(n_pos, -1, dtype=int)
            h2 = np.full(n_pos, -1, dtype=int)
            
            for j, pos in enumerate(self.all_positions):
                if pos in by_pos:
                    r = by_pos[pos]
                    h1[j] = r['a1']
                    h2[j] = r['a2']
            
            self.hap_matrix[(lims_id, 1)] = h1
            self.hap_matrix[(lims_id, 2)] = h2
    
    def _identify_carriers(self):
        """Identify which haplotype carries the VOI"""
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        
        voi_idx = (
            self.all_positions.index(self.voi_pos)
            if self.voi_pos in self.all_positions else None
        )
        
        print("\n  VOI carrier status:")
        for lims_id in sample_ids:
            h1 = self.hap_matrix[(lims_id, 1)]
            h2 = self.hap_matrix[(lims_id, 2)]
            
            h1_carrier = voi_idx is not None and h1[voi_idx] == 1
            h2_carrier = voi_idx is not None and h2[voi_idx] == 1
            
            self.voi_strand[lims_id] = {'h1': h1_carrier, 'h2': h2_carrier}
            
            if h1_carrier and h2_carrier:
                status = "HOM ALT *** AFFECTED ***"
            elif h1_carrier:
                status = "HET CARRIER (Hap1)"
            elif h2_carrier:
                status = "HET CARRIER (Hap2)"
            else:
                status = "WILDTYPE"
            
            print(f"    {lims_id:25s}  {status}")
    
    def _compute_ibs(self):
        """Compute pairwise IBS similarity between all haplotypes"""
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        
        self.strand_keys = []
        for lims_id in sample_ids:
            self.strand_keys.append((lims_id, 1))
            self.strand_keys.append((lims_id, 2))
        
        n = len(self.strand_keys)
        self.sim_matrix = np.zeros((n, n))
        
        for i, k1 in enumerate(self.strand_keys):
            for j, k2 in enumerate(self.strand_keys):
                self.sim_matrix[i, j] = haplotype_similarity(
                    self.hap_matrix[k1], self.hap_matrix[k2]
                )
        
        # Print carrier-to-carrier IBS
        carrier_strands = [
            (lims_id, st)
            for lims_id in sample_ids
            for st in [1, 2]
            if self.voi_strand[lims_id][f'h{st}']
        ]
        
        if len(carrier_strands) >= 2:
            print("\n  Carrier haplotype IBS similarity:")
            for i, ks_a in enumerate(carrier_strands):
                for ks_b in carrier_strands[i+1:]:
                    sim = haplotype_similarity(
                        self.hap_matrix[ks_a], self.hap_matrix[ks_b]
                    )
                    interp = (
                        "→ IDENTICAL (shared founder)"
                        if sim >= 0.95 else
                        "→ SIMILAR (shared ancestry)"
                        if sim >= 0.70 else
                        "→ DIFFERENT (independent)"
                    )
                    print(f"    {ks_a[0]} Hap{ks_a[1]} ↔ {ks_b[0]} Hap{ks_b[1]}  "
                          f"IBS={sim:.4f}  {interp}")
    
    def _create_results_tables(self):
        """Create TSV output tables"""
        rows = []
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        
        for lims_id in sample_ids:
            by_pos = {r['pos']: r for r in self.all_records[lims_id]}
            
            for pos in self.all_positions:
                meta = self.pos_meta.get(pos, {})
                rec = by_pos.get(pos)
                
                if rec:
                    gt = gt_string(rec['a1'], rec['a2'], rec['phased'])
                    dp = rec['dp']
                    gq = rec['gq']
                    arf = (
                        round(rec['ad_alt'] / rec['dp'], 4)
                        if rec['dp'] and rec['ad_alt'] is not None else None
                    )
                else:
                    gt, dp, gq, arf = './..', None, None, None
                
                rows.append({
                    'LIMS_ID': lims_id,
                    'CHROM': self.voi_chrom_full,
                    'POS': pos,
                    'RSID': meta.get('rsid', '.'),
                    'REF': meta.get('ref', '.'),
                    'ALT': meta.get('alt', '.'),
                    'GT': gt,
                    'DP': dp,
                    'GQ': gq,
                    'ALT_AF': arf,
                    'IS_VOI': meta.get('is_voi', False),
                })
        
        self.results_df = pd.DataFrame(rows)
        
        # Save full table
        full_tsv = self.results_dir / 'haplotype_results.tsv'
        self.results_df.to_csv(full_tsv, sep='\t', index=False)
        print(f"  Full results: {full_tsv}")
        
        # Save VOI-only table
        voi_tsv = self.results_dir / 'voi_calls.tsv'
        if 'IS_VOI' in self.results_df.columns and self.results_df['IS_VOI'].any():
            voi_df = self.results_df[self.results_df['IS_VOI']].copy()
            voi_df.to_csv(voi_tsv, sep='\t', index=False)
            print(f"  VOI calls:    {voi_tsv}")
        else:
            pd.DataFrame(columns=self.results_df.columns).to_csv(voi_tsv, sep='\t', index=False)
            print(f"  VOI calls:    {voi_tsv} [EMPTY]")
    
    def _create_visualizations(self):
        """Create all panel visualizations"""
        dpi = self.cfg.get_output_dpi()
        
        print("\n  Generating panels...")
        
        # Panel A
        self._panel_a_haplotype_grid(dpi)
        
        # Panel B
        self._panel_b_ibs_heatmap(dpi)
        
        # Panel C
        self._panel_c_alt_af(dpi)
        
        # Panel D
        self._panel_d_voi_summary(dpi)
    
    def _panel_a_haplotype_grid(self, dpi):
        """Panel A: Haplotype block grid"""
        fig, ax = plt.subplots(figsize=(24, 9), facecolor=self.colors['bg'])
        style_ax(ax, self.colors)
        
        n_pos = len(self.all_positions)
        n_rows = len(self.strand_keys)
        cell_h = 0.9
        row_gap = 1.25
        
        # Background
        row_bgs = self.colors['row_bg']
        for row_i in range(n_rows):
            y = (n_rows - 1 - row_i) * row_gap
            ax.barh(y, n_pos, height=cell_h, left=-0.5,
                   color=row_bgs[row_i % 2], zorder=0)
        
        # Plot haplotypes
        voi_idx = (
            self.all_positions.index(self.voi_pos)
            if self.voi_pos in self.all_positions else None
        )
        
        for row_i, (lims_id, strand) in enumerate(self.strand_keys):
            hap = self.hap_matrix[(lims_id, strand)]
            y = (n_rows - 1 - row_i) * row_gap
            
            for col_j, pos in enumerate(self.all_positions):
                a = hap[col_j]
                is_voi = pos == self.voi_pos
                
                if a == -1:
                    color = self.colors['missing']
                elif is_voi and a == 1:
                    color = self.colors['voi']
                elif a == 1:
                    color = self.colors['alt']
                else:
                    color = self.colors['ref']
                
                ax.bar(col_j, cell_h, bottom=y - cell_h/2,
                      width=0.85, color=color,
                      linewidth=0.25, edgecolor=self.colors['bg'], zorder=2)
            
            # VOI label
            if voi_idx is not None:
                h1 = self.hap_matrix[(lims_id, 1)]
                h2 = self.hap_matrix[(lims_id, 2)]
                allele_val = h1[voi_idx] if strand == 1 else h2[voi_idx]
                label = "ALT" if allele_val == 1 else "REF"
                label_color = self.colors['voi'] if allele_val == 1 else self.colors['ref']
                
                ax.text(n_pos - 0.3, y, label, va='center', ha='left',
                       fontsize=8, color=label_color, fontweight='bold')
        
        # VOI line
        if voi_idx is not None:
            ax.axvline(voi_idx, color=self.colors['voi'], lw=2.5,
                      linestyle='--', alpha=0.6, zorder=3)
            ax.text(voi_idx, (n_rows - 0.15) * row_gap,
                   f"★ {self.voi_rsid}", color=self.colors['voi'],
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Labels
        ax.set_xticks(range(n_pos))
        x_labels = []
        for pos in self.all_positions:
            m = self.pos_meta.get(pos, {})
            rid = m.get('rsid', str(pos))
            x_labels.append('★ ' + rid if m.get('is_voi') else rid)
        ax.set_xticklabels(x_labels, rotation=50, ha='right',
                          fontsize=8, color=self.colors['text'])
        
        y_ticks = [(n_rows - 1 - i) * row_gap for i in range(n_rows)]
        y_labels = [f"{s}\nHap {st}" for s, st in self.strand_keys]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, fontsize=8, color=self.colors['text'])
        ax.tick_params(colors=self.colors['text'], which='both', length=3)
        
        # Sample separators
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        for si in range(1, len(sample_ids)):
            sep_y = (n_rows - si * 2 - 1) * row_gap + row_gap / 2
            ax.axhline(sep_y, color=self.colors['grid'], lw=1.2, alpha=0.8)
        
        ax.set_xlim(-0.7, n_pos + 0.5)
        ax.set_ylim(-cell_h, n_rows * row_gap - 0.3)
        
        title = (f"Haplotype Block Grid  |  {self.gene}  "
                f"{self.voi_chrom_full}:{self.region_start:,}–{self.region_end:,}")
        ax.set_title(title, color=self.colors['text'], fontsize=13, fontweight='bold', pad=10)
        
        # Legend
        leg = [
            mpatches.Patch(facecolor=self.colors['ref'], label='REF allele'),
            mpatches.Patch(facecolor=self.colors['alt'], label='ALT allele'),
            mpatches.Patch(facecolor=self.colors['voi'],
                          label=f"Pathogenic ALT  {self.voi_rsid}  {self.voi_cdna}  {self.voi_protein}"),
            mpatches.Patch(facecolor=self.colors['missing'], label='Missing / no-call'),
        ]
        ax.legend(handles=leg, loc='upper left', framealpha=0.25,
                 labelcolor=self.colors['text'], fontsize=8.5,
                 facecolor=self.colors['panel_bg'], edgecolor=self.colors['grid'])
        
        out_file = self.results_dir / 'panelA_haplotype_grid.png'
        fig.savefig(out_file, dpi=dpi, bbox_inches='tight', facecolor=self.colors['bg'])
        plt.close(fig)
        print(f"    ✓ Panel A: {out_file.name}")
    
    def _panel_b_ibs_heatmap(self, dpi):
        """Panel B: IBS similarity heatmap"""
        fig, ax = plt.subplots(figsize=(10, 9), facecolor=self.colors['bg'])
        style_ax(ax, self.colors)
        
        n_strands = len(self.strand_keys)
        cmap = plt.cm.RdYlGn
        
        im = ax.imshow(self.sim_matrix, cmap=cmap, vmin=0.0, vmax=1.0, aspect='auto')
        
        labels = [f"{s}\nHap{st}" for s, st in self.strand_keys]
        ax.set_xticks(range(n_strands))
        ax.set_yticks(range(n_strands))
        ax.set_xticklabels(labels, rotation=45, ha='right',
                          fontsize=7, color=self.colors['text'])
        ax.set_yticklabels(labels, fontsize=7, color=self.colors['text'])
        ax.tick_params(colors=self.colors['text'], length=2)
        
        for i in range(n_strands):
            for j in range(n_strands):
                v = self.sim_matrix[i, j]
                tcol = 'black' if 0.3 < v < 0.75 else 'white'
                ax.text(j, i, f'{v:.2f}', ha='center', va='center',
                       fontsize=8, color=tcol, fontweight='bold')
        
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('IBS Similarity', color=self.colors['text'], fontsize=8)
        cbar.ax.yaxis.set_tick_params(color=self.colors['text'], labelsize=7)
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color=self.colors['text'])
        
        ax.set_title('Pairwise IBS Similarity\n(all haplotypes)',
                    color=self.colors['text'], fontsize=11, fontweight='bold', pad=10)
        
        out_file = self.results_dir / 'panelB_ibs_heatmap.png'
        fig.savefig(out_file, dpi=dpi, bbox_inches='tight', facecolor=self.colors['bg'])
        plt.close(fig)
        print(f"    ✓ Panel B: {out_file.name}")
    
    def _panel_c_alt_af(self, dpi):
        """Panel C: ALT Allele Frequency"""
        fig, ax = plt.subplots(figsize=(16, 7), facecolor=self.colors['bg'])
        style_ax(ax, self.colors)
        
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        n_pos = len(self.all_positions)
        n_samp = len(sample_ids)
        bar_w = 0.8 / n_samp
        x_idx = np.arange(n_pos)
        
        for si, lims_id in enumerate(sample_ids):
            by_pos = {r['pos']: r for r in self.all_records[lims_id]}
            afs = []
            for pos in self.all_positions:
                rec = by_pos.get(pos)
                if rec and rec['dp'] and rec['ad_alt'] is not None:
                    afs.append(rec['ad_alt'] / rec['dp'])
                else:
                    afs.append(0.0)
            
            ax.bar(x_idx + si * bar_w, afs, width=bar_w,
                  label=f"LIMS {lims_id}",
                  color=self.colors['samples'][si % len(self.colors['samples'])],
                  alpha=0.85, edgecolor=self.colors['bg'], linewidth=0.3)
        
        # VOI highlight
        voi_idx = (
            self.all_positions.index(self.voi_pos)
            if self.voi_pos in self.all_positions else None
        )
        
        if voi_idx is not None:
            ax.axvspan(voi_idx - 0.05, voi_idx + n_samp * bar_w + 0.05,
                      alpha=0.10, color=self.colors['voi'], zorder=0)
            ax.axvline(voi_idx + (n_samp - 1) * bar_w / 2,
                      color=self.colors['voi'], lw=1.5, linestyle='--', alpha=0.7)
        
        # Labels
        rsid_labels = [self.pos_meta[p]['rsid'] for p in self.all_positions]
        ax.set_xticks(x_idx + (n_samp - 1) * bar_w / 2)
        ax.set_xticklabels(rsid_labels, rotation=50, ha='right',
                          fontsize=7.5, color=self.colors['text'])
        
        ax.set_ylim(0, 1.08)
        ax.axhline(0.5, color='#7f8c8d', lw=0.7, linestyle=':', alpha=0.6)
        ax.set_ylabel('ALT Allele Frequency', color=self.colors['text'], fontsize=9.5)
        ax.tick_params(colors=self.colors['text'], length=2)
        ax.grid(axis='y', color=self.colors['grid'], lw=0.4, alpha=0.6)
        
        ax.legend(fontsize=8, facecolor=self.colors['panel_bg'],
                 edgecolor=self.colors['grid'], labelcolor=self.colors['text'])
        
        ax.set_title('ALT Allele Frequency per Variant',
                    color=self.colors['text'], fontsize=11, fontweight='bold', pad=10)
        
        out_file = self.results_dir / 'panelC_alt_af.png'
        fig.savefig(out_file, dpi=dpi, bbox_inches='tight', facecolor=self.colors['bg'])
        plt.close(fig)
        print(f"    ✓ Panel C: {out_file.name}")
    
    def _panel_d_voi_summary(self, dpi):
        """Panel D: VOI call summary table"""
        fig, ax = plt.subplots(figsize=(18, 5), facecolor=self.colors['bg'])
        style_ax(ax, self.colors)
        ax.axis('off')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        
        sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
        tbl_rows = []
        
        for lims_id in sample_ids:
            by_pos = {r['pos']: r for r in self.all_records[lims_id]}
            rec = by_pos.get(self.voi_pos)
            
            if rec:
                a1v, a2v = rec['a1'], rec['a2']
                ph = '|' if rec['phased'] else '/'
                gt = f"{a1v}{ph}{a2v}"
                
                if a1v == 0 and a2v == 0:
                    zyg = "Hom REF"
                elif a1v == 1 and a2v == 1:
                    zyg = "Hom ALT !!!"
                elif a1v == -1 or a2v == -1:
                    zyg = "No-call"
                else:
                    zyg = "Heterozygous"
                
                dp_ = str(rec['dp']) if rec['dp'] else "—"
                gq_ = str(rec['gq']) if rec['gq'] else "—"
                ad_ = (f"{rec['ad_ref']},{rec['ad_alt']}"
                      if rec['ad_ref'] is not None else "—")
            else:
                gt, zyg, dp_, gq_, ad_ = "./.", "Not called", "—", "—", "—"
            
            tbl_rows.append([
                lims_id, f"{self.voi_chrom_full}:{self.voi_pos}",
                self.voi_rsid, f"{self.voi_ref}>{self.voi_alt}",
                gt, zyg, dp_, gq_, ad_,
                f"{self.voi_cdna} {self.voi_protein}"
            ])
        
        col_hdrs = ["LIMS ID", "Position", "rsID", "Change", "GT",
                   "Zygosity", "DP", "GQ", "AD", "HGVS"]
        n_c = len(col_hdrs)
        col_w = 1.0 / n_c
        row_h = 0.22
        
        # Header
        for ci, hdr in enumerate(col_hdrs):
            rect = mpatches.FancyBboxPatch(
                (ci * col_w, 1.0 - row_h), col_w, row_h,
                boxstyle='square,pad=0', lw=0.5,
                edgecolor=self.colors['grid'],
                facecolor=self.colors['hdr_bg'],
                transform=ax.transAxes, clip_on=False, zorder=2
            )
            ax.add_patch(rect)
            ax.text((ci + 0.5) * col_w, 1.0 - row_h/2, hdr,
                   ha='center', va='center', fontsize=9,
                   fontweight='bold', color=self.colors['text'],
                   transform=ax.transAxes, zorder=3)
        
        # Rows
        zyg_col_map = {
            'Heterozygous': self.colors['voi'],
            'Hom ALT !!!': self.colors['alt'],
            'Hom REF': self.colors['ref'],
            'No-call': '#95a5a6',
            'Not called': '#95a5a6',
        }
        
        row_bgs = self.colors['row_bg']
        
        for ri, row_vals in enumerate(tbl_rows):
            yb = 1.0 - row_h * (ri + 2)
            bg = row_bgs[ri % 2]
            
            for ci, val in enumerate(row_vals):
                rect = mpatches.FancyBboxPatch(
                    (ci * col_w, yb), col_w, row_h,
                    boxstyle='square,pad=0', lw=0.4,
                    edgecolor=self.colors['grid'],
                    facecolor=bg,
                    transform=ax.transAxes, clip_on=False, zorder=2
                )
                ax.add_patch(rect)
                
                fc = zyg_col_map.get(val, self.colors['text']) if ci == 5 else self.colors['text']
                ax.text((ci + 0.5) * col_w, yb + row_h/2, val,
                       ha='center', va='center', fontsize=8.5,
                       color=fc, transform=ax.transAxes, zorder=3)
        
        title = (f"VOI Call Summary  |  {self.voi_chrom_full}:{self.voi_pos}  "
                f"{self.voi_ref}>{self.voi_alt}  {self.voi_cdna}  {self.voi_protein}")
        ax.set_title(title, color=self.colors['text'], fontsize=11, fontweight='bold', pad=10)
        
        out_file = self.results_dir / 'panelD_voi_summary.png'
        fig.savefig(out_file, dpi=dpi, bbox_inches='tight', facecolor=self.colors['bg'])
        plt.close(fig)
        print(f"    ✓ Panel D: {out_file.name}")
    
    def _write_summary(self):
        """Write text summary report"""
        report_file = self.results_dir / 'analysis_summary.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write(f"HAPLOTYPE ANALYSIS SUMMARY\n")
            f.write(f"Gene: {self.gene}\n")
            f.write(f"VOI: {self.voi_chrom_full}:{self.voi_pos}  {self.voi_ref}>{self.voi_alt}\n")
            f.write(f"HGVS: {self.voi_cdna}  {self.voi_protein}\n")
            f.write(f"rsID: {self.voi_rsid}\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"Analysis Parameters:\n")
            f.write(f"  Region: {self.voi_chrom_full}:{self.region_start:,}–{self.region_end:,}\n")
            f.write(f"  Variants in region: {len(self.all_positions)}\n")
            f.write(f"  Flank size: {self.flank:,} bp each side\n\n")
            
            sample_ids = [lims_id for lims_id, _ in self.cfg.get_samples()]
            
            f.write("Sample VOI Status:\n")
            for lims_id in sample_ids:
                info = self.voi_strand[lims_id]
                if info['h1'] and info['h2']:
                    call = "HOM ALT *** AFFECTED ***"
                elif info['h1']:
                    call = "HET CARRIER (Hap1)"
                elif info['h2']:
                    call = "HET CARRIER (Hap2)"
                else:
                    call = "WILDTYPE"
                f.write(f"  {lims_id:25s}  {call}\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("Output Files:\n")
            f.write(f"  {(self.results_dir / 'panelA_haplotype_grid.png').name}\n")
            f.write(f"  {(self.results_dir / 'panelB_ibs_heatmap.png').name}\n")
            f.write(f"  {(self.results_dir / 'panelC_alt_af.png').name}\n")
            f.write(f"  {(self.results_dir / 'panelD_voi_summary.png').name}\n")
            f.write(f"  haplotype_results.tsv\n")
            f.write(f"  voi_calls.tsv\n")
            f.write("="*80 + "\n")
        
        print(f"  Summary: {report_file}")

# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 2:
        config_file = 'config.ini'
        if not os.path.exists(config_file):
            print("\nUsage:")
            print("  python3 haplotype_analyzer.py config.ini\n")
            print("Creating template config.ini...\n")
            
            # Show template
            print(open(__file__).read().split('# ═══════════════════')[1].split('# ═══════════════════')[0])
            sys.exit(1)
    else:
        config_file = sys.argv[1]
    
    try:
        cfg_mgr = ConfigManager(config_file)
        analyzer = HaplotypeAnalyzer(cfg_mgr)
        analyzer.run()
    except Exception as e:
        print(f"\n[ERROR] {e}\n")
        sys.exit(1)

if __name__ == '__main__':
    main()
