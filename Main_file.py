#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genomic Sequence Analysis — Mouse Genome (mm9)
===============================================
Performs a series of analyses on selected chromosomes from the mm9 mouse
genome build, covering gene counting, sequence extraction, coding/non-coding
quantification, and TATA motif detection near transcription start sites.

Author : Benjamin Garcés Cifuentes
Email  : agarces2381@gmail.com
Python : 3.12
"""

# ── Imports ───────────────────────────────────────────────────────────────────
# Standard library
import csv
import gc
import io
import zipfile
from datetime import date

# Third-party  (pip install numpy matplotlib)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Local module
from LoadFASTA_Function import LoadFastaFile, LoadGene, TSSChroms


# ── 1. Load gene annotations ──────────────────────────────────────────────────
# LoadGene parses the UCSC knownGene table into a dictionary keyed by
# transcript identifier. Each entry holds chromosome, strand, and
# transcript start/end coordinates.
GENE_FILE  = 'mm9_sel_chroms_knownGene.txt'
FASTA_FILE = 'selChroms_mm9.fa.zip'
REPORT_OUT = 'genomic_analysis_report.txt'

gene_inf = LoadGene(GENE_FILE)

# Quick access example — chromosome of uc009auw.1:
gene_inf['uc009auw.1']['chr']


# ── 2. Identify chromosomes and count genes per chromosome ────────────────────
# Collect unique chromosome names, then tally transcripts per chromosome.
# Using a set avoids repeated membership checks on a growing list.
chroms = sorted(set(gene_inf[g]['chr'] for g in gene_inf))

gene_counts = {}
for chrom in chroms:
    gene_counts[chrom] = sum(1 for g in gene_inf if gene_inf[g]['chr'] == chrom)


# ── 3. Load raw chromosomal sequences ────────────────────────────────────────
# Reads the zipped FASTA file into a dictionary of strings.
# This is the most memory-intensive step; expect 1-3 minutes on a typical
# workstation depending on available RAM.
seq_dict = LoadFastaFile(FASTA_FILE)

# Length of chromosome 6 as a quick sanity check:
len(seq_dict['chr6'])


# ── 4. Extract a gene sequence — Cntn4 (uc009dcr.2) ─────────────────────────
# Retrieves the full transcript sequence of Cntn4 by slicing the chromosomal
# string with the annotated start and end coordinates.
CNTN4_ID = 'uc009dcr.2'

cntn4_chr   = gene_inf[CNTN4_ID]['chr']
cntn4_start = gene_inf[CNTN4_ID]['start']
cntn4_end   = gene_inf[CNTN4_ID]['end']
cntn4_seq   = seq_dict[cntn4_chr][cntn4_start:cntn4_end]

# Inspect a short window and locate the first ATG start codon:
cntn4_seq[5:200]
cntn4_atg = cntn4_seq.index('ATG')


# ── 5. Gene length distribution ───────────────────────────────────────────────
# Calculates the span (txEnd - txStart) for every transcript. Note that this
# reflects the full pre-mRNA length including introns, not just coding exons.
gene_lengths = {
    g: np.absolute(gene_inf[g]['end'] - gene_inf[g]['start'])
    for g in gene_inf
}

lengths_array = np.array(list(gene_lengths.values()))

# ── Histogram ─────────────────────────────────────────────────────────────────
# Gene lengths span several orders of magnitude, so the x-axis is plotted on
# a log10 scale. Bins are computed in log space for uniform visual density.
# The y-axis is also log-scaled to prevent the tall short-gene peak from
# crushing the long-tail signal. A vertical line marks the median.
fig, ax = plt.subplots(figsize=(9, 4.5))

log_lengths = np.log10(lengths_array[lengths_array > 0])

ax.hist(log_lengths, bins=80, color='#2E6F8E', alpha=0.85,
        edgecolor='white', linewidth=0.3)

ax.set_xlabel('Gene length (bp)', fontsize=11)
ax.set_ylabel('Number of transcripts', fontsize=11)
ax.set_title('Gene length distribution — mm9 selected chromosomes', fontsize=12, pad=12)

# Restore human-readable bp values on the log-scaled x-axis
ax.xaxis.set_major_formatter(
    ticker.FuncFormatter(lambda val, _: f'{10**val:,.0f}')
)
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

# Median annotation
median_log = np.median(log_lengths)
ax.axvline(median_log, color='#E05A2B', linewidth=1.4,
           linestyle='--', label=f'Median: {10**median_log:,.0f} bp')
ax.legend(fontsize=10)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('gene_length_distribution.png', dpi=150)
plt.show()

# Example: length difference between two transcripts
gene_lengths['uc012enb.1'] - len(cntn4_seq)


# ── 6. Coding vs. non-coding DNA on chromosome 6 ─────────────────────────────
# A boolean NumPy array marks every base position covered by at least one
# annotated transcript. Overlapping transcripts are handled correctly because
# setting an already-True index to True is idempotent — no double counting.
chr6_max     = max(gene_inf[g]['end'] for g in gene_inf if gene_inf[g]['chr'] == 'chr6')
ingene_numpy = np.zeros(chr6_max, dtype=bool)

for gene in gene_inf:
    if gene_inf[gene]['chr'] == 'chr6':
        ingene_numpy[gene_inf[gene]['start']:gene_inf[gene]['end']] = True

coding_bp    = ingene_numpy.sum()
total_bp     = len(ingene_numpy)
noncoding_bp = total_bp - coding_bp
coding_frac  = coding_bp / total_bp

print(f'Chr6 — coding bp     : {coding_bp:,}')
print(f'Chr6 — non-coding bp : {noncoding_bp:,}')
print(f'Chr6 — coding fraction: {coding_frac:.4f}')


# ── 7. TATA motif search near transcription start sites ───────────────────────
# A 40 bp window immediately upstream of each TSS (downstream for minus-strand
# genes) is searched for the canonical TATA box sequence. Position is recorded
# as the offset within the window, which approximates distance from the TSS.
chr6_starts = TSSChroms(gene_inf, 'chr6')
WINDOW      = 40   # bp to search upstream/downstream of each TSS

tata_dis = {}
for g in chr6_starts:
    tss    = chr6_starts[g]
    strand = gene_inf[g]['strand']
    if strand == '+':
        window_seq = seq_dict['chr6'][tss - WINDOW : tss].upper()
        if 'TATA' in window_seq:
            tata_dis[g] = window_seq.rindex('TATA')
    else:
        window_seq = seq_dict['chr6'][tss : tss + WINDOW].upper()
        if 'TATA' in window_seq:
            tata_dis[g] = window_seq.index('TATA')

print(f'Genes with TATA motif within {WINDOW} bp of TSS: {len(tata_dis)}')

if tata_dis:
    mean_tata_dist = sum(tata_dis.values()) / len(tata_dis)
    print(f'Mean TATA distance from TSS: {mean_tata_dist:.2f} bp')
else:
    print('No TATA motifs found within the search window — '
          'mean distance cannot be calculated.')


# ── 8. Report output ──────────────────────────────────────────────────────────
# Writes a plain-text summary of all computed results. Useful as a portable
# record attached to a notebook, thesis appendix, or repository README.

def write_report(path):
    """Write a structured plain-text summary of all analysis results."""
    sep = '─' * 60
    lines = [
        'GENOMIC SEQUENCE ANALYSIS REPORT',
        f'Organism : Mus musculus (mm9)',
        f'Date     : {date.today().isoformat()}',
        sep,
        '',
        '1. DATASET OVERVIEW',
        f'   Gene annotation file : {GENE_FILE}',
        f'   FASTA sequence file  : {FASTA_FILE}',
        f'   Total transcripts    : {len(gene_inf):,}',
        f'   Chromosomes included : {", ".join(chroms)}',
        '',
        '2. GENE COUNTS PER CHROMOSOME',
    ]
    for chrom, count in sorted(gene_counts.items()):
        lines.append(f'   {chrom:<8} {count:>6,} transcripts')

    lines += [
        '',
        '3. GENE LENGTH STATISTICS (full transcript span, all chromosomes)',
        f'   Mean   : {lengths_array.mean():>12,.0f} bp',
        f'   Median : {np.median(lengths_array):>12,.0f} bp',
        f'   Min    : {lengths_array.min():>12,} bp',
        f'   Max    : {lengths_array.max():>12,} bp',
        f'   Histogram saved to : gene_length_distribution.png',
        '',
        '4. CODING vs. NON-CODING DNA — chromosome 6',
        f'   Coding bp       : {coding_bp:>14,}',
        f'   Non-coding bp   : {noncoding_bp:>14,}',
        f'   Total span      : {total_bp:>14,}',
        f'   Coding fraction : {coding_frac:>14.4f}  ({coding_frac*100:.2f}%)',
        '',
        '5. TATA MOTIF ANALYSIS — chromosome 6',
        f'   Search window   : {WINDOW} bp upstream/downstream of TSS',
        f'   Genes with TATA : {len(tata_dis):,}  of  {len(chr6_starts):,} chr6 transcripts',
    ]
    if tata_dis:
        lines.append(f'   Mean TATA dist  : {mean_tata_dist:.2f} bp from TSS')
    else:
        lines.append('   Mean TATA dist  : N/A (no motifs found)')

    lines += [
        '',
        '6. EXAMPLE GENE — Cntn4 (uc009dcr.2)',
        f'   Chromosome  : {cntn4_chr}',
        f'   Start       : {cntn4_start:,}',
        f'   End         : {cntn4_end:,}',
        f'   Length      : {len(cntn4_seq):,} bp',
        f'   First ATG   : position {cntn4_atg} within transcript',
        '',
        sep,
        'End of report',
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'Report written to: {path}')

write_report(REPORT_OUT)
