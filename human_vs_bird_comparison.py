#!/usr/bin/env python3
"""
Compare bird vs human adaptive-immune receptor loci.

Loci
    Immunoglobulin : IGH, IGL
    T-cell receptor: TRA, TRB, TRG   (bird TRD has no human counterpart -> skipped)

Analyses
  1. V-gene spacing   -- distance between consecutive V genes within the same
                         haplotype, main contig only.  (all loci)
  2. V-gene packing   -- V-cluster span and V-gene density (V per Mb of span).
                         Controls for the fact that avian loci are physically
                         smaller: denser packing shows up as higher density even
                         when the locus is small.  (all loci)
  3. V-to-D distance  -- gap between the V-gene cluster and the D-gene cluster.
                         IGH only: birds have D genes for IGH only.
  4. Cysteine content -- cysteine content of V-gene sequences (all loci) and of
                         IGH D-gene sequences, birds vs human.

Main-contig selection
    IGH / IGL : from summary_features.csv (highest-NumV contig per haplotype).
    TCR loci  : not in summary_features -> the contig carrying the most V genes
                per haplotype x locus in gene_list.csv (reproduces the summary
                choice for 753/755 IGH/IGL haplotypes; the 2 exceptions are ties).

Bird inputs
    summary_features.csv                       main contig for IGH/IGL
    gene_list.csv                              all V genes (Pos, Sequence, Locus, ...)
    {Order}/{Species}/{Haplotype}/IGHD.csv     IGH D genes

Human inputs (single T2T haplotype, one contig per locus)
    human_{IGH,IGL,TRA,TRB,TRG}.csv            V/D/J genes with Pos, Sequence

Usage
    python human_vs_bird_comparison.py \
        -i /local/storage/kav67/clean_birds \
        -H /local/storage/kav67/primate_t2t_igtr/data_human_t2t_gene_positions \
        --outdir human_bird_comparison --figdir figures
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Ordered so Ig loci come first, then TCR loci. Used for V-based analyses.
V_LOCI = ['IGH', 'IGL', 'TRA', 'TRB', 'TRG']


# ---------------------------------------------------------------------------
# Translation helpers
# ---------------------------------------------------------------------------
_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_frame(seq, frame):
    seq = seq.upper()
    return ''.join(_CODON_TABLE.get(seq[i:i + 3], 'X')
                   for i in range(frame, len(seq) - 2, 3))


def cysteine_fraction(seq):
    """Fraction of translated residues that are cysteine, pooled over the three
    forward reading frames (frame-agnostic, since the coding frame is unknown)."""
    if not isinstance(seq, str) or len(seq) < 3:
        return np.nan
    n_cys = n_res = 0
    for f in (0, 1, 2):
        prot = translate_frame(seq, f)
        n_cys += prot.count('C')
        n_res += len(prot)
    return n_cys / n_res if n_res else np.nan


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def build_main_contigs(summ, genes):
    """One main contig per haplotype x locus, for every locus in gene_list.

    IGH/IGL take the summary_features choice; other loci take the contig with the
    most V genes."""
    ig = (summ.sort_values('NumV', ascending=False)
              .groupby(['Order', 'Species', 'Haplotype', 'Locus'], as_index=False)
              .first()
              .rename(columns={'Contig': 'MainContig'})
          [['Order', 'Species', 'Haplotype', 'Locus', 'MainContig']])

    cnt = genes.groupby(['Order', 'Species', 'Haplotype', 'Locus', 'Contig']).size()
    tcr = (cnt.rename('n').reset_index()
              .sort_values('n', ascending=False)
              .groupby(['Order', 'Species', 'Haplotype', 'Locus'], as_index=False)
              .first()
              .rename(columns={'Contig': 'MainContig'}))
    tcr = tcr[~tcr['Locus'].isin(['IGH', 'IGL'])]
    return pd.concat([ig, tcr[['Order', 'Species', 'Haplotype', 'Locus', 'MainContig']]],
                     ignore_index=True)


def load_bird_v_genes(input_dir):
    """V genes on the main contig of each haplotype x locus (all loci)."""
    summ = pd.read_csv(os.path.join(input_dir, 'summary_features.csv'))
    genes = pd.read_csv(os.path.join(input_dir, 'gene_list.csv'))
    src = genes['Source'].str.split('/', expand=True)
    genes['Order'], genes['Species'], genes['Haplotype'] = src[0], src[1], src[2]

    main = build_main_contigs(summ, genes)
    merged = genes.merge(main, on=['Order', 'Species', 'Haplotype', 'Locus'], how='inner')
    merged = merged[merged['Contig'] == merged['MainContig']].copy()
    merged = merged[merged['Locus'].isin(V_LOCI)]
    return merged, main


def load_bird_d_genes(input_dir, main):
    """IGH D genes per haplotype, restricted to the IGH main contig."""
    igh_main = main[main['Locus'] == 'IGH']
    rows = []
    for _, r in igh_main.iterrows():
        path = os.path.join(input_dir, r['Order'], r['Species'], r['Haplotype'], 'IGHD.csv')
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path)
        d = d[d['Contig'] == r['MainContig']].copy()
        if d.empty:
            continue
        d['Order'], d['Species'], d['Haplotype'] = r['Order'], r['Species'], r['Haplotype']
        rows.append(d)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def load_human(human_dir):
    out = {}
    for locus in V_LOCI:
        p = os.path.join(human_dir, f'human_{locus}.csv')
        if os.path.exists(p):
            out[locus] = pd.read_csv(p)
        else:
            print(f'[WARN] missing human file for {locus}: {p}')
    return out


# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------
def compute_vv_distances(bird_v, human):
    rec = []
    for (o, sp, hap, locus), sub in bird_v[bird_v.GeneType == 'V'].groupby(
            ['Order', 'Species', 'Haplotype', 'Locus']):
        pos = np.sort(sub['Pos'].values)
        for d in np.diff(pos):
            rec.append({'group': 'Bird', 'Locus': locus, 'id': f'{sp}/{hap}', 'dist': float(d)})
    for locus, df in human.items():
        pos = np.sort(df[df.GeneType == 'V']['Pos'].values)
        for d in np.diff(pos):
            rec.append({'group': 'Human', 'Locus': locus, 'id': 'human_T2T', 'dist': float(d)})
    return pd.DataFrame(rec)


def compute_density(bird_v, human):
    rec = []
    for (o, sp, hap, locus), sub in bird_v[bird_v.GeneType == 'V'].groupby(
            ['Order', 'Species', 'Haplotype', 'Locus']):
        n = len(sub); span = sub['Pos'].max() - sub['Pos'].min()
        if n < 2 or span <= 0:
            continue
        rec.append({'group': 'Bird', 'Locus': locus, 'id': f'{sp}/{hap}',
                    'n_v': n, 'span_bp': span, 'density_per_mb': n / (span / 1e6)})
    for locus, df in human.items():
        v = df[df.GeneType == 'V']; span = v['Pos'].max() - v['Pos'].min()
        if span <= 0:
            continue
        rec.append({'group': 'Human', 'Locus': locus, 'id': 'human_T2T',
                    'n_v': len(v), 'span_bp': span, 'density_per_mb': len(v) / (span / 1e6)})
    return pd.DataFrame(rec)


def cluster_gap(v_pos, d_pos):
    v_lo, v_hi = v_pos.min(), v_pos.max()
    d_lo, d_hi = d_pos.min(), d_pos.max()
    if d_lo > v_hi:
        return d_lo - v_hi
    if v_lo > d_hi:
        return v_lo - d_hi
    return 0.0


def compute_vd_gap(bird_v, bird_d, human):
    rec = []
    igh_v = bird_v[(bird_v.Locus == 'IGH') & (bird_v.GeneType == 'V')]
    for (o, sp, hap), vsub in igh_v.groupby(['Order', 'Species', 'Haplotype']):
        dsub = bird_d[(bird_d.Order == o) & (bird_d.Species == sp) & (bird_d.Haplotype == hap)]
        if dsub.empty:
            continue
        rec.append({'group': 'Bird', 'id': f'{sp}/{hap}',
                    'gap': cluster_gap(vsub['Pos'].values, dsub['Pos'].values)})
    hdf = human['IGH']
    hv, hd = hdf[hdf.GeneType == 'V'], hdf[hdf.GeneType == 'D']
    if len(hv) and len(hd):
        rec.append({'group': 'Human', 'id': 'human_T2T',
                    'gap': cluster_gap(hv['Pos'].values, hd['Pos'].values)})
    return pd.DataFrame(rec)


def compute_cysteine(bird_v, bird_d, human):
    rec = []
    for _, r in bird_v[bird_v.GeneType == 'V'].iterrows():
        cf = cysteine_fraction(r['Sequence'])
        if not np.isnan(cf):
            rec.append({'group': 'Bird', 'gene': 'V', 'Locus': r['Locus'], 'cys_frac': cf})
    for _, r in bird_d.iterrows():
        cf = cysteine_fraction(r['Sequence'])
        if not np.isnan(cf):
            rec.append({'group': 'Bird', 'gene': 'D', 'Locus': 'IGH', 'cys_frac': cf})
    for locus, df in human.items():
        for gene in ('V', 'D'):
            for _, r in df[df.GeneType == gene].iterrows():
                cf = cysteine_fraction(r['Sequence'])
                if not np.isnan(cf):
                    rec.append({'group': 'Human', 'gene': gene, 'Locus': locus, 'cys_frac': cf})
    return pd.DataFrame(rec)


# ---------------------------------------------------------------------------
# Plotting
#   Bird = distribution (violin + stat line); Human T2T = single reference point.
#   All summary numbers live in the x-axis tick labels, never on top of the data.
# ---------------------------------------------------------------------------
COLORS = {'Bird': '#8ab4d6', 'Human': '#e08a4e'}
BIRD_LINE = '#22425f'
plt.rcParams.update({'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.6,
                     'axes.axisbelow': True, 'font.size': 10})


def _mwu(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return mannwhitneyu(a, b, alternative='two-sided').pvalue


def _p_label(p):
    if np.isnan(p):
        return 'n/a'
    return 'p < 1e-4' if p < 1e-4 else f'p = {p:.3g}'


def _log_yticks(ax, lo, hi):
    ticks = list(range(lo, hi + 1))
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{10.0**t:g}' for t in ticks])


def _fmt(v, unit):
    if np.isnan(v):
        return 'n/a'
    if unit == 'kb':
        return f'{v:.1f} kb' if v < 100 else f'{v:.0f} kb'
    if unit == '%':
        return f'{v:.2f}%'
    return f'{v:.0f} {unit}'.strip()


def cell(ax, bird, human, stat, log, unit, ylo=None, yhi=None, rng=None):
    """One Bird-vs-Human comparison in a single Axes. `stat` is 'median'|'mean'.
    Returns the two-sided Mann-Whitney p between bird and human values."""
    rng = rng or np.random.default_rng(0)
    statf = np.median if stat == 'median' else np.mean
    bird = np.asarray(bird, float); bird = bird[np.isfinite(bird)]
    human = np.asarray(human, float); human = human[np.isfinite(human)]
    if log:
        bird = bird[bird > 0]; human = human[human > 0]
    tb = np.log10(bird) if log else bird

    parts = ax.violinplot([tb], positions=[1], widths=0.85, showextrema=False)
    parts['bodies'][0].set_facecolor(COLORS['Bird']); parts['bodies'][0].set_alpha(0.7)
    parts['bodies'][0].set_edgecolor(BIRD_LINE); parts['bodies'][0].set_linewidth(0.8)
    sb = statf(bird)
    ax.hlines(np.log10(sb) if log else sb, 0.62, 1.38, color=BIRD_LINE, lw=2.2, zorder=6)

    sh = statf(human) if len(human) else np.nan
    if len(human) == 1:                       # single T2T reference value
        y = np.log10(human[0]) if log else human[0]
        ax.scatter([2], [y], s=130, marker='D', color=COLORS['Human'],
                   edgecolor='k', linewidth=0.8, zorder=6)
        ax.hlines(y, 0.5, 2.5, color=COLORS['Human'], ls='--', lw=1, alpha=0.6, zorder=2)
    else:                                     # human distribution (>1 value)
        th = np.log10(human) if log else human
        p2 = ax.violinplot([th], positions=[2], widths=0.85, showextrema=False)
        p2['bodies'][0].set_facecolor(COLORS['Human']); p2['bodies'][0].set_alpha(0.7)
        ax.hlines(np.log10(sh) if log else sh, 1.62, 2.38, color='#7a3d12', lw=2.2, zorder=6)

    ax.set_xlim(0.4, 2.6)
    ax.set_xticks([1, 2])
    tag = 'med' if stat == 'median' else 'mean'
    ax.set_xticklabels([f'Bird\nn={len(bird)}\n{tag} {_fmt(sb, unit)}',
                        f'Human\nn={len(human)}\n{_fmt(sh, unit)}'], fontsize=8.5)
    if log:
        _log_yticks(ax, ylo, yhi)
    ax.tick_params(axis='x', length=0)
    return _mwu(bird, human)


def _save(fig, figdir, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(figdir, f'{name}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'[FIG] {name}.png / .pdf')


def plot_vv(vv, figdir):
    loci = [l for l in V_LOCI if l in vv.Locus.unique()]
    fig, axes = plt.subplots(1, len(loci), figsize=(3.0 * len(loci), 4.8), sharey=True)
    for ax, locus in zip(np.atleast_1d(axes), loci):
        b = vv[(vv.Locus == locus) & (vv.group == 'Bird')]['dist'] / 1e3
        h = vv[(vv.Locus == locus) & (vv.group == 'Human')]['dist'] / 1e3
        p = cell(ax, b, h, 'median', log=True, unit='kb', ylo=-3, yhi=4)
        ax.set_title(f'{locus}\n{_p_label(p)}', fontsize=11)
    np.atleast_1d(axes)[0].set_ylabel('Distance between consecutive\nV genes (kb, log scale)')
    fig.suptitle('V-gene spacing: birds vs human', fontweight='bold', fontsize=13)
    fig.tight_layout()
    _save(fig, figdir, 'compare_vv_spacing')


def plot_density(dens, figdir):
    loci = [l for l in V_LOCI if l in dens.Locus.unique()]
    fig, axes = plt.subplots(2, len(loci), figsize=(3.0 * len(loci), 8.2))
    rows = [('span_bp', 1e3, 'V-cluster span\n(kb, log scale)', 0, 4),
            ('density_per_mb', 1.0, 'V-gene density\n(V per Mb, log scale)', 1, 5)]
    for ri, (col, scale, ylab, ylo, yhi) in enumerate(rows):
        row_axes = axes[ri] if len(loci) > 1 else [axes[ri]]
        for ci, locus in enumerate(loci):
            ax = row_axes[ci]
            b = dens[(dens.Locus == locus) & (dens.group == 'Bird')][col] / scale
            h = dens[(dens.Locus == locus) & (dens.group == 'Human')][col] / scale
            unit = 'kb' if col == 'span_bp' else 'V/Mb'
            cell(ax, b, h, 'median', log=True, unit=unit, ylo=ylo, yhi=yhi)
            if ri == 0:
                ax.set_title(locus, fontsize=11, fontweight='bold')
            if ci == 0:
                ax.set_ylabel(ylab)
    fig.suptitle('V-gene packing: bird loci are smaller, so density is the fair '
                 'comparison', fontweight='bold', fontsize=13)
    fig.tight_layout()
    _save(fig, figdir, 'compare_vgene_density')


def plot_vd(vd, figdir):
    fig, ax = plt.subplots(figsize=(5, 5))
    b = np.maximum(vd[vd.group == 'Bird']['gap'].values / 1e3, 0.05)
    h = np.maximum(vd[vd.group == 'Human']['gap'].values / 1e3, 0.05)
    cell(ax, b, h, 'median', log=True, unit='kb', ylo=-2, yhi=5)
    ax.set_ylabel('V-cluster to D-cluster gap (kb, log scale)')
    ax.set_title('IGH V-to-D distance: birds vs human', fontweight='bold', fontsize=12)
    fig.tight_layout()
    _save(fig, figdir, 'compare_vd_gap')


def plot_cysteine(cys, figdir):
    panels = [('V', l) for l in V_LOCI if l in cys[cys.gene == 'V'].Locus.unique()]
    panels.append(('D', 'IGH'))
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.4 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (gene, locus) in zip(axes, panels):
        b = cys[(cys.gene == gene) & (cys.Locus == locus) & (cys.group == 'Bird')]['cys_frac'] * 100
        h = cys[(cys.gene == gene) & (cys.Locus == locus) & (cys.group == 'Human')]['cys_frac'] * 100
        p = cell(ax, b, h, 'mean', log=False, unit='%')
        ax.set_title(f'{gene} genes ({locus})\n{_p_label(p)}', fontsize=11)
        ax.set_ylabel('Cysteine content (% residues, 3-frame)')
    for ax in axes[len(panels):]:
        ax.axis('off')
    fig.suptitle('Cysteine content: birds vs human', fontweight='bold', fontsize=13)
    fig.tight_layout()
    _save(fig, figdir, 'compare_cysteine')


# ---------------------------------------------------------------------------
def summarize(vv, dens, vd, cys):
    print('\n================ SUMMARY ================')
    for locus in V_LOCI:
        b = vv[(vv.Locus == locus) & (vv.group == 'Bird')]['dist'] / 1e3
        h = vv[(vv.Locus == locus) & (vv.group == 'Human')]['dist'] / 1e3
        if len(b) and len(h):
            print(f'V-V spacing {locus}: bird med {b.median():.1f} kb (n={len(b)}), '
                  f'human med {h.median():.1f} kb (n={len(h)}), {_p_label(_mwu(b, h))}')
    for locus in V_LOCI:
        b = dens[(dens.Locus == locus) & (dens.group == 'Bird')]
        h = dens[(dens.Locus == locus) & (dens.group == 'Human')]
        if len(b) and len(h):
            print(f'V-cluster {locus}: bird span med {b.span_bp.median()/1e3:.0f} kb / '
                  f'{b.n_v.median():.0f} V / {b.density_per_mb.median():.0f} V/Mb  vs  '
                  f'human {h.span_bp.iloc[0]/1e3:.0f} kb / {h.n_v.iloc[0]:.0f} V / '
                  f'{h.density_per_mb.iloc[0]:.0f} V/Mb')
    b = vd[vd.group == 'Bird']['gap'] / 1e3
    h = vd[vd.group == 'Human']['gap'] / 1e3
    print(f'V-D gap IGH: bird med {b.median():.1f} kb (n={len(b)}), '
          f'human {h.iloc[0]:.1f} kb' if len(h) else '')
    for gene, locus in [('V', l) for l in V_LOCI] + [('D', 'IGH')]:
        b = cys[(cys.gene == gene) & (cys.Locus == locus) & (cys.group == 'Bird')]['cys_frac'] * 100
        h = cys[(cys.gene == gene) & (cys.Locus == locus) & (cys.group == 'Human')]['cys_frac'] * 100
        if len(b) and len(h):
            print(f'Cysteine {gene}/{locus}: bird mean {b.mean():.2f}% (n={len(b)}), '
                  f'human mean {h.mean():.2f}% (n={len(h)}), {_p_label(_mwu(b, h))}')
    print('=========================================\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input_dir', default='/local/storage/kav67/clean_birds')
    ap.add_argument('-H', '--human_dir',
                    default='/local/storage/kav67/primate_t2t_igtr/data_human_t2t_gene_positions')
    ap.add_argument('--outdir', default='human_bird_comparison')
    ap.add_argument('--figdir', default='figures')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    print('Loading bird V genes...')
    bird_v, main = load_bird_v_genes(args.input_dir)
    for locus in V_LOCI:
        sub = bird_v[bird_v.Locus == locus]
        print(f'  {locus}: {len(sub)} V genes, {sub.groupby(["Species","Haplotype"]).ngroups} haplotypes')
    print('Loading bird IGH D genes...')
    bird_d = load_bird_d_genes(args.input_dir, main)
    print(f'  {len(bird_d)} D genes across {bird_d.Haplotype.nunique()} haplotypes')
    print('Loading human genes...')
    human = load_human(args.human_dir)

    vv = compute_vv_distances(bird_v, human)
    dens = compute_density(bird_v, human)
    vd = compute_vd_gap(bird_v, bird_d, human)
    cys = compute_cysteine(bird_v, bird_d, human)

    vv.to_csv(os.path.join(args.outdir, 'vv_distances.csv'), index=False)
    dens.to_csv(os.path.join(args.outdir, 'vgene_density.csv'), index=False)
    vd.to_csv(os.path.join(args.outdir, 'vd_gap.csv'), index=False)
    cys.to_csv(os.path.join(args.outdir, 'cysteine.csv'), index=False)

    plot_vv(vv, args.figdir)
    plot_density(dens, args.figdir)
    plot_vd(vd, args.figdir)
    plot_cysteine(cys, args.figdir)

    summarize(vv, dens, vd, cys)


if __name__ == '__main__':
    main()
