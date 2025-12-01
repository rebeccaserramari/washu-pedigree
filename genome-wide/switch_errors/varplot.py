from __future__ import annotations
import math
import multiprocessing as mp
from dataclasses import dataclass
from functools import partial
from typing import Dict, Iterable, List, Optional, Tuple
import bisect
import itertools

import numpy as np
import pandas as pd
from matplotlib import gridspec
from tqdm import tqdm

import matplotlib.pyplot as plt
import seaborn as sns


def _parse_contig(contig: str) -> Tuple[str, str]:
    """
    Parse contig string of form PAN027.{chrom}.{hap} into (chrom, hap)
    """


def _normalize_kind(kind: Optional[str]) -> Optional[str]:
    """
    Normalize switch_error_kind to 'mother' or 'grandparent' or None
    """

    if kind is None:
        return None
    k = kind.lower()
    if 'mother' in k:
        return 'mother'
    if 'grandparent' in k or 'grand parent' in k or 'grand_parent' in k:
        return 'grandparent'
    if 'gp' in k or 'grand' in k:
        return 'grandparent'
    return k


@dataclass
class IntervalIndex:
    """
    Simple interval index per contig for fast membership queries.
    Builds sorted lists of (starts, ends) and we check membership via bisect.
    """
    starts: List[int]
    ends: List[int]

    def contains(self, pos: int) -> bool:
        # find rightmost interval with start <= pos
        i = bisect.bisect_right(self.starts, pos) - 1
        if i < 0:
            return False
        return self.ends[i] >= pos


def build_interval_index(clusters: Iterable[Tuple[str, int, int, int]]) -> Dict[str, IntervalIndex]:
    """
    Group clusters by contig and build IntervalIndex for each contig.
    Assumes clusters are non-overlapping or if overlapping they will be merged by
    sorted order to ensure correctness for membership queries.
    """
    grouped = {}
    for contig, start, end, count in clusters:
        grouped.setdefault(contig, []).append((int(start), int(end)))
    idx = {}
    for contig, ivs in grouped.items():
        ivs_sorted = sorted(ivs, key=lambda x: x[0])
        # optionally merge overlapping intervals for safety
        merged = []
        for s, e in ivs_sorted:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        starts = [s for s, e in merged]
        ends = [e for s, e in merged]
        idx[contig] = IntervalIndex(starts, ends)
    return idx


# ---------- Window counting logic ----------

def _sliding_windows_for_group(
    group_variants: List[Tuple[int, int, str, Optional[str]]],
    interval_index: Optional[IntervalIndex],
    window_size: int = 100_000,
    step: int = 1_000,
) -> pd.DataFrame:
    """Compute sliding windows for a single contig+haplotype group.

    group_variants: list of tuples (pos_from, pos_to, is_switch_kind, switch_kind)
      where is_switch_kind is 'none' or 'mother' or 'grandparent'

    Returns a DataFrame with columns:
      ['window_start', 'window_end', 'n_variants',
       'n_mother_in', 'n_mother_out', 'n_gp_in', 'n_gp_out']
    """
    if len(group_variants) == 0:
        return pd.DataFrame(
            columns=[
                'window_start','window_end','n_variants',
                'n_mother_in','n_mother_out','n_gp_in','n_gp_out'
            ]
        )

    positions = sorted([int(v[0]) for v in group_variants])
    min_pos = positions[0]
    max_pos = positions[-1]
    if max_pos - min_pos < 0:
        max_pos = min_pos

    pos_arr = np.array([int(v[0]) for v in group_variants], dtype=int)
    kinds = np.array([_normalize_kind(v[3]) for v in group_variants], dtype=object)
    in_cluster_arr = np.array([
        (interval_index.contains(p) if interval_index is not None else False) for p in pos_arr
    ], dtype=bool)

    starts = list(range(min_pos, max_pos - window_size + step + 1, step))
    rows = []

    for s in starts:
        e = s + window_size
        mask = (pos_arr >= s) & (pos_arr < e)
        n_variants = int(mask.sum())
        if n_variants == 0:
            rows.append((s, e, 0, 0, 0, 0, 0))
            continue
        selected_kinds = kinds[mask]
        selected_in_cluster = in_cluster_arr[mask]
        # mother
        mother_mask = selected_kinds == 'mother'
        gp_mask = selected_kinds == 'grandparent'
        n_mother_in = int((mother_mask & selected_in_cluster).sum())
        n_mother_out = int((mother_mask & ~selected_in_cluster).sum())
        n_gp_in = int((gp_mask & selected_in_cluster).sum())
        n_gp_out = int((gp_mask & ~selected_in_cluster).sum())
        rows.append((s, e, n_variants, n_mother_in, n_mother_out, n_gp_in, n_gp_out))

    df = pd.DataFrame(rows, columns=[
        'window_start','window_end','n_variants',
        'n_mother_in','n_mother_out','n_gp_in','n_gp_out'
    ])

    df['window_mid'] = ((df['window_start'] + df['window_end']) / 2).astype(int)
    return df


def _worker_compute(args):
    """
    Worker wrapper for multiprocessing.
    args = (chr_hap, group_variants, interval_index, window_size, step)
    """
    chr_hap, group_variants, interval_index, window_size, step = args
    df = _sliding_windows_for_group(group_variants, interval_index, window_size, step)
    df['chrom_hap'] = chr_hap
    return df



class VariantSwitchErrorPlotter:
    def __init__(self, variants: Dict, sw_error_prone_variants: Dict, clusters: Iterable[Tuple[str,int,int,int]]):
        """
        Initialize with raw data structures

        Parameters
        ----------
        variants                - dict of id -> dict (contig, pos_from, pos_to...)
        sw_error_prone_variants - dict subset of variants with switch_error_kind field
        clusters                - iterable of (contig, start, end, count)
        """

        self.variants = variants
        self.sw_variants = sw_error_prone_variants
        self.clusters = list(clusters)
        self.cluster_index = build_interval_index(self.clusters)

        # Precompute mapping from variant id to kind if present
        self.variant_kind = {}
        for vid, data in self.sw_variants.items():
            k = _normalize_kind(data.get('switch_error_kind'))
            self.variant_kind[vid] = k

    def _group_variants(self, chromosomes: Optional[Iterable[str]] = None, haplotypes: Optional[Iterable[str]] = None):
        """
        Yield groups by (chromosome,haplotype) as key and list of tuples

        Parameters
        ----------
        chromosomes
        haplotypes

        Returns
        -------
        (pos_from, pos_to, variant_id, switch_error_kind)
        """

        groups = {}
        for vid, data in self.variants.items():
            contig = data.get('contig')
            if contig is None:
                continue
            chrom, hap = _parse_contig(contig)
            if chromosomes is not None and chrom not in chromosomes:
                continue
            if haplotypes is not None and hap not in haplotypes:
                continue
            pos = int(data.get('pos_from') if data.get('pos_from') is not None else data.get('pos_to'))
            kind = self.variant_kind.get(vid)
            key = f"{chrom}.{hap}"
            groups.setdefault((contig, key), []).append((pos, int(data.get('pos_to', pos)), vid, kind))
        return groups

    def compute_all(self,
                    chromosomes: Optional[Iterable[str]] = None,
                    haplotypes: Optional[Iterable[str]] = None,
                    window_size: int = 100_000,
                    step: int = 1_000,
                    workers: int = 1) -> pd.DataFrame:
        """
        Compute sliding-window counts for all contig+haplotype groups
        Parameters
        ----------
        chromosomes
        haplotypes
        window_size
        step
        workers

        Returns
        -------
        concatenated DataFrame
        """

        groups = self._group_variants(chromosomes, haplotypes)
        tasks = []
        for (contig, chromhap), group_variants in groups.items():
            interval_index = self.cluster_index.get(contig)
            tasks.append((chromhap, group_variants, interval_index, window_size, step))

        results = []
        if workers and workers > 1:
            with mp.Pool(processes=workers) as pool:
                for df in tqdm(pool.imap_unordered(_worker_compute, tasks), total=len(tasks), desc='groups'):
                    results.append(df)
        else:
            for t in tqdm(tasks, desc='groups'):
                df = _worker_compute(t)
                results.append(df)

        if not results:
            return pd.DataFrame()
        combined = pd.concat(results, ignore_index=True)
        return combined

    def _prepare_scatter_dataframe(self, windows_df: pd.DataFrame) -> pd.DataFrame:
        """
        Vectorized long-form transformation.
        Creates a 4× longer DataFrame with columns:
          n_variants, n_unreliable, category
        Categories:
          - grandparent_in
          - grandparent_out
          - mother_in
          - mother_out
        """

        df = windows_df
        mapping = {
            "n_gp_in": "grandparent_in",
            "n_gp_out": "grandparent_out",
            "n_mother_in": "mother_in",
            "n_mother_out": "mother_out",
        }

        frames = []
        for col, cat in mapping.items():
            sub = df[["n_variants", col]].rename(columns={col: "n_unreliable"})
            sub["category"] = cat
            frames.append(sub)

        return pd.concat(frames, ignore_index=True)

    def plot(self, windows_df: pd.DataFrame, out_file: Optional[str] = None, figsize=(10,10)):
        """
        Create scatter + KDE using a seaborn JointGrid

        Parameters
        ----------
        windows_df  - dataframe to plot
        out_file    - output file
        figsize     - figure size

        Returns
        -------
`       None
        """

        import seaborn as sns
        import matplotlib.pyplot as plt

        if windows_df.empty:
            raise ValueError("windows_df is empty. Run compute_all() first.")

        plot_df = self._prepare_scatter_dataframe(windows_df)

        sns.set(style="white", rc={"axes.grid": False})

        palette: dict[str, str] = {
            'grandparent_in': 'blue',
            'grandparent_out': 'green',
            'mother_in': 'red',
            'mother_out': 'orange'
        }

        category_labels: dict[str, str] = {
            'grandparent_in': 'Switch error in grandparent (inside cluster)',
            'grandparent_out': 'Switch error in grandparent (outside cluster)',
            'mother_in': 'Switch error in mother (inside cluster)',
            'mother_out': 'Switch error in mother (outside cluster)'
        }

        fig = plt.figure(figsize=figsize)

        gs = gridspec.GridSpec(
            4, 4,
            height_ratios=[0.12, 1, 1, 1],
            width_ratios=[1, 1, 1, 0.12],
            hspace=0.05,
            wspace=0.05
        )

        ax_main = fig.add_subplot(gs[1:, :3])
        ax_top = fig.add_subplot(gs[0, :3], sharex=ax_main)
        ax_right = fig.add_subplot(gs[1:, 3], sharey=ax_main)

        for cat, sub in plot_df.groupby('category'):
            ax_main.scatter(
                sub['n_variants'], sub['n_unreliable'],
                s=5, alpha=0.3, color=palette[cat],
                label=category_labels[cat]
            )

        ax_main.set_xlabel("Number of variants (per 100kb)")
        ax_main.set_ylabel("Number of unreliable variants (per 100kb)")

        ax_main.legend(
            loc='center left',
            bbox_to_anchor=(0.05, 0.9),
            fontsize='small'
        )

        # ========== TOP KDE ==========
        sns.kdeplot(
            data=plot_df,
            x="n_variants",
            ax=ax_top,
            fill=True,
            bw_adjust=1.3
        )
        ax_top.set_ylabel("")
        ax_top.set_xlabel("")
        ax_top.set_yticks([])
        sns.despine(ax=ax_top, left=True, bottom=True)

        # ========== RIGHT KDE ==========
        sns.kdeplot(
            data=plot_df,
            y="n_unreliable",
            ax=ax_right,
            fill=True,
            bw_adjust=1.3
        )
        ax_right.set_xlabel("")
        ax_right.set_ylabel("")
        ax_right.set_xticks([])
        sns.despine(ax=ax_right, left=True, bottom=True)

        fig.suptitle("Variants vs Switch Error Variants (clusters) [cluster=20kb dist, min 8 SNVs]", fontsize=14)
        fig.tight_layout()

        if out_file:
            fig.savefig(out_file, dpi=1200)

        return fig
