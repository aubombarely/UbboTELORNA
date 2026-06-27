#!/usr/bin/env python3
"""Generate the UbboTELORNA benchmark report figures.

Produces a multi-panel PDF/PNG with three sections:
  Panel A — Correctness: F1 heatmap (tools × genomes) for rRNA and tRNA
  Panel B — Robustness: bar chart comparing native vs. adversarial vs.
             fragmented feature counts for UbboTELORNA and barrnap
  Panel C — Performance: runtime and memory scaling curves

Usage:
    python3 benchmark_report.py \\
        --metrics_dir    results/metrics \\
        --robustness_dir results/robustness \\
        --perf_dir       results/performance \\
        --output         results/figures/benchmark_report
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

matplotlib.rcParams.update({
    "font.size":        11,
    "axes.titlesize":   12,
    "axes.labelsize":   11,
    "figure.dpi":       150,
    "savefig.dpi":      150,
    "figure.facecolor": "white",
})

# colour palette consistent with UbboTELORNA
_C_UBBO   = "#4C9BE8"   # blue  — UbboTELORNA
_C_COMP   = "#E8604C"   # red   — comparator tools
_C_ACCENT = "#F5A623"   # amber — highlights
_C_GREY   = "#888888"

_TOOL_COLOURS = {
    "ubbotelorna_nhmmer":   _C_UBBO,
    "ubbotelorna_cmsearch": "#2E6FAA",
    "barrnap":              _C_COMP,
    "trnascan":             "#C0392B",
    "tidk":                 _C_GREY,
}

_TOOL_LABELS = {
    "ubbotelorna_nhmmer":   "UbboTELORNA\n(nhmmer)",
    "ubbotelorna_cmsearch": "UbboTELORNA\n(cmsearch)",
    "barrnap":              "barrnap",
    "trnascan":             "tRNAscan-SE",
}


# ── data loading ──────────────────────────────────────────────────────────────

def _load_metrics(metrics_dir: Path, feature: str) -> dict:
    """Return {(genome, tool): {subtype: {sensitivity, precision, F1}}}."""
    data = {}
    suffix = f"_{feature.lower()}_metrics.tsv"
    for tsv in sorted(metrics_dir.glob(f"*/*{suffix}")):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                key = (row["genome"], row["tool"])
                if key not in data:
                    data[key] = {}
                data[key][row["subtype"]] = {
                    "sensitivity": float(row["sensitivity"]),
                    "precision":   float(row["precision"]),
                    "F1":          float(row["F1"]),
                }
    return data


def _load_robustness(robustness_dir: Path) -> list:
    rows = []
    for tsv in sorted(robustness_dir.glob("*/robustness_summary.tsv")):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(reader))
    return rows


def _load_performance(perf_dir: Path) -> list:
    rows = []
    for tsv in sorted(perf_dir.glob("*/perf_summary.tsv")):
        with open(tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(reader))
    return rows


# ── Panel A: F1 heatmap ───────────────────────────────────────────────────────

def _f1_heatmap(ax, metrics: dict, feature: str, title: str) -> None:
    tools   = ["ubbotelorna_nhmmer", "ubbotelorna_cmsearch", "barrnap"] \
              if feature == "rRNA" else ["ubbotelorna_nhmmer", "trnascan"]
    genomes = sorted({g for g, _ in metrics})

    matrix = np.full((len(tools), len(genomes)), np.nan)
    for j, genome in enumerate(genomes):
        for i, tool in enumerate(tools):
            key = (genome, tool)
            if key in metrics and "ALL" in metrics[key]:
                matrix[i, j] = metrics[key]["ALL"]["F1"]

    im = ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)

    ax.set_xticks(range(len(genomes)))
    ax.set_xticklabels(genomes, rotation=40, ha="right", fontsize=9)
    ax.set_yticks(range(len(tools)))
    ax.set_yticklabels([_TOOL_LABELS.get(t, t) for t in tools], fontsize=9)
    ax.set_title(title)

    # annotate cells
    for i in range(len(tools)):
        for j in range(len(genomes)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if v > 0.7 else "black")

    plt.colorbar(im, ax=ax, shrink=0.7, label="F1 score")


# ── Panel B: robustness bar chart ─────────────────────────────────────────────

def _robustness_bars(ax, rows: list, genome: str) -> None:
    scenarios = ["native", "adversarial", "fragmented"]
    tools     = ["ubbotelorna", "barrnap"]
    width     = 0.35
    x         = np.arange(len(scenarios))

    for idx, tool in enumerate(tools):
        counts = []
        for sc in scenarios:
            match = [r for r in rows
                     if r["genome"] == genome
                     and r["scenario"] == sc
                     and r["tool"] == tool]
            counts.append(int(match[0]["n_rrna"]) if match else 0)
        colour = _C_UBBO if tool == "ubbotelorna" else _C_COMP
        label  = "UbboTELORNA" if tool == "ubbotelorna" else "barrnap"
        ax.bar(x + idx * width, counts, width, label=label, color=colour)

    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(scenarios)
    ax.set_ylabel("rRNA features detected")
    ax.set_title(f"Robustness — {genome}")
    ax.legend(fontsize=9)


# ── Panel C: performance curves ───────────────────────────────────────────────

def _perf_runtime(ax, rows: list) -> None:
    """Runtime vs thread count for one representative genome."""
    tools   = ["ubbotelorna_nhmmer", "barrnap", "trnascan"]
    genomes = sorted({r["genome"] for r in rows})
    genome  = genomes[len(genomes) // 2]  # pick median-sized genome

    for tool in tools:
        pts = sorted(
            [(int(r["threads"]), float(r["s"]) / 60)
             for r in rows
             if r["genome"] == genome and r["tool"] == tool and r["s"]],
            key=lambda x: x[0]
        )
        if not pts:
            continue
        xs, ys = zip(*pts)
        c = _TOOL_COLOURS.get(tool, _C_GREY)
        ax.plot(xs, ys, "o-", color=c, label=_TOOL_LABELS.get(tool, tool), linewidth=1.5)

    ax.set_xlabel("Threads")
    ax.set_ylabel("Wall-clock time (min)")
    ax.set_title(f"Thread scaling — {genome}")
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.legend(fontsize=9)


def _perf_memory(ax, rows: list) -> None:
    """Peak RSS vs genome size across all genomes at default threads."""
    tools = ["ubbotelorna_nhmmer", "barrnap"]

    # sizes in MB — approximate from config; read from data
    genome_order = sorted({r["genome"] for r in rows})

    for tool in tools:
        pts = []
        for i, genome in enumerate(genome_order):
            match = [r for r in rows
                     if r["genome"] == genome and r["tool"] == tool and r["max_rss"]]
            if match:
                rss_mb = float(match[0]["max_rss"])
                pts.append((i, rss_mb / 1024))  # MB → GB
        if pts:
            xs, ys = zip(*pts)
            c = _TOOL_COLOURS.get(tool, _C_GREY)
            ax.plot(xs, ys, "s-", color=c,
                    label=_TOOL_LABELS.get(tool, tool), linewidth=1.5)

    ax.set_xticks(range(len(genome_order)))
    ax.set_xticklabels(genome_order, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Peak RSS (GB)")
    ax.set_title("Memory usage vs. genome")
    ax.legend(fontsize=9)


# ── main figure ───────────────────────────────────────────────────────────────

def build_report(metrics_dir: Path, robustness_dir: Path,
                 perf_dir: Path, out_base: Path,
                 formats: list) -> None:

    rrna_metrics = _load_metrics(metrics_dir, "rRNA")
    trna_metrics = _load_metrics(metrics_dir, "tRNA")
    rob_rows     = _load_robustness(robustness_dir)
    perf_rows    = _load_performance(perf_dir)

    fig = plt.figure(figsize=(18, 22))
    # layout: 3 rows × 2 columns
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.35,
                  height_ratios=[1.2, 0.9, 0.9])

    # row 0: F1 heatmaps
    ax_rrna = fig.add_subplot(gs[0, 0])
    ax_trna = fig.add_subplot(gs[0, 1])
    if rrna_metrics:
        _f1_heatmap(ax_rrna, rrna_metrics, "rRNA", "A  rRNA — F1 score (vs. RefSeq)")
    else:
        ax_rrna.text(0.5, 0.5, "No rRNA metrics available",
                     ha="center", va="center", transform=ax_rrna.transAxes)
        ax_rrna.set_title("A  rRNA — F1 score (vs. RefSeq)")

    if trna_metrics:
        _f1_heatmap(ax_trna, trna_metrics, "tRNA", "B  tRNA — F1 score (vs. RefSeq)")
    else:
        ax_trna.text(0.5, 0.5, "No tRNA metrics available",
                     ha="center", va="center", transform=ax_trna.transAxes)
        ax_trna.set_title("B  tRNA — F1 score (vs. RefSeq)")

    # row 1: robustness (one panel per adversarial genome, side by side)
    adv_genomes = sorted({r["genome"] for r in rob_rows})
    ax_rob1 = fig.add_subplot(gs[1, 0])
    ax_rob2 = fig.add_subplot(gs[1, 1])
    rob_axes = [ax_rob1, ax_rob2]
    for i, genome in enumerate(adv_genomes[:2]):
        _robustness_bars(rob_axes[i], rob_rows, genome)
        rob_axes[i].set_title(
            f"{'C' if i == 0 else 'D'}  Robustness — {genome}")

    # row 2: performance
    ax_rt  = fig.add_subplot(gs[2, 0])
    ax_mem = fig.add_subplot(gs[2, 1])
    if perf_rows:
        _perf_runtime(ax_rt, perf_rows)
        ax_rt.set_title("E  Runtime scaling (thread sweep)")
        _perf_memory(ax_mem, perf_rows)
        ax_mem.set_title("F  Peak memory vs. genome size")
    else:
        for ax, lbl in [(ax_rt, "E"), (ax_mem, "F")]:
            ax.text(0.5, 0.5, "No performance data available",
                    ha="center", va="center", transform=ax.transAxes)

    fig.suptitle("UbboTELORNA benchmark", fontsize=14, fontweight="bold", y=1.01)

    out_base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        out_path = out_base.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"[report] saved {out_path}", file=sys.stderr)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate benchmark report figures")
    ap.add_argument("--metrics_dir",    required=True, type=Path)
    ap.add_argument("--robustness_dir", required=True, type=Path)
    ap.add_argument("--perf_dir",       required=True, type=Path)
    ap.add_argument("--output",         required=True, type=Path,
                    help="Output base path (extension added per format)")
    ap.add_argument("--format", default="pdf,png",
                    help="Comma-separated formats: pdf, png, svg (default: pdf,png)")
    args = ap.parse_args()

    formats = [f.strip().lstrip(".") for f in args.format.split(",")]
    build_report(args.metrics_dir, args.robustness_dir,
                 args.perf_dir, args.output, formats)


if __name__ == "__main__":
    main()
