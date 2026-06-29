# UbboTELORNA benchmark

This directory contains the reproducible benchmark comparing UbboTELORNA
against barrnap, tRNAscan-SE, and tidk across 14 genomes spanning all
major branches of life. The plant set includes *Nicotiana tabacum*
(GCF_000715075.1, allotetraploid ~4.5 Gb) as the large-genome stress test
and *Citrus sinensis* (GCF_022201045.2, ~360 Mb) for citrus diversity.

See [docs/benchmarking.md](../docs/benchmarking.md) for the full benchmark
design, genome descriptions, and rationale.

---

## Requirements

```bash
conda env create -f envs/benchmark.yaml
conda activate ubbotelorna_bench
```

The benchmark environment (`ubbotelorna_bench`) installs: `snakemake`,
`ncbi-genome-download`, `barrnap`, `trnascan-se`, `tidk`, `matplotlib`,
`pandas`, `scipy`, and `seaborn`.

UbboTELORNA itself runs inside its own environment (`ubbotelorna`), which
Snakemake activates automatically via the `conda:` directives in the rules.

---

## Quick start

```bash
cd benchmark/

# Download all genomes (run once; large files)
snakemake --cores 4 --use-conda download_all

# Run the full benchmark (replace 32 with your core count)
snakemake --cores 32 --use-conda all

# Run a single genome end-to-end
snakemake --cores 8 --use-conda \
    results/metrics/athaliana/ubbotelorna_nhmmer_rrna_metrics.tsv
```

The final report is written to `results/figures/benchmark_report.pdf`.

---

## Large genome note (Nicotiana tabacum)

*N. tabacum* (GCF_000715075.1) is an allotetraploid at ~4.5 Gb — the
largest genome in the benchmark. Allow extra time and memory for the
download and annotation steps. Consider running it last or on a dedicated
high-memory node.

---

## Output directory layout

```
benchmark/results/
├── genomes/
│   └── {genome}/
│       ├── {genome}.fasta          downloaded assembly
│       └── {genome}_ref.gff3       RefSeq reference annotations
├── annotations/
│   └── {genome}/
│       ├── ubbotelorna_nhmmer/     UbboTELORNA run (nhmmer)
│       ├── ubbotelorna_cmsearch/   UbboTELORNA run (cmsearch)
│       ├── barrnap/                barrnap output (normalised)
│       ├── trnascan/               tRNAscan-SE output (normalised)
│       └── tidk/                   tidk telomere TSV
├── robustness/
│   └── {genome}/
│       ├── {genome}_adversarial.fasta   telomere-prefix sequences
│       ├── {genome}_fragmented.fasta    10 kb scaffold simulation
│       └── robustness_summary.tsv       feature counts per scenario
├── metrics/
│   └── {genome}/
│       ├── ubbotelorna_nhmmer_rrna_metrics.tsv
│       ├── ubbotelorna_cmsearch_rrna_metrics.tsv
│       ├── barrnap_rrna_metrics.tsv
│       ├── ubbotelorna_trna_metrics.tsv
│       ├── trnascan_trna_metrics.tsv
│       └── telomere_comparison.tsv
├── performance/
│   └── {genome}/
│       ├── {tool}_t{N}.benchmark.tsv   Snakemake benchmark file
│       └── perf_summary.tsv            combined thread-sweep table
└── figures/
    ├── benchmark_report.pdf
    └── benchmark_report.png
```

---

## Running only selected benchmark dimensions

```bash
# Correctness only (all genomes with RefSeq GFF3)
snakemake --cores 16 --use-conda \
    $(snakemake --list | grep _rrna_metrics.tsv | tr '\n' ' ')

# Robustness only (Arabidopsis + Nicotiana)
snakemake --cores 16 --use-conda \
    results/robustness/athaliana/robustness_summary.tsv \
    results/robustness/ntabacum/robustness_summary.tsv

# Performance sweep only
snakemake --cores 32 --use-conda \
    results/performance/athaliana/perf_summary.tsv
```

---

## Script reference

| Script | Purpose |
|---|---|
| `scripts/download_genomes.py` | Download FASTA + GFF3 from NCBI via ncbi-genome-download |
| `scripts/normalise_gff3.py` | Convert barrnap / tRNAscan-SE output to UbboTELORNA name conventions |
| `scripts/create_adversarial.py` | Create telomere-prefix and fragmented assemblies |
| `scripts/compare_annotations.py` | Compute TP/FP/FN, sensitivity, precision, F1 |
| `scripts/benchmark_report.py` | Generate multi-panel PDF/PNG report |

---

## Interpreting the results

### Metrics TSV columns

| Column | Description |
|---|---|
| `genome` | Genome key from config.yaml |
| `tool` | Tool name |
| `feature` | `rRNA` or `tRNA` |
| `subtype` | Rfam model name (rRNA) or tRNA amino acid (tRNA); `ALL` = total |
| `TP` | Predicted features matched to a reference feature |
| `FP` | Predicted features with no reference match |
| `FN` | Reference features not covered by any prediction |
| `n_ref` | Total reference features |
| `n_pred` | Total predicted features |
| `sensitivity` | TP / (TP + FN) |
| `precision` | TP / (TP + FP) |
| `F1` | Harmonic mean of sensitivity and precision |

A feature is counted as TP if ≥ 80% of the shorter feature overlaps
(configurable via `overlap_threshold` in config.yaml).

### Benchmark TSV columns (Snakemake benchmark)

Snakemake writes these automatically for every `benchmark:` rule:

| Column | Description |
|---|---|
| `s` | Wall-clock seconds |
| `h:m:s` | Human-readable time |
| `max_rss` | Peak RSS memory in MB |
| `cpu_time` | Total CPU seconds across all threads |

### Robustness summary columns

| Column | Description |
|---|---|
| `scenario` | `native`, `adversarial`, or `fragmented` |
| `tool` | `ubbotelorna` or `barrnap` |
| `n_rrna` | rRNA features detected in this scenario |
| `n_trna` | tRNA features detected |
| `n_tel` | Telomere regions detected |

A robust tool shows equal `n_rrna` across all three scenarios.
A failing tool shows lower counts in the adversarial or fragmented scenario.
