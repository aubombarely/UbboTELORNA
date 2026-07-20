# UbboTELORNA — Green Computing Benchmark

This page documents a focused, fast companion to the main
[correctness/performance benchmark](benchmarking.md): measuring **real
energy consumption and carbon emissions**, not just wall-clock time, for
UbboTELORNA against the tools it replaces.

---

## Motivation

The main benchmark already shows UbboTELORNA is the fastest tool at every
thread count tested against barrnap and tRNAscan-SE, and reaches near-
minimum wall time at just 4 threads (see
[benchmarking.md](benchmarking.md) and the
[results summary](../benchmark/results/UbboTELORNA_benchmark_summary.md)).
Wall-clock time and peak memory are reasonable proxies for computational
efficiency, but they're still proxies — this benchmark measures the thing
"green computing" actually refers to directly: **kg CO2eq emitted per
run**, via [`codecarbon`](https://github.com/mlco2/codecarbon).

## Methodology

Every tool invocation is wrapped in its own fresh `codecarbon`
`EmissionsTracker` instance, so all tools are measured with the same
methodology — no tool's self-reported number is mixed with another's
externally observed one (UbboTELORNA's own built-in `codecarbon` support
is explicitly disabled for this benchmark via `--disable_co2_tracking`,
for that reason).

UbboTELORNA is run once per genome as the **full pipeline** (Modules 0-4:
telomere detection, masking, rRNA annotation, tRNA annotation,
integration), and compared against the **sum** of barrnap + tRNAscan-SE +
tidk run separately for the same genome — the same "full pipeline" framing
already established in the main benchmark ("UbboTELORNA (Modules 0-4)"
vs. "barrnap + tRNAscan-SE"), extended here to include tidk so the
telomere-detection cost is counted on both sides.

## Genomes

Five genomes already used in the main 40-genome benchmark, chosen to span
~3 orders of magnitude in size while keeping the whole benchmark fast to
rerun:

| Genome key | Organism | Size |
|---|---|---|
| `ecoli_k12` | *Escherichia coli* K-12 MG1655 | 4.6 Mb |
| `scerevisiae` | *Saccharomyces cerevisiae* S288C | 12 Mb |
| `celegans` | *Caenorhabditis elegans* WS285 | 100 Mb |
| `athaliana` | *Arabidopsis thaliana* TAIR10.1 | 135 Mb |
| `osativa` | *Oryza sativa* ssp. japonica IRGSP-1.0 | 375 Mb |

## Running it

```bash
conda activate ubbotelorna_bench

# Download the 5 genomes if not already present (or use the main
# Snakefile's `download_all` target to fetch the full 40-genome set,
# which includes these five):
python3 scripts/download_genomes.py --accession GCF_000005845.2 \
    --outdir results/genomes/ecoli_k12 --genome ecoli_k12
# ... repeat for scerevisiae, celegans, athaliana, osativa

python3 scripts/green_benchmark.py --outdir results/green/ --threads 8
```

Full script: [`benchmark/scripts/green_benchmark.py`](../benchmark/scripts/green_benchmark.py).

## Output

- `results/green/green_benchmark_summary.tsv` — one row per genome × tool:
  `genome`, `tool`, `ok`, `elapsed_s`, `emissions_kg_co2eq`, `error`.
- `results/green/emissions/emissions.csv` — raw `codecarbon` output
  (per-invocation detail: duration, energy, CO2eq, hardware/region info).
- A per-genome "full pipeline" comparison (UbboTELORNA vs. the summed
  emissions of barrnap + tRNAscan-SE + tidk) printed to stdout at the end
  of the run.

## Results

Run 2026-07-20 on Salvia (`ubbotelorna_bench` env, Python 3.14, 8 threads).
`codecarbon` fell back to CPU constant-power mode on this hardware (no
RAPL/NVML available), so absolute kg CO2eq values are TDP-based estimates
rather than direct hardware energy measurements — the *relative* comparison
between tools, measured identically for all of them, is the reliable part
of this result, not the absolute numbers in isolation.

| Genome | UbboTELORNA (s / kg CO2eq) | barrnap+tRNAscan-SE+tidk (s / kg CO2eq) | Speedup | Emissions reduction |
|---|---|---|---|---|
| ecoli_k12 (4.6 Mb) | 9.88 / 0.000180 | 23.44 / 0.000428 | 2.37x | 2.38x |
| scerevisiae (12 Mb) | 15.98 / 0.000285 | 79.81 / 0.001436 | 4.99x | 5.04x |
| athaliana (135 Mb) | 60.04 / 0.001278 | 182.60 / 0.003330 | 3.04x | 2.61x |
| celegans (100 Mb) | 51.70 / 0.000929 | 195.00 / 0.003539 | 3.77x | 3.81x |
| osativa (375 Mb) | 167.49 / 0.002967 | 357.43 / 0.006438 | 2.13x | 2.17x |
| **Total (all 5)** | **305.09 / 0.005639** | **838.28 / 0.015170** | **2.75x** | **2.69x** |

(`ecoli_k12`'s tidk step is excluded — circular chromosome, no
`telomere_repeat` configured, correctly skipped rather than run.)

**Interpretation:** the wall-clock speed advantage already shown in the
main benchmark (Section 6, [summary](../benchmark/results/UbboTELORNA_benchmark_summary.md))
translates directly into a real emissions advantage of comparable
magnitude — UbboTELORNA is never slower or higher-emitting than the
combined barrnap+tRNAscan-SE+tidk pipeline for any genome tested, and the
advantage holds across the full size range (4.6 Mb to 375 Mb), not just at
one end. The largest advantage (~5x) is on the smallest eukaryotic genome
tested (*S. cerevisiae*); the smallest advantage (~2.1-2.4x) is at both
size extremes (*E. coli* and *O. sativa*) — plausibly because per-process
overhead (three separate tool startups vs. one combined pipeline) matters
proportionally more on very small genomes, while on the largest genome the
absolute compute cost of all four tools starts to converge. Worth
re-testing on a larger genome set (or the full 40) before treating the
~2-5x range as a general claim rather than a result specific to these five
genomes.
