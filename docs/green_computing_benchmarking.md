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

*Not yet run.* The script has been written and smoke-tested (argument
handling, help text, graceful failure when tools are missing), but a real
run requires `barrnap`/`tRNAscan-SE`/`tidk`/`codecarbon`, which weren't
available in the environment it was developed in. Needs to run on a
machine with the `ubbotelorna_bench` environment (e.g. Salvia) before any
numbers can be reported here.

<!--
Once run, replace this section with the actual results, e.g.:

| Genome | UbboTELORNA (kg CO2eq) | barrnap + tRNAscan-SE + tidk (kg CO2eq) | Delta |
|---|---|---|---|
| ecoli_k12 | | | |
| scerevisiae | | | |
| celegans | | | |
| athaliana | | | |
| osativa | | | |

plus a short interpretation paragraph (does the wall-clock speed
advantage shown in the main benchmark translate into a real emissions
advantage? does it hold across the whole size range, or only for small/
large genomes?).
-->
