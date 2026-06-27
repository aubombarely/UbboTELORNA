# UbboTELORNA — Benchmarking Plan

This document describes the planned benchmark of UbboTELORNA against
comparable tools, covering correctness, robustness, and performance across
a diverse set of genomes from prokaryotes to large polyploid plants.

---

## Tools compared

Each UbboTELORNA module is benchmarked against the tools most commonly used
for the same task.  The full-pipeline comparison (last row) reflects the
workflow that UbboTELORNA is designed to replace.

| Module | UbboTELORNA | Comparators |
|---|---|---|
| 0 — Telomere detection | Custom k-mer scan | tidk, quarTeT |
| 2 — rRNA annotation | nhmmer + Rfam HMMs (default) | barrnap, RNAmmer, cmsearch + Rfam |
| 3 — tRNA annotation | ARAGORN | tRNAscan-SE |
| Full pipeline | UbboTELORNA (Modules 0–4) | barrnap + tRNAscan-SE |

---

## Benchmark genome set

Twelve genomes are selected to cover the major branches of life, a range
of genome sizes (12 Mb – 3.1 Gb), and assembly qualities from complete
chromosomal references to highly fragmented polyploid assemblies.

### Prokaryotes

| # | Organism | Assembly | Size | rRNA operons | Source |
|---|---|---|---|---|---|
| 1 | *Escherichia coli* K-12 MG1655 | GCF_000005845.2 | 4.6 Mb | 7 (exact) | RefSeq |
| 2 | *Bacillus subtilis* 168 | GCF_000009045.1 | 4.2 Mb | 10 (exact) | RefSeq |

Both genomes are complete single-chromosome sequences with manually curated
RefSeq annotations.  The exact rRNA operon count is known from decades of
experimental work, making them the most reliable ground truth in the
benchmark.

### Fungi

| # | Organism | Assembly | Size | rDNA copies | Source |
|---|---|---|---|---|---|
| 3 | *Saccharomyces cerevisiae* S288C | GCF_000146045.2 | 12 Mb | ~150 (chr XII) | SGD / RefSeq |

*S. cerevisiae* harbours ~150 tandem rDNA copies in a single array on
chromosome XII and 274 nuclear tRNA genes curated to single-copy resolution
in the Saccharomyces Genome Database (SGD).  This makes it the ideal
organism for validating Module 6c (tandem array detection) against a known
array structure.

### Invertebrates

| # | Organism | Assembly | Size | Annotations | Source |
|---|---|---|---|---|---|
| 4 | *Caenorhabditis elegans* WS285 | GCF_000002985.6 | 100 Mb | WormBase | RefSeq |
| 5 | *Drosophila melanogaster* dm6 | GCF_000001215.4 | 140 Mb | FlyBase | RefSeq |

*D. melanogaster* is a particularly interesting test case: rDNA is located
in the heterochromatic regions of the X and Y chromosomes, which are often
collapsed or incomplete in genome assemblies.  This tests how tools handle
fragmented or absent rDNA loci.

### Plants

| # | Organism | Assembly | Size | Ploidy | Annotations | Source |
|---|---|---|---|---|---|---|
| 6 | *Selaginella moellendorffii* v1.0 | GCF_000143415.1 | 100 Mb | Diploid | Phytozome | RefSeq |
| 7 | *Arabidopsis thaliana* TAIR10.1 | GCF_000001735.4 | 135 Mb | Diploid | TAIR | RefSeq |
| 8 | *Oryza sativa* ssp. japonica IRGSP-1.0 | GCF_001433935.1 | 375 Mb | Diploid | RAP-DB | RefSeq |
| 9 | *Actinidia arguta* | — | ~700 Mb | 4× polyploid | — | — |
| 10 | *Zea mays* B73 RefGen_v4 | GCF_000005005.2 | 2.1 Gb | Diploid | MaizeGDB | RefSeq |

*Selaginella* is a lycophyte — an early-diverging land plant outside the
seed plant lineage.  Its rRNA sequences are more divergent from the Rfam
consensus than angiosperm sequences, testing whether nhmmer and cmsearch
maintain sensitivity on distant homologues.

*Actinidia arguta* is the validation genome used during tool development
(~700 Mb, 4× polyploid, 3,961 sequences; 71 telomeres, 21,871 rRNA copies,
3,546 tRNA genes confirmed).

*Zea mays* (~85% repetitive) is the largest and most repetitive plant
genome in the set and is used primarily as a performance stress test.

### Vertebrates

| # | Organism | Assembly | Size | Annotations | Source |
|---|---|---|---|---|---|
| 11 | *Gallus gallus* GRCg7b | GCF_016699485.2 | 1.0 Gb | Ensembl | RefSeq |
| 12 | *Danio rerio* GRCz11 | GCF_000002035.6 | 1.4 Gb | ZFIN / Ensembl | RefSeq |
| 13 | *Homo sapiens* T2T-CHM13v2.0 | GCF_009914755.1 | 3.1 Gb | RefSeq / T2T | RefSeq |

The human T2T-CHM13 assembly is preferred over GRCh38 because it resolves
the five acrocentric rDNA arrays (NORs on chromosomes 13, 14, 15, 21, 22)
that are collapsed or missing in GRCh38.  This makes it directly applicable
to Module 6c (array detection) and gives a biologically grounded expected
copy number.

---

## Reference annotations

Ground-truth annotations are drawn from three sources:

### RefSeq GFF3

RefSeq provides manually curated rRNA and tRNA feature annotations for
most assemblies listed above.  These are used as the primary reference for
sensitivity and precision calculations.

### GtRNAdb (Genomic tRNA Database)

[GtRNAdb](http://gtrnadb.ucsc.edu/) provides tRNAscan-SE predictions for
all 13 genomes in this set, computed with consistent parameters.  Using
GtRNAdb as the tRNA reference provides a direct, consistent comparison
between UbboTELORNA (ARAGORN) and tRNAscan-SE across all organisms.

### SILVA and Rfam

SILVA (release 138.2) provides curated rRNA sequences for all kingdoms.
Rfam (release 14.10) CM-based annotations serve as an independent rRNA
reference, particularly for lineages with divergent rRNA sequences
(Selaginella, Drosophila rDNA heterochromatin).

---

## Benchmark dimensions

### 1. Correctness

Annotations produced by each tool are compared against the reference GFF3
using reciprocal overlap (≥ 80% of the shorter feature must overlap).

Metrics computed per tool, per genome, per feature class:

| Metric | Formula |
|---|---|
| Sensitivity (recall) | TP / (TP + FN) |
| Precision | TP / (TP + FP) |
| F1 score | 2 × (precision × recall) / (precision + recall) |
| Copy count | n detected vs. n in reference |

### 2. Robustness — the key differentiator

UbboTELORNA was created because `barrnap` (and tools that call it
internally) fails on sequences that begin with telomeric repeats or other
low-complexity patterns.  The robustness benchmark tests this directly:

1. Run all tools on the native assemblies.
2. Create adversarial versions of Arabidopsis and Actinidia by prepending
   1,000 copies of the telomeric repeat unit (`TTTAGGG`) to each
   chromosome sequence.
3. Compare rRNA hit counts, particularly on the first and last 500 kb of
   each chromosome, between native and adversarial assemblies.

Expected result: barrnap misses or fails on sequences with terminal
repeats; UbboTELORNA (with Module 1 masking) is unaffected.

Additional robustness tests:

- **Fragmented assembly simulation** — split each genome into 10 kb
  scaffolds and compare hit counts vs. the intact assembly.
- **Thread scaling** — run at 1, 4, 8, 16, 32, 48 threads and verify
  output is identical (not just similar) across thread counts.

### 3. Performance

Wall-clock time and peak RSS memory are recorded for every tool on every
genome.  Scaling curves are plotted for:

- Runtime vs. genome size (at fixed thread count)
- Runtime vs. thread count (at fixed genome size)
- Peak memory vs. genome size

---

## What each lineage adds

| Lineage | Primary contribution |
|---|---|
| Prokaryotes | Exact ground truth; `--kingdom bacteria` validation |
| Fungi | Single-chromosome tandem rDNA array; Module 6c validation |
| Invertebrates | Heterochromatic rDNA (Drosophila); compact well-annotated genomes |
| Selaginella | Divergent rRNA sequences; sensitivity on non-angiosperm plants |
| Arabidopsis / rice | Plant gold standards with curated annotations |
| Actinidia | Polyploid challenge; real validation dataset from tool development |
| Maize | Large repetitive genome; performance stress test |
| Chicken / zebrafish | Mid-size vertebrates with good reference annotations |
| Human T2T | Resolved rDNA arrays; reviewer-facing reference; largest genome |

---

## Expected output table

The benchmark will produce a single comparison table of the form:

```
Tool              Genome         Module  Sensitivity  Precision    F1   Time(min)  Mem(Gb)
UbboTELORNA(nhm)  Arabidopsis    rRNA       0.99        0.98      0.98      2.1       0.8
barrnap           Arabidopsis    rRNA       0.95        0.99      0.97      1.8       0.6
RNAmmer           Arabidopsis    rRNA       0.91        0.97      0.94     12.3       2.1
cmsearch+Rfam     Arabidopsis    rRNA       0.99        0.99      0.99      8.4       4.2
UbboTELORNA(cms)  Arabidopsis    rRNA       0.99        0.99      0.99      9.1       4.5
UbboTELORNA       Arabidopsis    tRNA       0.97        0.98      0.97      0.4       0.3
tRNAscan-SE       Arabidopsis    tRNA       0.99        0.99      0.99      1.1       0.5
...
```

---

## Planned implementation

The benchmark will live in a dedicated `benchmark/` directory:

```
benchmark/
├── config.yaml                  genomes, tools, parameters, thread counts
├── Snakefile                    orchestrates downloads, runs, and comparisons
├── envs/
│   └── benchmark.yaml           conda environment with all comparator tools
├── scripts/
│   ├── download_genomes.py      fetch assemblies and reference GFF3 from NCBI
│   ├── compare_annotations.py   compute TP/FP/FN vs. reference; write metrics TSV
│   └── benchmark_report.py      generate comparison figures (matplotlib)
└── README.md                    quickstart for reproducing the benchmark
```

The `compare_annotations.py` script will accept any two GFF3 files
(predicted vs. reference) and produce a metrics TSV, making it reusable
outside this benchmark.

---

## Citation of comparator tools

| Tool | Reference |
|---|---|
| **barrnap** | Seemann T (2013) barrnap — BAsic Rapid Ribosomal RNA Predictor. GitHub |
| **RNAmmer** | Lagesen K et al. (2007) *Nucleic Acids Res* 35:3100–3108 |
| **tRNAscan-SE** | Chan PP et al. (2021) *Nucleic Acids Res* 49:D421–D430 |
| **ARAGORN** | Laslett D, Canback B (2004) *Nucleic Acids Res* 32:11–16 |
| **tidk** | Hall MB (2022) tidk — Telomere Identification Toolkit. GitHub |
| **quarTeT** | Li H et al. (2022) *Bioinformatics* 38:5511–5513 |
| **HMMER3 / nhmmer** | Eddy SR (2011) *PLoS Comput Biol* 7:e1002195 |
| **Infernal / cmsearch** | Nawrocki EP, Eddy SR (2013) *Bioinformatics* 29:2933–2935 |
| **Rfam** | Kalvari I et al. (2021) *Nucleic Acids Res* 49:D192–D200 |
| **GtRNAdb** | Chan PP, Lowe TM (2016) *Nucleic Acids Res* 44:D184–D189 |
| **SILVA** | Quast C et al. (2013) *Nucleic Acids Res* 41:D590–D596 |
