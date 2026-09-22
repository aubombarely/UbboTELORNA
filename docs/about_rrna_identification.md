# About rRNA identification in UbboTELORNA

UbboTELORNA supports two search tools for rRNA annotation (Module 2):
**nhmmer** (default) and **cmsearch**.  Both use profiles from the
[Rfam database](https://rfam.org/) and produce equivalent GFF3 output,
but they differ substantially in speed, memory usage, and sensitivity.

---

## How each tool works

### nhmmer (HMMER3)

nhmmer searches a genome with **profile Hidden Markov Models (HMMs)**.
A profile HMM captures positional conservation along the rRNA sequence
but does not model RNA secondary structure.  The search algorithm is
O(L × M) where L is the genome length and M is the model length —
essentially a single pass through the sequence per model.

Rfam seed alignments (Stockholm format) are downloaded from
`https://rfam.org/family/{acc}/alignment?format=stockholm` and used to
build HMMER3 HMM profiles with `hmmbuild --rna`.  Built profiles are cached
in `~/.ubbotelorna/rfam/` as `{acc}.hmm` and concatenated into
`ubbotelorna_{kingdom}.hmm`.  The Rfam REST API does not expose standalone
HMM files, so this build step is required on first use.

### cmsearch (Infernal)

cmsearch uses **Covariance Models (CMs)**, which extend profile HMMs to
capture RNA secondary structure (base-pair covariation).  The CM search
algorithm is O(L × W²) where W is the CM state space — much more
computationally intensive, but structurally aware.

By default UbboTELORNA runs cmsearch with `--rfam` (activates the Rfam
HMM pre-filter that identifies candidate windows before running the full
CM) and `--mxsize` (caps the DP matrix per thread to avoid OOM).  Even
with these flags, cmsearch is substantially slower and more
memory-intensive than nhmmer.

Rfam CM profiles are downloaded from `https://rfam.org/family/{acc}/cm`
and cached as `ubbotelorna_{kingdom}.cm`.

---

## Speed comparison

nhmmer is typically **5–20× faster** than cmsearch on the same genome.
The table below shows approximate wall-clock times at 48 threads for
plant genomes of different sizes (estimates; actual times depend on
sequence complexity and hardware):

| Genome size | nhmmer | cmsearch (`--rfam --mxsize 512`) |
|---|---|---|
| 500 Mb | ~5–10 min | ~30–60 min |
| 1 Gb | ~10–20 min | ~60–120 min |
| 3 Gb | ~30–60 min | ~3–6 h |

---

## Memory comparison

This is the most important practical difference for large plant genomes.

nhmmer's memory footprint scales with **thread count, not sequence or
chromosome length** — but earlier estimates in this doc understated how
much: it is not a small fixed cost per thread, it is a substantial one.
Real profiling (see below) shows each thread adds roughly **270 Mb**, so
high-thread-count runs on HPC nodes (e.g. 48 threads) can reach into the
tens of Gb even though the genome itself barely matters.

cmsearch allocates a **CM DP matrix proportional to sequence length**,
one per thread, on top of its own per-thread base cost. On a genome with
large chromosomes (100–800 Mb each) and 48 threads, this can exceed
128 Gb even with `--rfam --mxsize 512`.

| Scenario | nhmmer | cmsearch |
|---|---|---|
| 4 threads, any genome size | ~1.1 Gb (measured) | ~4–8 Gb (estimate) |
| 8 threads, any genome size | ~2.2 Gb (measured) | — |
| 16 threads, any genome size | ~4.3 Gb (measured) | — |
| 24 threads, any genome size | ~6.6 Gb (measured) | — |
| 48 threads, 1 Gb genome | ~13 Gb (extrapolated) | ~30–60 Gb (estimate) |
| 48 threads, 3 Gb genome | ~13 Gb (extrapolated) | OOM likely |

### nhmmer real-world memory profiling (2026-09-22)

Measured with `/usr/bin/time -v nhmmer` on Salvia (real HPC node), single
combined Rfam euka HMM (`ubbotelorna_euka.hmm`), against the hard-masked
*Zea mays* genome (`zeaMays.hardmasked.fasta`, ~2.1 Gb):

| `--cpu` | Peak RSS | Wall clock |
|---|---|---|
| 4  | 1.13 Gb | 6:48 |
| 8  | 2.17 Gb | 3:30 |
| 16 | 4.29 Gb | 2:05 |
| 24 | 6.58 Gb | 1:27 |

Fit (least-squares across all four points, each prediction within ~0.3%
of measured): **peak RSS ≈ 41 Mb + ~272 Mb × `--cpu`**. Total CPU-seconds
is roughly constant across thread counts (~1700–1830s), so the tradeoff
is a clean one: lower `--cpu` buys lower, more predictable memory at the
cost of wall-clock time, not a change in total compute.

**Practical implication:** when annotating many large genomes in
parallel (e.g. batch mode across 100+ genomes), the nhmmer step's memory
budget should be planned from `--cpu` (via the formula above), not from
genome size — a smaller genome at 24 threads uses the same memory as a
larger one at 24 threads. Capping the nhmmer invocation's own thread
count independently of the pipeline's overall `--threads` is the
most direct lever for keeping memory bounded and predictable; a
windowed/chunked execution mode is also being considered as a
complementary fix. Not yet implemented as of this writing — track
progress in the CHANGELOG.

---

## Sensitivity comparison

cmsearch is more sensitive because it models RNA secondary structure via
base-pair covariation.  However, for **ribosomal RNA genes** — which are
among the most conserved sequences in biology — the practical difference
is small in most cases:

| rRNA gene | Length | nhmmer vs cmsearch |
|---|---|---|
| 5S rRNA | ~121 bp | Nearly identical — extremely conserved |
| 5.8S rRNA | ~160 bp | Nearly identical |
| SSU (18S) | ~1800 bp | cmsearch slightly better on divergent copies |
| LSU (28S) | ~3400 bp | cmsearch noticeably better on highly divergent copies |

The sensitivity gap is most relevant for:

- **Highly divergent lineages** (e.g., early-diverging plants, algae,
  fungi) where rRNA sequences deviate substantially from the Rfam consensus
- **Fragmented assemblies** where partial rRNA copies may be missed by
  the simpler HMM model
- **Detection of pseudogenes** or highly degenerate copies

For well-assembled, standard plant genomes (angiosperms, gymnosperms),
nhmmer finds the same copies as cmsearch in the vast majority of cases.

---

## Which tool to use

| Situation | Recommendation |
|---|---|
| Routine annotation of plant genomes | `nhmmer` (default) |
| >100 genomes in HPC batch mode | `nhmmer` (default) |
| Memory-limited nodes | `nhmmer` (default) |
| Highly divergent or early-diverging lineages | `--search_tool cmsearch` |
| Verification of suspiciously low rRNA counts | Rerun Module 2 with `--search_tool cmsearch` |
| Maximum sensitivity needed for publication | `--search_tool cmsearch` |

---

## Running cmsearch selectively

Because UbboTELORNA checkpoints each module independently, you can rerun
only Module 2 with cmsearch on a genome that has already been processed
with nhmmer — masking (Module 1) is reused:

```bash
python3 scripts/UbboTELORNA.py \
    --fasta           genome.fasta \
    --output          annotation_run/ \
    --search_tool     cmsearch \
    --cmsearch_mxsize 1024 \
    --skip_module     0,1,3,4,5 \
    --force
```

This overwrites only `results/mod02_rRNA_{prefix}.gff3` and
`workdir/cmsearch_rRNA.tblout`.

---

## References

| Tool | Reference |
|---|---|
| **HMMER3 / nhmmer** | Eddy SR (2011) *PLoS Comput Biol* 7:e1002195 |
| **Infernal / cmsearch** | Nawrocki EP, Eddy SR (2013) *Bioinformatics* 29:2933–2935 |
| **Rfam** | Kalvari I et al. (2021) *Nucleic Acids Res* 49:D192–D200 |
