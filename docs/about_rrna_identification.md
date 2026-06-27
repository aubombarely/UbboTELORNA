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

Rfam HMM profiles are downloaded from `https://rfam.org/family/{acc}/hmm`
and cached in `~/.ubbotelorna/rfam/` as `ubbotelorna_{kingdom}.hmm`.

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

nhmmer's memory footprint is **essentially constant per thread** regardless
of sequence or chromosome length — the HMM DP matrix is small and fixed.

cmsearch allocates a **CM DP matrix proportional to sequence length**,
one per thread.  On a genome with large chromosomes (100–800 Mb each)
and 48 threads, this can exceed 128 Gb even with `--rfam --mxsize 512`.

| Scenario | nhmmer | cmsearch |
|---|---|---|
| 4 threads, 500 Mb genome | ~200 Mb | ~4–8 Gb |
| 48 threads, 1 Gb genome | ~500 Mb | ~30–60 Gb |
| 48 threads, 3 Gb genome | ~500 Mb | OOM likely |

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
    --fasta        genome.fasta \
    --output       annotation_run/ \
    --search_tool  cmsearch \
    --cmsearch_mxsize 1024 \
    --skip_module0 \
    --skip_module1 \
    --skip_module3 \
    --skip_integration \
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
