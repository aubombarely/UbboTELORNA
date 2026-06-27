<p align="center">
<img src="assets/UbboTELORNA_logo.svg" width="260"/>
</p>

<p align="center">
<img src="https://img.shields.io/badge/version-v0.1.0-teal"/>
<img src="https://img.shields.io/badge/python-3.10%2B-blue"/>
<img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey"/>
<a href="CHANGELOG.md"><img src="https://img.shields.io/badge/changelog-v0.1.0-orange"/></a>
</p>

---

## UbboTELORNA

**UbboTELORNA** annotates three classes of ancient, conserved genomic elements
in genome assemblies: telomeres, ribosomal RNA (rRNA) genes, and transfer RNA
(tRNA) genes — with built-in low-complexity masking that prevents tool failures
on telomeric and repetitive sequences.

The name blends **Ubbo-Sathla** (Clark Ashton Smith, 1933) — the primordial
protoplasmic source of all terrestrial life in the Lovecraft Mythos — with
**TELO**mere and **RNA**.  Just as Ubbo-Sathla is the origin from which all
biological complexity descends, telomeres, rDNA, and tRNA represent the most
ancient, conserved elements of the eukaryotic genome.

### Why not barrnap?

`barrnap` (internally `nhmmer`) fails when a sequence starts with a
low-complexity pattern such as a telomeric repeat (`ACACAC…`, `TTTAGGG…`).
UbboTELORNA solves this by:

1. **Identifying telomeres first** (Module 0) — pure-Python k-mer scan, no external dependencies.
2. **Masking low-complexity regions** (Module 1) — `tantan` soft-mask before any search tool runs.
3. **Flexible rRNA search** (Module 2) — `nhmmer` (default, fast, low memory) or Infernal
   `cmsearch` (higher sensitivity) via `--search_tool`.

See [docs/about_rrna_identification.md](docs/about_rrna_identification.md) for a detailed
comparison of nhmmer vs cmsearch: speed, memory, sensitivity, and when to use each.

---

## Overview

| Module | Step | Tool | Output |
|---|---|---|---|
| 0 | Telomere identification | Custom k-mer scan | `mod00_telomeres_{prefix}.gff3` |
| 1 | Low-complexity masking | tantan | `workdir/masked_soft.fasta`, `workdir/masked_hard.fasta` |
| 2 | rRNA annotation | nhmmer (default) or cmsearch + Rfam profiles | `mod02_rRNA_{prefix}.gff3` |
| 3 | tRNA annotation | ARAGORN | `mod03_tRNA_{prefix}.gff3` |
| 4 | Integration | — | `mod04_annotation_{prefix}.gff3`, `mod04_summary_{prefix}.tsv` |

### Rfam models used

| Kingdom | Model name | Rfam accession | Target |
|---|---|---|---|
| `euka` | 5S_rRNA | RF00001 | 5S ribosomal RNA |
| `euka` | 5_8S_rRNA | RF00002 | 5.8S ribosomal RNA |
| `euka` | SSU_rRNA_eukarya | RF01960 | 18S small subunit rRNA |
| `euka` | LSU_rRNA_eukarya | RF02543 | 28S large subunit rRNA |
| `bacteria` | 5S_rRNA | RF00001 | 5S rRNA |
| `bacteria` | SSU_rRNA_bacteria | RF00177 | 16S rRNA |
| `bacteria` | LSU_rRNA_bacteria | RF02541 | 23S rRNA |
| `archaea` | 5S_rRNA | RF00001 | 5S rRNA |
| `archaea` | SSU_rRNA_archaea | RF01959 | 16S archaeal rRNA |
| `archaea` | LSU_rRNA_archaea | RF02540 | 23S archaeal rRNA |

On first use, models are fetched automatically from `https://rfam.org/`:

- **nhmmer** — downloads the Rfam seed alignment (Stockholm) for each
  family and builds a HMMER3 `.hmm` profile with `hmmbuild --rna`.
- **cmsearch** — downloads the Rfam `.cm` covariance model directly.

All files are cached in `~/.ubbotelorna/rfam/`.  On **air-gapped HPC
nodes**, build or download the files once on a machine with internet
access and supply the directory with `--rfam_dir /path/to/rfam/`.

---

## Requirements

```bash
conda env create -f envs/UbboTELORNA.yaml
conda activate ubbotelorna
```

| Dependency | Role | Install |
|---|---|---|
| Python ≥ 3.10 | Runtime | included in conda env |
| `tantan` | Low-complexity masking (Module 1) | `conda install -c bioconda tantan` |
| `hmmer` (`nhmmer`) | rRNA annotation — default (Module 2) | `conda install -c bioconda hmmer` |
| `infernal` (`cmsearch`) | rRNA annotation — alternative (Module 2) | `conda install -c bioconda infernal` |
| `aragorn` | tRNA annotation (Module 3) | `conda install -c bioconda aragorn` |
| `codecarbon` | Carbon footprint tracking (optional) | `conda install -c conda-forge codecarbon` |

---

## Quick start

```bash
python3 scripts/UbboTELORNA.py \
    --fasta  genome.fasta \
    --output annotation_run/
```

---

## Usage

```
UbboTELORNA.py --fasta FASTA --output DIR [options]
```

---

## Options

### Required

| Flag | Description |
|---|---|
| `--fasta` | Input genome assembly FASTA |
| `--output` | Output directory (created if absent) |

### Module 0 — Telomere identification

| Flag | Default | Description |
|---|---|---|
| `--telomere_repeat` | auto-detect | Known repeat unit (e.g. `TTTAGGG` plants, `TTAGGG` vertebrates) |
| `--telomere_window` | 10000 | bp to scan at each contig end |
| `--telomere_density` | 0.5 | Minimum repeat density (0–1) to call a telomere |
| `--telomere_min_len` | 100 | Minimum telomere length to report (bp) |

### Module 2 — rRNA annotation

| Flag | Default | Description |
|---|---|---|
| `--kingdom` | `euka` | Organism kingdom: `euka`, `bacteria`, or `archaea` |
| `--evalue` | 1e-5 | E-value threshold for rRNA search |
| `--rfam_dir` | `~/.ubbotelorna/rfam/` | Directory with pre-downloaded Rfam `.hmm`/`.cm` files |
| `--search_tool` | `nhmmer` | rRNA search tool: `nhmmer` (fast, low memory) or `cmsearch` (higher sensitivity) |
| `--cmsearch_mxsize` | 512 | Max DP matrix per cmsearch thread (Mb) — only used with `--search_tool cmsearch` |

See [docs/about_rrna_identification.md](docs/about_rrna_identification.md) for a full
comparison of the two tools.

### General

| Flag | Default | Description |
|---|---|---|
| `--threads` | 4 | CPU threads for the rRNA search tool |
| `--skip_module0` | — | Skip Module 0 (telomere identification) |
| `--skip_module1` | — | Skip Module 1 (low-complexity masking) |
| `--skip_module2` | — | Skip Module 2 (rRNA annotation) |
| `--skip_module3` | — | Skip Module 3 (tRNA annotation) |
| `--skip_integration` | — | Skip Module 4 (do not write merged GFF3) |
| `--force` | — | Rerun all steps even if outputs already exist |
| `--dry_run` | — | Validate inputs, print steps, exit |
| `--disable_co2_tracking` | — | Disable codecarbon carbon tracking |
| `--version` | — | Print version and exit |

---

## Output directory layout

```
{output}/
├── results/
│   ├── mod00_telomeres_{prefix}.gff3       Telomere features (Module 0)
│   ├── mod02_rRNA_{prefix}.gff3            rRNA features (Module 2)
│   ├── mod03_tRNA_{prefix}.gff3            tRNA features (Module 3)
│   ├── mod04_annotation_{prefix}.gff3      Combined GFF3 (Module 4)
│   ├── mod04_summary_{prefix}.tsv          Detailed feature summary (count, length, % genome)
│   └── {prefix}.run_summary.json           Run metadata and resource usage
├── workdir/
│   ├── masked_soft.fasta                   Soft-masked FASTA (tantan lowercase)
│   ├── masked_hard.fasta                   Hard-masked FASTA (N's, used by search tools)
│   ├── nhmmer_rRNA.tblout                  Raw nhmmer tabular output (or cmsearch_rRNA.tblout)
│   └── aragorn.txt                         Raw ARAGORN output
└── logs/
    ├── Run_UbboTELORNA.log                 Full timestamped run log
    └── {prefix}.emissions.csv              Carbon footprint (codecarbon)
```

### GFF3 attribute fields

**Telomere features:**
```
ID=tel_{seq}_{end}_{n};Name=telomere_{5prime|3prime};repeat_unit=TTTAGGG;density=0.923
```

**rRNA features:**
```
ID=rRNA_{n};Name=5S_rRNA;model=RF00001;score=142.3;E-value=1.2e-42
```

**tRNA features:**
```
ID=tRNA_{n};Name=tRNA-Phe(GAA);anticodon=GAA
```

### Summary TSV

`results/mod04_summary_{prefix}.tsv` contains one row per subtype plus a
`TOTAL` row for each feature class.  All GFF3 outputs are sorted by SeqID
(natural order, so `Chr2 < Chr10`) then by start coordinate.

```
feature_type    subtype              count    total_length_bp    pct_genome
telomere        TOTAL                71       710000             0.1014
telomere        3prime               36       360000             0.0514
telomere        5prime               35       350000             0.0500
rRNA            TOTAL                21871    43742000           6.2489
rRNA            5S_rRNA              10233    1238193            0.1769
rRNA            5_8S_rRNA            3372     539520             0.0771
rRNA            SSU_rRNA_eukarya     3177     5890476            0.8415
rRNA            LSU_rRNA_eukarya     3279     36074311           5.1534
tRNA            TOTAL                3546     265950             0.0380
tRNA            tRNA-Ala             312      23400              0.0033
tRNA            tRNA-Gly             289      21675              0.0031
...
```

- `total_length_bp` — sum of (end − start + 1) for all features of that subtype
- `pct_genome` — `total_length_bp / genome_size × 100` (4 decimal places; `NA` if genome size unavailable)
- rRNA subtypes are written in biological order (5S → 5.8S → SSU → LSU)
- tRNA subtypes are sorted by count descending, then alphabetically

---

## Run log

The full timestamped log is written to `logs/Run_UbboTELORNA.log`:

```
==============================================================
  UbboTELORNA v0.1.0  —  Run Log
==============================================================
Date      : 2026-06-27 10:14:53
User      : jdoe
Server    : hpc-node-01
OS        : Linux 5.15.0 (x86_64)
Directory : /home/jdoe/projects/apple_genome
Command   : scripts/UbboTELORNA.py --fasta genome.fasta --output annotation_run/
```

---

## Run summary JSON

`results/{prefix}.run_summary.json` is written at the end of every run:

```json
{
  "date": "2026-06-27 10:22:11",
  "version": "v0.1.0",
  "input_fasta": "/data/apple/genome.fasta",
  "genome_size_bp": 699868534,
  "kingdom": "euka",
  "telomere_repeat": "TTTAGGG",
  "parameters": {
    "telomere_window": 10000,
    "telomere_density": 0.5,
    "telomere_min_len": 100,
    "evalue": 1e-05,
    "search_tool": "nhmmer",
    "threads": 48
  },
  "feature_counts": {
    "telomere": {
      "total": 71,
      "total_length_bp": 710000,
      "by_end":     { "5prime": 35, "3prime": 36 },
      "len_by_end": { "5prime": 350000, "3prime": 360000 }
    },
    "rRNA": {
      "total": 21871,
      "total_length_bp": 43742000,
      "by_type":    { "5S_rRNA": 10233, "5_8S_rRNA": 3372, "SSU_rRNA_eukarya": 3177, "LSU_rRNA_eukarya": 3279 },
      "len_by_type": { "5S_rRNA": 1238193, "5_8S_rRNA": 539520, "SSU_rRNA_eukarya": 5890476, "LSU_rRNA_eukarya": 36074311 }
    },
    "tRNA": {
      "total": 3546,
      "total_length_bp": 265950,
      "by_type":    { "tRNA-Ala": 312, "tRNA-Gly": 289, "...": "..." },
      "len_by_type": { "tRNA-Ala": 23400, "tRNA-Gly": 21675, "...": "..." }
    }
  },
  "resource_usage": {
    "wall_clock_s": 312.4,
    "peak_mem_mb": 1842.3,
    "emissions_kg_CO2eq": 0.000021
  }
}
```

---

## Carbon footprint

When `codecarbon` is installed, emissions are tracked automatically and saved to
`logs/{prefix}.emissions.csv`.  Disable with `--disable_co2_tracking`.

---

## Recommended workflow

```bash
# 1. Activate environment
conda activate ubbotelorna

# 2. Dry run — validate inputs and preview steps
python3 scripts/UbboTELORNA.py \
    --fasta genome.fasta \
    --output annotation_run/ \
    --dry_run

# 3. Run full pipeline (eukaryote)
python3 scripts/UbboTELORNA.py \
    --fasta   genome.fasta \
    --output  annotation_run/ \
    --kingdom euka \
    --threads 8

# 4. Inspect telomere detection
grep "telomere" annotation_run/results/mod00_telomeres_*.gff3 | wc -l

# 5. Check rRNA counts
grep -v "^#" annotation_run/results/mod02_rRNA_*.gff3 \
    | cut -f9 | grep -oP 'Name=\K[^;]+' | sort | uniq -c | sort -rn

# 6. Resume from checkpoint (if run was interrupted)
python3 scripts/UbboTELORNA.py \
    --fasta  genome.fasta \
    --output annotation_run/
    # Modules already completed are skipped automatically

# 7. Force full rerun
python3 scripts/UbboTELORNA.py \
    --fasta  genome.fasta \
    --output annotation_run/ \
    --force

# 8. Supply known telomere repeat (skip auto-detection)
python3 scripts/UbboTELORNA.py \
    --fasta            genome.fasta \
    --output           annotation_run/ \
    --telomere_repeat  TTAGGG          # vertebrate

# 9. Air-gapped HPC: use local Rfam CM directory
python3 scripts/UbboTELORNA.py \
    --fasta    genome.fasta \
    --output   annotation_run/ \
    --rfam_dir /shared/databases/rfam_cms/
```

---

## Quicktest

```bash
conda activate ubbotelorna

python3 scripts/UbboTELORNA.py \
    --fasta   test/test_genome.fasta \
    --output  test_run/ \
    --kingdom euka \
    --threads 2

# Offline (no internet): skip Module 2 rRNA annotation
python3 scripts/UbboTELORNA.py \
    --fasta        test/test_genome.fasta \
    --output       test_run_offline/ \
    --skip_module2 \
    --threads 2
```

See `test/README.md` for expected outputs.

---

## Integration with YuggASMoth

UbboTELORNA is a drop-in replacement for the barrnap + tRNAscan-SE step in
[YuggASMoth](../YuggASMoth/) (Module 1).  Use the combined GFF3 output
(`mod04_annotation_{prefix}.gff3`) directly as the rDNA/tRNA annotation input
to YuggASMoth's filtering step.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full version history.

**v0.1.0** _(2026-06-27)_ — initial release: telomere k-mer scan, tantan
masking, nhmmer (default) / cmsearch rRNA annotation, ARAGORN tRNA annotation,
GFF3 integration.

---

## Third-party tools and citations

| Tool | Reference |
|---|---|
| **HMMER3 / nhmmer** | Eddy SR (2011) *PLoS Comput Biol* 7:e1002195 |
| **Infernal / cmsearch** | Nawrocki EP, Eddy SR (2013) *Bioinformatics* 29:2933–2935 |
| **Rfam** | Kalvari I et al. (2021) *Nucleic Acids Res* 49:D192–D200 |
| **ARAGORN** | Laslett D, Canback B (2004) *Nucleic Acids Res* 32:11–16 |
| **tantan** | Frith MC (2011) *Nucleic Acids Res* 39:e23 |

---

## FAIR compliance

- **License:** MIT (see [LICENSE](LICENSE))
- **Citation:** see [CITATION.cff](CITATION.cff)
- **Zenodo DOI:** _mint after first public release_
- **bio.tools:** _register after first public release_
