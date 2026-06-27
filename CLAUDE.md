# UbboTELORNA — Project Notes

Annotates three classes of ancient, conserved genomic elements in genome
assemblies: telomeres (Module 0), rRNA genes (Module 2), and tRNA genes
(Module 3). Includes low-complexity masking (Module 1) to prevent annotation
tool failures on telomeric and repetitive sequences.

Named after Ubbo-Sathla (Clark Ashton Smith / Lovecraft Mythos) — the
primordial source of all terrestrial life.

**Current version:** v0.1.0 — `scripts/UbboTELORNA.py`

This project fully follows the shared coding blueprint at `../CLAUDE.md`.
Apply those standards to any changes or additions here.

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/UbboTELORNA.py` | Main pipeline (Modules 0–4) |

## External tools required

| Tool | Module | Install |
|---|---|---|
| `tantan` | 1 — Low-complexity masking | `conda install -c bioconda tantan` |
| `cmsearch` (Infernal) | 2 — rRNA annotation | `conda install -c bioconda infernal` |
| `aragorn` | 3 — tRNA annotation | `conda install -c bioconda aragorn` |
| `codecarbon` | optional | `conda install -c conda-forge codecarbon` |

## Key design notes

- **Telomere detection** uses a pure-Python k-mer density scan (no external
  dependencies). The repeat unit is auto-detected from terminal windows or
  supplied via `--telomere_repeat`.
- **Rfam CMs** are downloaded on first use from `https://rfam.org/family/{acc}/cm`
  and cached in `~/.ubbotelorna/rfam/`. Provide `--rfam_dir` to use local copies
  (required on air-gapped HPC nodes).
- **Masking strategy**: tantan produces a soft-masked FASTA (lowercase).
  A hard-masked version (lowercase → N) is created in `workdir/` for use
  by cmsearch and ARAGORN, preventing search failures on repetitive sequences.
- **Kingdom flag** (`--kingdom euka/bacteria/archaea`) selects the appropriate
  Rfam model set for rRNA annotation.

---

## FAIR compliance status

- [x] `LICENSE` — MIT, 2026
- [x] `CITATION.cff` — author, ORCID, version, keywords, repository URL
- [ ] Zenodo DOI — mint after first public release
- [ ] bio.tools registration — register after first public release
