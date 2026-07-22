# Changelog — UbboTELORNA

## [v0.3.0] — 2026-07-22

### Added
- `--chunk_max_bp` (default 50,000,000): Module 2 rRNA chunking now also
  flushes a chunk once its cumulative sequence length reaches this many bp,
  in addition to the existing `--chunk_size` sequence-count threshold

### Fixed
- Module 2 chunking (`--chunk_size`) only flushed on sequence count.
  Chromosome-scale assemblies with few, very large sequences (e.g. maize)
  never reached the count threshold, so the whole genome still landed in a
  single nhmmer/cmsearch call despite chunking being enabled — the v0.2.0
  changelog's claim of eliminating the maize memory spike was incomplete.
  `--chunk_max_bp` (see Added) is the actual fix; verified against a
  synthetic multi-chromosome genome and confirmed backward-compatible when
  disabled (`--chunk_max_bp 0`)
- Benchmark comparator (`benchmark/scripts/compare_annotations.py`) ignored
  the `array_member`/`array_id` metadata that `--flag_5s_arrays` already
  produces, so every member of a correctly-detected 5S tandem array was
  scored as an independent false positive against a sparse reference.
  Array members sharing an `array_id` are now collapsed into one
  representative prediction before matching; verified on a synthetic
  20-copy array (F1: 0.174 → 1.0)
- Benchmark perf sweep (`benchmark/Snakefile`, `perf_ubbotelorna_sweep`) had
  no stale-output guard; reruns against an existing `_run` directory
  silently hit UbboTELORNA's own module-level checkpoints and reported
  near-instant, near-zero-memory results instead of a genuine measurement.
  Output directory is now cleared and `--force` is passed before every
  sweep run
- Benchmark report's peak-memory panel (`benchmark/scripts/benchmark_report.py`)
  read only the `threads=1` row per genome instead of the max across the
  full thread sweep; a single under-sampled reading (Snakemake's benchmark
  RSS polling can miss the true peak on very short jobs) could report a
  false near-zero outlier for an entire genome

### Changed
- Benchmark report's F1 panels (rRNA/tRNA vs. RefSeq) replaced with a
  connected scatter plot (genomes sorted by F1, one line per tool) instead
  of an annotated heatmap, for readability at 40-genome scale

## [v0.2.0] — 2026-07-04

### Added
- Module 2: 5S rRNA tandem array flagging (`--flag_5s_arrays`, on by default);
  consecutive 5S hits within `--5s_array_max_gap` bp are grouped into arrays,
  tagged with `array_member=true;array_id=arrNNNN` in GFF3 attributes, and
  summarised in `mod02_5s_arrays_{prefix}.tsv`
- Module 2: optional copy-number cap (`--cap_5s N`); retains only the N
  highest-scoring 5S predictions per sequence (off by default)
- New `--chunk_size` argument (default 200): feeds the genome to nhmmer/cmsearch
  in batches to keep peak RAM proportional to chunk size rather than genome size
- Streaming FASTA helpers (`_iter_fasta`, `_write_hard_masked_streaming`) that
  hold only one sequence in memory at a time — eliminates the ~1.2 GB spike
  at large genomes (e.g. maize)

### Changed
- `run_module0_telomeres` and `run_module1_masking` switched to streaming FASTA
  passes; full-genome dicts no longer held in RAM
- `run_module5_plot` split into four panel helpers; `run_module6_evolution` split
  into three analysis helpers; `_open_run_log` and `_setup_co2_tracker` extracted
  from `main()`
- Duplicate `_rev_comp` / `_revcomp` and double `_COMP` definitions consolidated
- Module 6 TSV outputs now sorted by seqname / start before writing

### Fixed
- Removed dead variables (`in_tel`, `tel_start_rel`) and duplicate `Counter`
  imports inside functions
- ARAGORN output file was opened twice per run; redundant open removed

## [v0.1.0] — 2026-06-27

### Added
- Module 0: telomere identification via k-mer density scan at contig ends
  with auto-detection of the repeat unit
- Module 1: low-complexity masking with `tantan` (soft-mask → hard-mask
  for search tools)
- Module 2: rRNA annotation using Infernal `cmsearch` with Rfam covariance
  models; supports eukaryote, bacteria, and archaea model sets; auto-downloads
  and caches CMs from Rfam
- Module 3: tRNA annotation using ARAGORN
- Module 4: integration of all GFF3 outputs into a combined annotation file
  with a summary TSV
- Blueprint-compliant logging, checkpoint/resume (`--force`), dry-run
  (`--dry_run`), run summary JSON, carbon footprint tracking (codecarbon),
  and resource usage reporting
- H.P. Lovecraft / Clark Ashton Smith Ubbo-Sathla quote banner
- `envs/UbboTELORNA.yaml` conda environment
- `LICENSE` (MIT), `CITATION.cff`
