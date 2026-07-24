# Changelog — UbboTELORNA

## [v0.4.0] — 2026-07-24

### Added
- **Module 7 — Subtelomeric tandem repeats**: scans the terminal window of
  every scaffold end with TRF (Tandem Repeats Finder) and classifies each
  end into a telomere-completeness tier: Tier 1 (Module 0 confirmed a
  telomere), Tier 2 (no confirmed telomere, but a genuine tandem repeat
  array was found nearby — some evidence of proximity to a chromosome end),
  or Tier 3 (neither). Motivation: subtelomeric satellite/tandem-repeat
  arrays are commonly (not universally) found adjacent to true telomeres in
  eukaryotic genomes, and can still be present even where the assembly
  stops just short of a fully resolved canonical telomere array, or where
  that array is too degraded for Module 0's strict k-mer scan to call
  confidently — giving a coarse, non-authoritative signal of how close an
  unresolved scaffold end may be to a true chromosome terminus.
  - New outputs: `results/mod07_subtelomeric_{prefix}.gff3`
    (`subtelomeric_tandem_repeat` features) and
    `results/mod07_completeness_{prefix}.tsv` (one row per scaffold end:
    seqname, end, length, tier, confirmed_telomere, best TRF hit period/copy
    number).
  - New `run_summary.json` section `telomere_completeness`, with a
    genome-wide rollup of scaffold ends per tier (count and %).
  - New flags: `--subtelomeric_window_bp` (default 20,000 — larger than
    `--telomere_window` since satellite arrays can sit further from the
    true terminus than the telomere repeat itself), `--subtelomeric_min_copies`
    (default 3.0 — minimum TRF-reported tandem copy number to count as a hit).
  - Runs immediately after Module 0 (its only dependency — needs Module 0's
    GFF3 to know which scaffold ends already have a confirmed telomere) and
    is included in `--skip_module`/`--dry_run` like every other module.
  - New dependency: `trf` (added to `envs/UbboTELORNA.yaml`). TRF is invoked
    directly rather than through the shared `_run()` helper, since it
    returns a non-zero exit code on ordinary success (a known TRF quirk) —
    success is instead detected by the presence of its `.dat` output file.
  - Validated: window-extraction/offset logic, `.dat` parsing, and the full
    tiering orchestration were unit-tested directly (including an
    end-to-end run with the TRF subprocess call stubbed out, confirming
    correct GFF3/TSV/rollup output for all three tiers). The real TRF
    subprocess invocation itself is **not yet validated against the actual
    binary** — TRF is not available in this development environment; test
    on a real genome before relying on it.

## [v0.3.1] — 2026-07-24

### Fixed
- Module 0 telomere repeat-unit auto-detection (`_detect_repeat_unit`) could
  return a degenerate homopolymer (e.g. `AAAAA`) instead of the real
  telomeric motif, causing Module 0 to report zero telomeres on genomes that
  genuinely have them (observed on a real Citrus sinensis assembly:
  `--telomere_repeat TTTAGGG` found real, high-density hits — including a
  clean 2,233bp array at the very start of one scaffold — that auto-detection
  had missed entirely). Root cause: the previous implementation pooled raw
  k-mer counts across all sequences and k-mer lengths with no complexity
  filter, so a trivial homopolymer could structurally outcompete a real but
  sparser telomeric signal. Fixed by (1) excluding low-complexity/homopolymer
  k-mers from candidacy outright (`_is_low_complexity`), and (2) requiring
  the winning k-mer to show a genuine unbroken tandem run
  (`_kmer_max_run`, >= `min_tandem_copies` consecutive copies) in at least
  `min_seq_support` distinct sequences, rather than picking whatever is
  globally most frequent when pooled. Validated against synthetic scenarios
  with realistic nucleotide composition, including a direct reproduction of
  the original bug (homopolymer-padded sequence ends, no real telomere
  anywhere) and a sparse-signal case (telomere present in only 2 of 60
  sequences).

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
