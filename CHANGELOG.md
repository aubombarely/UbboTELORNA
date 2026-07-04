# Changelog — UbboTELORNA

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
