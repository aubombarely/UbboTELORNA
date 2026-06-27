# Changelog — UbboTELORNA

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
