# Changelog — UbboTELORNA

## [v0.10.0] — 2026-08-19

### Added
- `mod05_summary_{prefix}.tsv` now includes a `centromere` row once Module 6
  (Centromere detection) has run. Module 6 runs after Module 5 (Integration)
  writes this file, so its results weren't previously reflected in the
  genome-wide feature census at all -- a real gap noticed by a user
  inspecting the summary after their v0.9.2 Phillyrea run. The row counts
  only the primary `centromere_candidate` per sequence (mirroring what
  Module 7's plot shows); secondary `satellite_array_candidate` and
  `te_cluster_candidate` rows remain in `mod06_centromere_summary_*.tsv`
  only, to avoid double-counting overlapping bp in the genome-wide total.
  Appending is idempotent -- rerunning against an existing summary with a
  `centromere` row already present does not duplicate it.

### Fixed
- The module list in `--help`'s description was stale, omitting Module 6
  (Centromere detection) entirely and mislabeling Modules 7/8 as 6/7 --
  left over from before Module 6 was added. Corrected to match the
  top-of-file docstring and actual execution order.

## [v0.9.3] — 2026-08-19

### Added
- Panel labels (A, B, C, D) in the upper-left corner of each subplot in
  `mod07_plot_{prefix}` (ideogram, rRNA bars, tRNA bars, composition donut)
  for easier reference in figure legends/manuscripts.

### Fixed
- **Plot filenames silently truncated the run prefix whenever it contained
  a dot before the version-like suffix** (e.g. `..._v0.9.2.pdf` was written
  as `..._v0.9.pdf`). Confirmed on a real run (*Phillyrea angustifolia*,
  prefix `PhangAGP1Chromosomes.ubbotelorna_v0.9.2`): both
  `mod07_plot_*` and `mod08_evolution_*` outputs lost the `.2`. Root cause:
  `out_base.with_suffix(f".{fmt}")` replaces whatever follows the *last*
  dot in the path's name, and `out_base`'s name embeds the run prefix
  verbatim, so a prefix ending in something like `v0.9.2` has its own `.2`
  mistaken for an existing suffix and overwritten. Fixed by adding a
  `_with_format()` helper that appends the format suffix via
  `Path.with_name()` instead, and switching all four plot-output call
  sites (Module 7 and Module 8) to use it. This did not affect any
  non-plot output file, since those are built via f-string concatenation
  rather than `.with_suffix()`.

## [v0.9.2] — 2026-08-19

### Fixed
- **TE-family cluster selection ranked by raw copy count, letting a
  ubiquitous, uniformly-distributed family dominate the "best cluster"
  slot on nearly every chromosome.** Confirmed on a real genome
  (*Phillyrea angustifolia*, v0.9.0 run): a single family
  (`RND-3_FAMILY-2147`, ~8,800 copies genome-wide) "won" as the top TE
  cluster on 18 of 23 chromosomes, but at only 0.7-2.3% concentration
  each -- i.e. only a tiny fraction of that family's copies were
  actually in any one cluster, the opposite of the localized-
  amplification signature a real centromeric retrotransposon shows.
  This happened because a family common enough clears
  `--centromere_te_min_copies` in *some* window on almost every
  chromosome from sheer genome-wide abundance alone, and ranking by
  raw `n_copies` let that noise consistently beat genuinely
  concentrated (but lower absolute copy number) families.
  Fixed by computing `concentration_pct` for every qualifying cluster
  up front and ranking by concentration first, copy count as tiebreak,
  instead of the reverse. Also adds an optional
  `--centromere_te_min_concentration_pct` (default 0.0, off) to
  exclude low-concentration clusters outright rather than just
  deprioritizing them. Verified with a controlled synthetic case (a
  1,000-copy genome-wide family with only 20 copies, 2%, in one
  cluster vs. a 30-copy family entirely, 100%, in one cluster): the
  truly localized family now correctly wins, both by default ranking
  and with the new minimum-concentration filter.

## [v0.9.1] — 2026-08-19

### Fixed
- **`KeyError` crash in Module 6 when `--centromere_te_gff` covers a
  different sequence set than the scanned FASTA.** Real crash on
  *Phillyrea angustifolia*: an EarlGrey GFF3 generated against the full
  assembly (chromosomes + unplaced contigs) was supplied while running
  Module 6 against a chromosomes-only FASTA subset; an unplaced contig
  name from the GFF3 (`PhangAGP1UCTG024`) had no entry in `seq_lengths`
  (built from the scanned FASTA), and indexing into it directly crashed
  with `KeyError: 'PhangAGP1UCTG024'`. Fixed by filtering EarlGrey hits
  to sequences actually present in the scanned FASTA before clustering,
  logging how many GFF3 sequences were skipped for not matching (rather
  than silently dropping them with no indication). This is a legitimate
  use case, not user error -- TE annotation is often run once on a full
  assembly and reused across different downstream sequence subsets.
  Verified by reproducing the exact crash (a mixed EarlGrey GFF3 with
  both a real chromosome and an unplaced-contig-style sequence name
  scanned against a chromosomes-only FASTA) and confirming it no longer
  errors, correctly processes the matching chromosome's TE cluster, and
  logs the one skipped sequence.

## [v0.9.0] — 2026-08-19

### Fixed
- **Module 6 (centromere detection, formerly numbered 8) could report a
  large subtelomeric satellite array as the primary centromere
  candidate.** Confirmed on a real genome (*Phillyrea angustifolia*,
  v0.8.0 run): a 609 kb array on one chromosome ran to the literal last
  base of the sequence, and Module 1 had already independently flagged
  that same terminus as a Tier-2 subtelomeric repeat -- Module 8 had no
  way to know that and reported it as the primary candidate anyway.
  Large subtelomeric arrays produce exactly the same TRF signature
  (long, high-copy, tandem) as a real centromere.
  Fixed with `_flag_island_context()`: every candidate array is now
  checked against (a) proximity to either sequence end
  (`--centromere_end_buffer_bp`, default 100000) and (b) overlap with
  an already-confirmed telomere (Module 0) or subtelomeric repeat
  (Module 1) region for the same sequence. Flagged candidates are
  deprioritized for primary-candidate status (a clean, non-flagged
  island is preferred even if smaller) but never silently discarded --
  everything found is still reported, now with `is_suspect`/
  `near_chromosome_end`/`overlaps_telomere`/`overlaps_subtelomeric`
  attributes (GFF3) and columns (summary TSV) so the caveat is visible
  rather than something the user has to manually cross-check.

### Added
- **Optional EarlGrey TE-family clustering as a second, independent
  centromere signal.** Many plant genomes have retrotransposon-based
  centromeres rather than simple tandem satellite DNA, which the
  TRF-only scan cannot detect under any parameters. `--centromere_te_gff`
  accepts an EarlGrey repeat-annotation GFF3 (standard 9-column GFF3,
  `ID=RND-.._FAMILY-..` and `KIMURA80=..` in column 9); TE instances are
  grouped by family and clustered per sequence
  (`_merge_te_family_clusters()`, `--centromere_te_merge_gap_bp`,
  default 300000 -- more generous than the satellite merge gap, since
  retrotransposon insertions are more sparsely spaced than a tandem
  array). A family whose copies concentrate tightly in one narrow
  window (`--centromere_te_min_copies`, default 15) rather than being
  scattered genome-wide is the classic retrotransposon-centromere
  signature. Reported alongside the TRF-based candidate as a
  cross-validating signal (new `te_*` summary TSV columns plus a
  `te_cluster_candidate` GFF3 feature), including a `concentration_pct`
  (this cluster's share of that family's genome-wide copy count) and a
  `combined_evidence` flag (do the TRF and TE candidates for this
  sequence actually overlap -- the strongest possible signal when they
  do). `--centromere_te_max_kimura` optionally prioritizes recently-
  active (low-divergence) clusters, off by default so an older but real
  cluster isn't silently excluded.

### Changed
- **Module numbers now match execution order.** Previously, Module 8
  (centromere detection) executed between Modules 5 and 6 (so its
  output was ready in time for Module 6's plot), which meant module
  *numbers* no longer matched execution order for the first time in
  this tool. Renumbered so execution order and numbering agree exactly:
  `0`=telomere, `1`=subtelomeric, `2`=masking, `3`=rRNA, `4`=tRNA,
  `5`=integration, `6`=centromere detection (was 8),
  `7`=visualization (was 6), `8`=evolutionary analysis (was 7). No code
  actually moved -- centromere detection already ran right after
  integration and before visualization; only the labels and output
  filenames changed. **This renames output files**:
  `mod06_plot_*` → `mod07_plot_*`; `mod07_rrna_scores_*`/
  `mod07_trna_class_*`/`mod07_arrays_*`/`mod07_evolution_*` →
  `mod08_*`; `mod08_centromere_*` → `mod06_centromere_*`. Functions
  renamed to match: `run_module6_plot`→`run_module7_plot`,
  `run_module7_evolution`→`run_module8_evolution`,
  `run_module8_centromere`→`run_module6_centromere`. Also caught and
  fixed a `--skip_module` example in the README that would have run the
  wrong module after the renumbering (skipping 0-6 no longer isolates
  just the evolutionary-analysis module on its own, since visualization
  is now module 7, not 6).

## [v0.8.0] — 2026-08-18

### Added
- **Module 8: Centromere detection.** Unlike Modules 0/1 (telomere/
  subtelomeric), which only scan short terminal windows, centromere
  position is not known a priori and centromeric satellites are not
  conserved across species (no fixed motif to search for, unlike
  telomere repeats) -- so this scans each full sequence with TRF for
  candidate satellite arrays: long, high-copy tandem repeats in the
  typical plant centromeric monomer period range (default 50-200 bp,
  configurable via `--centromere_min_period`/`--centromere_max_period`).
  Overlapping/nearby hits are merged into contiguous "islands"
  (`--centromere_merge_gap_bp`, default 10000); islands below
  `--centromere_min_array_bp` (default 50000) are dropped. The single
  largest remaining island per sequence is flagged as the primary
  candidate (`centromere_candidate`); other qualifying islands are also
  reported (`satellite_array_candidate`), since a genome can carry more
  than one satellite family and the true centromere is not always the
  single largest array. `--centromere_min_seq_length` (default 0, off)
  restricts scanning to sequences at least that long, avoiding wasted
  effort on thousands of small unplaced scaffolds on fragmented
  assemblies. TRF has no built-in multithreading, so whole-sequence
  scans are dispatched one per sequence to a thread pool (same
  rationale as Module 3's chunk parallelisation in v0.6.1) rather than
  run sequentially -- each parallel TRF invocation gets its own
  subdirectory so concurrent runs can't race on the `*.dat` glob inside
  `_run_trf`.
  Output: `mod08_centromere_{prefix}.gff3` +
  `mod08_centromere_summary_{prefix}.tsv`.
- **Centromere track in the Module 6 plot**, as requested directly:
  the primary candidate per sequence is now shown on the ideogram (new
  purple `_C_CEN` layer, alongside the existing telomere/rRNA/tRNA
  tracks) and in the genome-composition donut chart. Kept as a direct
  GFF3 read in `run_module6_plot` (a new `_parse_gff3_positions_by_type`
  helper, filtering to just the primary `centromere_candidate` rows so
  secondary satellite arrays don't clutter the overview) rather than
  routed through Module 5's combined annotation/summary, so Module 8
  only needs to execute before Module 6, not disturb Modules 0-7's
  existing execution order -- it runs right after Module 5 for exactly
  that reason (the same reason Module 5 itself runs before Module 6).
  Verified with a full synthetic run: correct merge-island math (nearby
  hits merged, distant/small ones correctly separated/dropped),
  correct primary/secondary ranking by span, correct exclusion of
  sequences below `--centromere_min_seq_length`, and an actual rendered
  PNG confirming the new ideogram track and donut wedge both appear
  correctly.

## [v0.7.3] — 2026-08-18

### Added
- **`best_length_bp` column in `mod01_completeness_*.tsv`**: the exact
  length (bp) of the best (highest-copy) subtelomeric tandem repeat
  hit at each scaffold end, computed as `end - start + 1` from TRF's
  own reported span -- mirroring Module 0's `mod00_summary_*.tsv`,
  which already has an equivalent `length_bp` column for telomeres.
  Previously, the exact span was only recoverable by cross-referencing
  `mod01_subtelomeric_*.gff3` coordinates (columns 4-5); the
  completeness TSV itself only had `best_period_bp`/`best_copy_number`,
  from which length can only be approximated (period x copies is not
  exact, since real tandem arrays commonly have imperfect terminal
  copies). Verified against a synthetic case with two different
  spans (1-210bp and 5-62bp), confirming exact length is reported
  (210 and 58 respectively), not an approximation.

## [v0.7.2] — 2026-08-18

### Added
- **`motif_is_telomeric` column in `mod01_completeness_*.tsv`** (and an
  `is_telomeric_motif` GFF3 attribute per feature in
  `mod01_subtelomeric_*.gff3`), distinguishing a genuinely distinct
  subtelomeric satellite from Module 1 independently rediscovering the
  telomere array itself. Module 1's `best_motif` is selected purely
  from its own TRF scan of a wider terminal window, with no awareness
  of Module 0's telomere call — so on Tier 1 (confirmed-telomere) rows,
  `best_motif` is very often just a rotation/strand variant of the
  telomere repeat itself (the telomere array is typically the highest-
  copy tandem repeat in that window by a wide margin), not a separate
  finding. More importantly, this was also ambiguous on **Tier 2**
  rows: `best_motif` matching the telomere repeat family there means
  something biologically different (a likely degraded/partial telomere
  that fell just short of Module 0's confirmation threshold) than
  `best_motif` being an unrelated satellite (a genuine, distinct
  subtelomeric repeat) -- previously indistinguishable without manually
  comparing the motif string by eye. Comparison
  (`_motif_matches_telomere()`) is strand- and rotation-aware (reuses
  the same rotation-set approach as `_kmer_density()`), since e.g.
  `GGTTTAG` is a cyclic rotation of `TTTAGGG`, not simply its reverse
  complement -- a naive `_canonical()` (reverse-complement-only)
  comparison would have missed that case. Verified against a 3-scaffold
  synthetic case covering all three outcomes: confirmed telomere
  (`motif_is_telomeric=yes`), genuine distinct satellite
  (`=no`), and unconfirmed-but-telomere-family motif (`=yes` despite
  Tier 2) -- each classified correctly.

## [v0.7.1] — 2026-08-18

### Fixed
- **Module 1's "best" subtelomeric tandem repeat could be a meaningless
  homopolymer run.** `best_hit = max(hits, key=lambda h: h["copies"])`
  selected purely by raw TRF copy count, with no period filtering —
  and a period-1 homopolymer (e.g. a 29 bp run of `A`) trivially
  racks up a high copy count relative to its tiny period, so it could
  out-rank a real, longer-period satellite repeat present in the same
  window. Observed directly in a real run's
  `mod01_completeness_*.tsv`: `best_period_bp=1, best_motif=A` reported
  as the representative subtelomeric repeat at a scaffold end that also
  had a real period-5 candidate (`CGAAC`) with a lower (but still
  well-supported) copy count. Fixed by adding
  `--subtelomeric_min_period` (default: 2), filtering out period-1
  homopolymer hits before the max-copy-number selection, alongside the
  existing `--subtelomeric_min_copies` filter. Verified by reproducing
  the exact scenario (a period-1/29-copy homopolymer competing against
  a period-5/11.6-copy real satellite) and confirming the homopolymer
  is excluded and the real motif is correctly selected once the
  default filter is applied.

## [v0.7.0] — 2026-08-18

### Added
- **`--telomere_min_seq_length`** to guard telomere repeat-unit
  auto-detection against small-fragment noise on non-chromosome-scale
  assemblies. Auto-detection (`_select_telomere_candidate_trf`) picks
  whichever repeat unit is supported by TRF hits at the most *distinct
  sequences'* termini -- on a whole, fragmented genome assembly with far
  more small unplaced scaffolds than real chromosomes, most of those
  scaffold ends are simply assembly breakpoints (frequently inside
  unrelated repetitive DNA, since assemblers commonly stall there), and
  a non-telomeric repeat that happens to be common among many such
  breakpoints can outnumber and mask the true, narrowly-distributed
  telomere signal from the real chromosomes. Observed on a real genome:
  running on the full assembly (206 sequences) auto-detected
  `CCAGGACATGG` and found 0 confirmed telomeres, while running on the
  same genome's chromosome-scale pseudomolecules only (24 sequences)
  correctly detected the canonical plant telomere repeat and confirmed
  telomeres at 81.2% of ends. `--telomere_min_seq_length` restricts
  which sequences are eligible to vote on the auto-detected motif to
  those at least as long as the given threshold (default: 0, i.e. no
  filtering -- opt-in, since no single default works across organisms
  with very different chromosome sizes, e.g. ~5 Mb+ for most plant/
  animal genomes vs. ~200 kb+ for *S. cerevisiae*). Has no effect on
  `--telomere_repeat` (explicit repeat unit) or on the main telomere
  scan itself -- every sequence, regardless of length, is still scanned
  and reported once the repeat unit is known. Verified by unit-testing
  `run_module0_telomeres()` end-to-end (TRF's own binary mocked out)
  against a synthetic mix of 3 large "chromosome" and 10 small
  "fragment" sequences: unfiltered, the fragments' bogus motif won the
  vote exactly as in the real-genome case; with the threshold set above
  the fragment size, the real telomere motif won correctly.

## [v0.6.2] — 2026-08-18

### Fixed
- **Subtelomeric tandem repeat consensus motif was parsed but never
  written to any output file.** `_parse_trf_dat()`'s own docstring said
  the parsed consensus sequence was "used... for Module 1's GFF3
  attributes", but `run_module1_subtelomeric()` built the GFF3 `attrs`
  dict and TSV summary rows without it — so the actual repeat motif
  (e.g. `TTTAGGG`) was unrecoverable from `mod01_subtelomeric_*.gff3`
  or `mod01_completeness_*.tsv`, only the period length and copy number
  were reported. Fixed by adding `motif=` to the GFF3 attributes (one
  per hit) and a `best_motif` column to the completeness TSV (the
  motif of the highest-copy-number hit at that scaffold end, matching
  the existing `best_period_bp`/`best_copy_number` fields). Verified
  by unit-testing `_parse_trf_dat()` and the output-row construction
  directly against a synthetic TRF `.dat` record, since `trf` itself
  wasn't available in this environment to run Module 1 end-to-end.

## [v0.6.1] — 2026-07-26

### Fixed
- **Module 3 (rRNA) chunked search did not scale with `--threads` on
  large/many-scaffold genomes.** Real benchmark data (40-genome
  correctness/performance run) showed barrnap scaling near-linearly with
  threads (e.g. `gallus`: ~20.7x speedup from 1->32 threads) while
  UbboTELORNA's chunked nhmmer search barely scaled at all on the same
  genomes (`dmelanogaster`: only ~1.7x speedup 1->32 threads) — at 32
  threads, barrnap ended up **7-15x faster** than UbboTELORNA on the
  larger genomes tested (`gallus`, `osativa`, `dmelanogaster`), a reversal
  of the advantage UbboTELORNA has at low thread counts. Root cause,
  confirmed by reading the code: chunks were processed in a plain
  sequential Python loop, each chunk's `nhmmer`/`cmsearch` subprocess call
  blocking before the next chunk started — `--threads` only sped up the
  search *within* one chunk, never let chunk N+1 overlap with chunk N.
  - Fixed by dispatching chunks to a `concurrent.futures.ThreadPoolExecutor`
    (chunk search is I/O-bound on an external subprocess, so threads
    parallelise real wall-clock work without needing multiprocessing) —
    new `--chunk_workers` (default 4) controls how many chunks run
    concurrently, with `--threads` split `threads // chunk_workers` ways
    per chunk rather than all handed to one sequential chunk at a time.
  - Validated (mocked `_run_single_chunk`, no real nhmmer/cmsearch
    needed): all chunks processed exactly once with correct threads/chunk
    math, hits correctly aggregated across chunks into the final GFF3,
    genuine measured concurrency (6 simulated 0.3s chunks completed in
    ~0.6s with 3 concurrent workers, vs. ~1.8s if run serially), and a
    simulated chunk failure still correctly propagates as a fatal error
    (`SystemExit`) to the caller, matching the original fail-fast
    behaviour. **Not yet validated against the real 40-genome benchmark**
    — the numbers above motivated the fix but this specific change has
    not itself been re-benchmarked yet.

## [v0.6.0] — 2026-07-25

### Changed
- **Renumbered all modules to match their actual execution order in
  `main()`.** Module 7 (Subtelomeric tandem repeats, added in v0.4.0) has
  always run immediately after Module 0 (it depends on Module 0's telomere
  GFF3), but kept the number "7" — confusing anyone reading the module list
  top-to-bottom against the log output. The numbers now read 0-7 in the
  order the modules actually execute:
  - `0` Telomere identification — unchanged
  - `7` → `1` Subtelomeric tandem repeats
  - `1` → `2` Low-complexity masking
  - `2` → `3` rRNA annotation
  - `3` → `4` tRNA annotation
  - `4` → `5` Integration
  - `5` → `6` Visualization
  - `6` → `7` Evolutionary analysis
- **BREAKING: output filenames changed for every module except Module 0.**
  Anyone scripting against specific `results/mod0N_*` filenames from a
  prior run must update those references:
  `mod07_subtelomeric_*.gff3` / `mod07_completeness_*.tsv` →
  `mod01_subtelomeric_*.gff3` / `mod01_completeness_*.tsv`;
  `mod02_rRNA_*.gff3` → `mod03_rRNA_*.gff3`;
  `mod02_5s_arrays_*.tsv` → `mod03_5s_arrays_*.tsv`;
  `mod03_tRNA_*.gff3` → `mod04_tRNA_*.gff3`;
  `mod04_annotation_*.gff3` → `mod05_annotation_*.gff3`;
  `mod04_summary_*.tsv` → `mod05_summary_*.tsv`;
  `mod05_plot_*.{pdf,png,svg}` → `mod06_plot_*.{pdf,png,svg}`;
  `mod06_rrna_scores_*.tsv` → `mod07_rrna_scores_*.tsv`;
  `mod06_trna_class_*.tsv` → `mod07_trna_class_*.tsv`;
  `mod06_arrays_*.tsv` → `mod07_arrays_*.tsv`;
  `mod06_evolution_*.{pdf,png,svg}` → `mod07_evolution_*.{pdf,png,svg}`.
  `mod00_telomeres_*.gff3` and `mod00_summary_*.tsv` (Module 0) are
  unchanged, as are the unprefixed `workdir/masked_soft.fasta` /
  `workdir/masked_hard.fasta` (Module 2 masking). Existing output
  directories from prior runs are **not** migrated automatically — rerun
  with `--force`, or manually rename files, if you depend on the new
  naming scheme.
- **`--skip_module` numbering changed to match**: what was
  `--skip_module 7` (subtelomeric) is now `--skip_module 1`, and so on per
  the mapping above. Update any saved commands or scripts accordingly.
- Renamed internal functions to match: `run_module7_subtelomeric` →
  `run_module1_subtelomeric`, `run_module1_masking` → `run_module2_masking`,
  `run_module2_rrna` → `run_module3_rrna`, `run_module3_trna` →
  `run_module4_trna`, `run_module4_integration` →
  `run_module5_integration`, `run_module5_plot` → `run_module6_plot`,
  `run_module6_evolution` → `run_module7_evolution`. No functional changes
  — labels, filenames, banners, and argparse group titles only.

## [v0.5.0] — 2026-07-25

### Changed
- **Module 0 telomere repeat-unit auto-detection rewritten around TRF**,
  replacing the from-scratch Python k-mer scanner (low-complexity filtering
  + contiguous-run counting). Observed on a real Citrus sinensis genome:
  the k-mer scanner picked a compositionally-generic-but-spurious motif
  (`AATAA`) over the true telomere repeat (`TTTAGGG`) — "most sequences
  with some support" doesn't distinguish a real, terminus-concentrated
  signal from AT-rich background noise diffused across the (10,000bp)
  scan window — and took over two hours to do it (exhaustively checking
  every complexity-filtered candidate's contiguous run against every
  sequence does not scale). TRF (compiled C, purpose-built for
  unknown-period tandem repeat detection, already a dependency via
  Module 7) is both more robust to imperfect/degenerate repeat copies
  and dramatically faster than an exact-match k-mer scan.
  - New `_select_telomere_candidate_trf()`: among TRF hits whose period
    falls in the canonical telomere range (4-12bp, distinguishing it from
    longer-period subtelomeric satellites), picks the consensus sequence
    supported by the most distinct sequences.
  - New `--telomere_detect_window` (default 5000bp) controls the
    auto-detection scan window, deliberately smaller than
    `--telomere_window` (used for extent-scanning once the repeat unit is
    known) — detection close to the true terminus avoids diffuse
    subtelomeric/intergenic sequence outcompeting the real signal.
  - Removed now-dead code: `_detect_repeat_unit`, `_is_low_complexity`,
    `_kmer_max_run` (the from-scratch scanner and its helpers).
  - `_parse_trf_dat` (shared with Module 7) now also captures the TRF
    consensus repeat sequence.

### Added
- **Module 0 summary table** (`mod00_summary_{prefix}.tsv`): one row per
  scaffold end (seqname, end, found, start, end_pos, length_bp, density,
  repeat_unit) — previously Module 0 only wrote a GFF3, with no
  standalone aggregate view (unlike Module 7's completeness table).
  Works independently of Module 7.

### Validated
- `_select_telomere_candidate_trf`: correctly prefers a short-period,
  multi-sequence-supported telomere motif over a long-period subtelomeric
  satellite in the same hit set, and correctly returns `None` when support
  is below `min_seq_support`.
- Full `run_module0_telomeres` orchestration (TRF stubbed, same pattern as
  Module 7's tests): correctly detects the repeat unit from synthetic TRF
  hits, writes correct GFF3 + summary TSV, and correctly reports
  `found=no` for a scaffold with no telomere.
- Re-ran Module 7's existing orchestration test to confirm the shared TRF
  parsing helpers (now returning an additional `consensus` field) didn't
  regress Module 7's behavior.
- **Not yet validated against the real TRF binary for Module 0
  specifically** (only Module 7's real-TRF run has been confirmed so far) —
  test on real data before trusting auto-detection results.

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
