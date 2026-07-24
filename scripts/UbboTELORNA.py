#!/usr/bin/env python3
"""
UbboTELORNA.py  —  Telomere, rRNA, and tRNA annotation for genome assemblies.

Annotates three classes of ancient, conserved genomic elements:

  Module 0  Telomere identification  — k-mer density scan at contig ends
  Module 1  Low-complexity masking   — tantan soft-mask (prevents search failures
                                       on telomeric and repetitive sequences)
  Module 2  rRNA annotation          — nhmmer (default) or cmsearch + Rfam profiles
  Module 3  tRNA annotation          — ARAGORN
  Module 4  Integration              — merged GFF3 + summary table
  Module 5  Visualization            — ideogram, subtype bars, composition donut
  Module 6  Evolutionary analysis    — rRNA scores, tRNA pseudogenes, tandem arrays
  Module 7  Subtelomeric tandem repeats — TRF scan of scaffold-end windows;
                                       corroborates Module 0 telomere calls and
                                       classifies every scaffold end into a
                                       telomere-completeness tier

Named after Ubbo-Sathla (Clark Ashton Smith / H.P. Lovecraft Mythos), the
primordial source of all terrestrial life — mirroring telomeres, rDNA, and
tRNA as the most ancient conserved elements of eukaryotic genomes.
"""

import argparse
import getpass
import json
import os
import platform
import re
import resource
import shutil
import subprocess
import sys
import time
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

matplotlib.rcParams.update({
    "font.size":        10,
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "figure.dpi":       150,
    "savefig.dpi":      150,
    "figure.facecolor": "white",
})

VERSION = "v0.4.0"

# ── Rfam covariance model registry ────────────────────────────────────────────

_RFAM_CMS = {
    "euka": {
        "5S_rRNA":        "RF00001",
        "5_8S_rRNA":      "RF00002",
        "SSU_rRNA_eukarya": "RF01960",
        "LSU_rRNA_eukarya": "RF02543",
    },
    "bacteria": {
        "5S_rRNA":        "RF00001",
        "SSU_rRNA_bacteria": "RF00177",
        "LSU_rRNA_bacteria": "RF02541",
    },
    "archaea": {
        "5S_rRNA":        "RF00001",
        "SSU_rRNA_archaea": "RF01959",
        "LSU_rRNA_archaea": "RF02540",
    },
}
_RFAM_CM_URL  = "https://rfam.org/family/{acc}/cm"
_RFAM_STO_URL = "https://rfam.org/family/{acc}/alignment?format=stockholm&download=1"
_DEFAULT_CACHE = Path.home() / ".ubbotelorna" / "rfam"

# ── Telomere repeat defaults ───────────────────────────────────────────────────

_PLANT_TELOMERE     = "TTTAGGG"
_VERTEBRATE_TELOMERE = "TTAGGG"

# ── Logging infrastructure ─────────────────────────────────────────────────────

_LOG_FH = None


def _log(msg: str) -> None:
    ts   = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr)
    if _LOG_FH is not None:
        print(line, file=_LOG_FH, flush=True)


def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


_QUOTE_LINES = [
    "\"For Ubbo-Sathla is the source and the end. Before the coming of",
    " Zhothaqquah or Yok-Zothoth or Kthulhut from the stars, Ubbo-Sathla",
    " dwelt in the steaming fens of the new-made Earth: a mass without head",
    " or members, spawning the grey, formless efts of the prime and the",
    " grisly prototypes of terrene life.\"",
    "              — Clark Ashton Smith, Ubbo-Sathla (1933)",
]


def _print_quote() -> None:
    width = max(len(l) for l in _QUOTE_LINES) + 4
    border = "─" * width
    _log(f"┌{border}┐")
    for line in _QUOTE_LINES:
        padding = width - len(line) - 1
        _log(f"│ {line}{' ' * padding}│")
    _log(f"└{border}┘")


# ── External tool helpers ──────────────────────────────────────────────────────

def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(f"ERROR: '{name}' not found in PATH.\n"
              f"       Install with:  conda install -c bioconda {name}",
              file=sys.stderr)
        sys.exit(1)
    return tool


def _run(cmd: list, capture_stdout: bool = False,
         env: dict = None, cwd: Path = None) -> subprocess.CompletedProcess:
    _log(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd,
                            stdout=subprocess.PIPE if capture_stdout else None,
                            stderr=subprocess.PIPE, text=True,
                            env=env, cwd=cwd)
    if result.returncode != 0:
        print(f"ERROR: command failed (exit {result.returncode}):\n"
              f"{result.stderr[-3000:]}", file=sys.stderr)
        sys.exit(1)
    return result


def _checkpoint(path: Path, label: str, force: bool) -> bool:
    if not force and path.exists() and path.stat().st_size > 0:
        _log(f"  [checkpoint] {label} — {path.name} already exists, skipping")
        return True
    return False


# ── Input validation ──────────────────────────────────────────────────────────

def _validate_inputs(pairs: list) -> None:
    ok = True
    for flag, path in pairs:
        if not path.exists():
            print(f"ERROR: {flag} not found: {path}", file=sys.stderr)
            ok = False
    if not ok:
        sys.exit(1)


# ── FASTA utilities ───────────────────────────────────────────────────────────

def _fasta_total_length(fasta: Path) -> int:
    """Count total nucleotides in a FASTA file without loading sequences."""
    total = 0
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.rstrip())
    return total


def _fasta_seq_lengths(fasta: Path) -> dict:
    """Return {seq_name: length} preserving FASTA order."""
    lengths: dict = {}
    name = None
    bp   = 0
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = bp
                name = line[1:].split()[0]
                bp   = 0
            else:
                bp += len(line)
    if name is not None:
        lengths[name] = bp
    return lengths


def _parse_gff3_positions(gff: Path | None) -> dict:
    """Return {seq_name: [(start, end), ...]} from a GFF3, skipping comment lines."""
    groups: dict = {}
    if gff is None or not gff.exists():
        return groups
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            try:
                seqn = cols[0]
                groups.setdefault(seqn, []).append((int(cols[3]), int(cols[4])))
            except ValueError:
                continue
    return groups


def _parse_summary_tsv(tsv: Path | None) -> dict:
    """Return {feature_type: {subtype: (count, length_bp)}} from summary TSV."""
    data: dict = {}
    if tsv is None or not tsv.exists():
        return data
    with open(tsv) as fh:
        fh.readline()   # skip header
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue
            ftype, subtype = cols[0], cols[1]
            try:
                count  = int(cols[2])
                length = int(cols[3])
            except ValueError:
                continue
            data.setdefault(ftype, {})[subtype] = (count, length)
    return data


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header, sequence) from a FASTA file."""
    records = []
    header  = None
    seqs    = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seqs)))
                header = line[1:].split()[0]
                seqs   = []
            else:
                seqs.append(line)
    if header is not None:
        records.append((header, "".join(seqs)))
    return records


def _write_fasta(records: list[tuple[str, str]], path: Path,
                 line_width: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


def _iter_fasta(path: Path):
    """Yield (header, sequence) one record at a time — O(1) peak memory per sequence.

    Use instead of _read_fasta when the genome is large and only a single forward
    pass is needed (e.g. telomere scanning, hard-masking).  Can be iterated more
    than once by calling _iter_fasta(path) again for each pass.
    """
    header = None
    parts: list = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].split()[0]
                parts  = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def _write_hard_masked_streaming(soft_fasta: Path, hard_fasta: Path) -> None:
    """Stream a soft-masked FASTA and write the hard-masked (lowercase → N) version.

    Processes one line at a time so the entire genome is never held in memory.
    """
    with open(soft_fasta) as src, open(hard_fasta, "w") as dst:
        for line in src:
            if line.startswith(">"):
                dst.write(line)
            else:
                dst.write(re.sub(r"[acgt]", "N", line.rstrip("\n")) + "\n")


# ── GFF3 utilities ────────────────────────────────────────────────────────────

_GFF3_HEADER = "##gff-version 3\n"


def _gff3_length(line: str) -> int:
    """Return feature length (end - start + 1) from a GFF3 data line."""
    cols = line.split("\t")
    try:
        return int(cols[4]) - int(cols[3]) + 1
    except (IndexError, ValueError):
        return 0


def _gff3_name(line: str) -> str:
    """Extract Name= attribute from a GFF3 line; strips trailing anticodon parentheses."""
    attrs = line.split("\t")[8] if line.count("\t") >= 8 else ""
    for field in attrs.split(";"):
        if field.startswith("Name="):
            name = field[5:].strip()
            return re.sub(r"\([^)]+\)$", "", name)
    return ""


def _gff3_sort_key(line: str) -> tuple:
    """Natural-sort key for a GFF3 line: (seqid digit/alpha parts, start)."""
    cols = line.split("\t")
    seqid = cols[0] if cols else ""
    start = int(cols[3]) if len(cols) > 3 and cols[3].isdigit() else 0
    parts = [int(c) if c.isdigit() else c.lower()
             for c in re.split(r"(\d+)", seqid)]
    return (parts, start)


def _sort_hits(hits: list[dict]) -> list[dict]:
    """Sort a hit-dict list by (seqname natural, start)."""
    def key(h):
        parts = [int(c) if c.isdigit() else c.lower()
                 for c in re.split(r"(\d+)", h["seqname"])]
        return (parts, h["start"])
    return sorted(hits, key=key)


def _gff3_record(seqname: str, source: str, feature: str,
                 start: int, end: int, score: float | str,
                 strand: str, frame: str, attrs: dict) -> str:
    score_str = f"{score:.2f}" if isinstance(score, (int, float)) else "."
    attr_str  = ";".join(f"{k}={v}" for k, v in attrs.items())
    return f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score_str}\t{strand}\t{frame}\t{attr_str}"


# ── Telomere detection helpers ────────────────────────────────────────────────

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _canonical(kmer: str) -> str:
    """Lexicographically smaller of kmer and its reverse complement."""
    rc = _revcomp(kmer.upper())
    return min(kmer.upper(), rc)


def _all_rotations(kmer: str) -> set[str]:
    """All cyclic rotations of a k-mer (for repeat density counting)."""
    k = len(kmer)
    ku = kmer.upper()
    return {ku[i:] + ku[:i] for i in range(k)}


def _kmer_density(seq: str, repeat_unit: str) -> float:
    """Fraction of bases in seq covered by the repeat unit (any rotation/strand)."""
    ku   = repeat_unit.upper()
    rots = _all_rotations(ku) | _all_rotations(_revcomp(ku))
    k    = len(ku)
    n    = len(seq)
    if n < k:
        return 0.0
    covered = 0
    i = 0
    su = seq.upper()
    while i <= n - k:
        if su[i:i + k] in rots:
            covered += k
            i += k
        else:
            i += 1
    return covered / n


def _kmer_max_run(seq: str, repeat_unit: str) -> int:
    """Length (bp) of the longest *unbroken* run of consecutive repeat-unit
    matches (any rotation/strand) in seq. Unlike _kmer_density (which sums
    coverage across the whole window regardless of clustering), this is the
    real signature of a genuine tandem array as opposed to scattered
    coincidental matches: compositional noise can rack up nontrivial
    cumulative density from isolated hits spread across a window, but it
    essentially never produces a long *unbroken* run of consecutive copies
    the way an actual telomeric repeat array does."""
    ku   = repeat_unit.upper()
    rots = _all_rotations(ku) | _all_rotations(_revcomp(ku))
    k    = len(ku)
    n    = len(seq)
    if n < k:
        return 0
    su = seq.upper()
    i = 0
    best_run = 0
    cur_run = 0
    while i <= n - k:
        if su[i:i + k] in rots:
            cur_run += k
            best_run = max(best_run, cur_run)
            i += k
        else:
            cur_run = 0
            i += 1
    return best_run


def _is_low_complexity(kmer: str) -> bool:
    """True if kmer is dominated by a single-character run (homopolymer or
    near-homopolymer, e.g. AAAAA or AAAAT). These trivially over-count in
    AT-rich terminal sequence regardless of whether a genuine tandem repeat
    is present, and would otherwise structurally outcompete a real but
    shorter-period telomeric motif in the k-mer count."""
    k = len(kmer)
    max_run = 1
    run = 1
    for i in range(1, k):
        if kmer[i] == kmer[i - 1]:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 1
    return max_run >= k - 1


def _detect_repeat_unit(seqs,
                        window: int,
                        k_range: tuple[int, ...] = (5, 6, 7, 8),
                        min_seq_support: int = 2,
                        min_tandem_copies: int = 4) -> str | None:
    """Return the dominant telomere repeat unit from terminal k-mer analysis,
    or None.

    A naive "most frequent k-mer" approach is vulnerable to compositional
    bias: in AT-rich terminal sequence, a short, non-homopolymer k-mer can
    accumulate a high raw occurrence count -- or even high *cumulative*
    density -- purely by chance, from scattered hits spread across a window,
    without reflecting a genuine tandemly repeated array. The real signature
    of a telomeric repeat is an *unbroken run* of several consecutive copies
    in a row, which coincidental compositional matches essentially never
    produce. So every candidate that survives the low-complexity/homopolymer
    filter is ranked by its longest contiguous run of consecutive matches
    (`_kmer_max_run`), and the winner must reach at least
    `min_tandem_copies` consecutive copies (i.e. an unbroken run of
    `min_tandem_copies * len(kmer)` bp) in at least `min_seq_support`
    distinct sequences -- so a telomere present on only a handful of
    scaffolds (the common case in draft-quality assemblies) isn't
    outcompeted by compositional noise pooled across many non-telomeric
    sequence ends.
    """
    # Stage 1 (cheap): track which distinct sequences each complexity-filtered
    # candidate appears in at all (raw presence, not yet checked for tandem
    # structure). This is a necessary-but-not-sufficient pre-filter that
    # prunes the vast majority of one-off candidates -- with realistic
    # (diverse) sequence composition, the number of distinct candidate
    # k-mers surviving just the complexity filter can be very large, and
    # running the expensive contiguous-run scan (stage 2) against every one
    # of them, against every sequence, does not scale.
    candidate_presence: dict[str, set[str]] = {}
    # Bounded memory: only the (small) terminal windows are retained per
    # sequence, not the full sequence, preserving the streaming/one-genome-
    # in-memory-at-a-time design intent of the caller.
    terminal_regions: list[tuple[str, list[str]]] = []

    for name, seq in seqs:
        w    = min(window, len(seq))
        ends = [seq[:w]]
        if len(seq) > w:
            ends.append(seq[-w:])
        terminal_regions.append((name, ends))
        for region in ends:
            ru = region.upper()
            for k in k_range:
                for i in range(len(ru) - k + 1):
                    kmer = ru[i:i + k]
                    if set(kmer) <= set("ACGT") and not _is_low_complexity(kmer):
                        candidate_presence.setdefault(_canonical(kmer), set()).add(name)

    if not candidate_presence:
        return None

    support_floor = min(min_seq_support, len(terminal_regions))
    # Stage 2 (expensive, only for candidates that cleared stage 1): confirm
    # a genuine unbroken tandem run, not just scattered presence.
    shortlist = [k for k, names in candidate_presence.items() if len(names) >= support_floor]

    best_kmer, best_support = None, 0
    for kmer in shortlist:
        min_run_bp = min_tandem_copies * len(kmer)
        support = sum(
            1 for _name, ends in terminal_regions
            if any(_kmer_max_run(region, kmer) >= min_run_bp for region in ends)
        )
        if support > best_support:
            best_kmer, best_support = kmer, support

    if best_kmer is None or best_support < support_floor:
        return None
    return best_kmer


# ── Module 0: Telomere identification ─────────────────────────────────────────

def run_module0_telomeres(fasta: Path, repeat_unit: str | None,
                          tel_window: int, tel_density: float,
                          tel_min_len: int,
                          results: Path, workdir: Path,
                          prefix: str, force: bool) -> tuple[Path, str | None]:
    """Scan contig ends for telomeric repeats; return (gff3_path, repeat_unit_used)."""
    out_gff3 = results / f"mod00_telomeres_{prefix}.gff3"
    if _checkpoint(out_gff3, "telomere-scan", force):
        # Try to recover the repeat unit from existing GFF3
        with open(out_gff3) as fh:
            for line in fh:
                if line.startswith("##repeat-unit"):
                    return out_gff3, line.split()[-1]
        return out_gff3, repeat_unit

    # Count sequences with a single cheap pass (no sequence data loaded)
    n_seqs = sum(1 for line in open(fasta) if line.startswith(">"))
    _log(f"  Found {n_seqs} sequences")

    if repeat_unit is None:
        _log("  Auto-detecting telomere repeat unit …")
        # Pass 1: streaming auto-detection — one sequence in memory at a time
        repeat_unit = _detect_repeat_unit(_iter_fasta(fasta), tel_window)
        if repeat_unit is None:
            _log("  WARNING: could not auto-detect repeat unit. "
                 "Provide --telomere_repeat explicitly to enable telomere annotation.")
            out_gff3.write_text(_GFF3_HEADER)
            return out_gff3, None
        _log(f"  Detected repeat unit: {repeat_unit}")
    else:
        repeat_unit = repeat_unit.upper()
        _log(f"  Using user-supplied repeat unit: {repeat_unit}")

    records   = []
    tel_count = 0

    # Pass 2 (or only pass): stream sequences for the telomere scan
    for name, seq in _iter_fasta(fasta):
        length = len(seq)
        # Scan 5' end
        for end_label, region_seq, offset in [
            ("5prime", seq[:tel_window], 0),
            ("3prime", seq[max(0, length - tel_window):], max(0, length - tel_window)),
        ]:
            region_seq = region_seq if len(region_seq) >= len(repeat_unit) else ""
            if not region_seq:
                continue

            # Sliding window to find telomere extent
            step  = max(1, len(repeat_unit))
            wsize = len(repeat_unit) * 10   # window = 10 repeats wide
            rlen  = len(region_seq)

            densities = []
            for i in range(0, rlen - wsize + 1, step):
                d = _kmer_density(region_seq[i:i + wsize], repeat_unit)
                densities.append((i, d))

            # Merge contiguous windows above threshold
            tel_spans = []
            span_start = None
            for i, d in densities:
                if d >= tel_density:
                    if span_start is None:
                        span_start = i
                    span_end = i + wsize
                elif span_start is not None:
                    tel_spans.append((span_start, span_end))
                    span_start = None
            if span_start is not None:
                tel_spans.append((span_start, rlen))

            for s, e in tel_spans:
                abs_s = offset + s + 1   # 1-based GFF3
                abs_e = offset + e
                if abs_e - abs_s + 1 < tel_min_len:
                    continue
                tel_count += 1
                avg_d = _kmer_density(seq[abs_s - 1:abs_e], repeat_unit)
                attrs = {
                    "ID":          f"tel_{name}_{end_label}_{tel_count}",
                    "Name":        f"telomere_{end_label}",
                    "repeat_unit": repeat_unit,
                    "density":     f"{avg_d:.3f}",
                }
                records.append(
                    _gff3_record(name, "UbboTELORNA", "telomere",
                                 abs_s, abs_e, avg_d, ".", ".", attrs)
                )

    records.sort(key=_gff3_sort_key)
    with open(out_gff3, "w") as fh:
        fh.write(_GFF3_HEADER)
        fh.write(f"##repeat-unit {repeat_unit}\n")
        for r in records:
            fh.write(r + "\n")

    _log(f"  Telomeres found: {tel_count}  →  {out_gff3.name}")
    return out_gff3, repeat_unit


# ── Module 1: Low-complexity masking ─────────────────────────────────────────

def run_module1_masking(fasta: Path, workdir: Path, force: bool) -> Path:
    """Soft-mask low-complexity regions with tantan; return path to masked FASTA."""
    masked = workdir / "masked_soft.fasta"
    if _checkpoint(masked, "tantan-masking", force):
        return masked

    tantan = _require_tool("tantan")

    # Stream tantan stdout directly to disk — avoids holding the full genome
    # as a Python string in result.stdout (which can be 2-3 GB for large genomes).
    _log(f"  $ {tantan} {fasta.resolve()} > {masked}")
    with open(masked, "w") as _out_fh:
        _proc = subprocess.Popen(
            [tantan, str(fasta.resolve())],
            stdout=_out_fh,
            stderr=subprocess.PIPE,
            text=True,
        )
        _, _stderr = _proc.communicate()
    if _proc.returncode != 0:
        print(f"ERROR: tantan failed (exit {_proc.returncode}):\n"
              f"{_stderr[-3000:]}", file=sys.stderr)
        sys.exit(1)
    _log(f"  Soft-masked FASTA: {masked.name}")

    # Create hard-masked version (lowercase → N) for search tools.
    # Uses streaming to avoid a second full-genome load into memory.
    hard = workdir / "masked_hard.fasta"
    _log(f"  Creating hard-masked FASTA for search tools: {hard.name}")
    _write_hard_masked_streaming(masked, hard)

    return masked


def _hard_masked(workdir: Path) -> Path:
    return workdir / "masked_hard.fasta"


# ── Module 2: rRNA annotation (Infernal + Rfam) ───────────────────────────────

def _ensure_cms(kingdom: str, rfam_dir: Path | None) -> Path:
    """Download (if needed) and return path to a concatenated .cm file for the kingdom."""
    cache = rfam_dir or _DEFAULT_CACHE
    cache.mkdir(parents=True, exist_ok=True)

    models = _RFAM_CMS[kingdom]
    combined = cache / f"ubbotelorna_{kingdom}.cm"

    if combined.exists():
        _log(f"  Using cached Rfam CMs: {combined}")
        return combined

    _log(f"  Downloading Rfam CMs for kingdom '{kingdom}' …")
    cm_parts = []
    for name, acc in models.items():
        cm_path = cache / f"{acc}.cm"
        if not cm_path.exists():
            url = _RFAM_CM_URL.format(acc=acc)
            _log(f"    Downloading {acc} ({name}) from {url}")
            try:
                urlretrieve(url, cm_path)
            except Exception as exc:
                print(f"ERROR: failed to download {url}: {exc}\n"
                      f"       Check internet access or provide --rfam_dir with "
                      f"pre-downloaded .cm files.", file=sys.stderr)
                sys.exit(1)
        cm_parts.append(cm_path.read_text())

    combined.write_text("\n".join(cm_parts))
    _log(f"  CMs written to {combined}")
    return combined


def _ensure_hmms(kingdom: str, rfam_dir: Path | None) -> Path:
    """Build (if needed) and return path to a concatenated HMMER3 .hmm file for the kingdom."""
    cache = rfam_dir or _DEFAULT_CACHE
    cache.mkdir(parents=True, exist_ok=True)

    models   = _RFAM_CMS[kingdom]
    combined = cache / f"ubbotelorna_{kingdom}.hmm"

    if combined.exists():
        _log(f"  Using cached Rfam HMMs: {combined}")
        return combined

    hmmbuild = _require_tool("hmmbuild")

    _log(f"  Building Rfam HMMs for kingdom '{kingdom}' …")
    hmm_parts = []
    for name, acc in models.items():
        hmm_path = cache / f"{acc}.hmm"
        if not hmm_path.exists():
            sto_path = cache / f"{acc}.sto"
            if not sto_path.exists():
                url = _RFAM_STO_URL.format(acc=acc)
                _log(f"    Downloading {acc} ({name}) seed alignment from Rfam …")
                try:
                    urlretrieve(url, sto_path)
                except Exception as exc:
                    print(f"ERROR: failed to download {url}: {exc}\n"
                          f"       Check internet access or provide --rfam_dir with "
                          f"pre-built .hmm files.", file=sys.stderr)
                    sys.exit(1)
            _log(f"    Building HMM for {acc} ({name}) …")
            _run([hmmbuild, "--rna", "--cpu", "1",
                  str(hmm_path), str(sto_path)],
                 capture_stdout=True)
        hmm_parts.append(hmm_path.read_text())

    combined.write_text("\n".join(hmm_parts))
    _log(f"  HMMs written to {combined}")
    return combined


def _parse_nhmmer_tblout(tblout: Path, evalue: float) -> list[dict]:
    """Parse nhmmer --tblout into a list of hit dicts, filtering by evalue."""
    hits = []
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 15:
                continue
            try:
                hit_evalue = float(cols[12])
            except ValueError:
                continue
            if hit_evalue > evalue:
                continue
            start  = int(cols[6])
            end    = int(cols[7])
            strand = cols[11]
            if start > end:
                start, end = end, start
            hits.append({
                "seqname":    cols[0],
                "model_name": cols[2],
                "model_acc":  cols[3],
                "start":      start,
                "end":        end,
                "strand":     strand,
                "score":      float(cols[13]),
                "evalue":     hit_evalue,
            })
    return hits


def _parse_cmsearch_tblout(tblout: Path, evalue: float) -> list[dict]:
    """Parse cmsearch --tblout output into a list of hit dicts."""
    hits = []
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 17:
                continue
            try:
                hit_evalue = float(cols[15])
            except ValueError:
                continue
            if hit_evalue > evalue:
                continue
            seq_from = int(cols[7])
            seq_to   = int(cols[8])
            strand   = cols[9]
            # Ensure start < end regardless of strand
            if seq_from > seq_to:
                seq_from, seq_to = seq_to, seq_from
            hits.append({
                "seqname":    cols[0],
                "model_name": cols[2],
                "model_acc":  cols[3],
                "start":      seq_from,
                "end":        seq_to,
                "strand":     strand,
                "score":      float(cols[14]),
                "evalue":     hit_evalue,
                "inc":        cols[16],
            })
    return hits


def _run_single_chunk(chunk_records: list, chunk_idx: int,
                      tool: str, search_tool: str, model_file: Path,
                      threads: int, evalue: float, mxsize: float,
                      workdir: Path) -> list:
    """Write a temporary chunk FASTA, run nhmmer/cmsearch, return parsed hits.

    Temporary files are removed after parsing so workdir does not accumulate
    per-chunk debris.  The chunk index is zero-padded to four digits so that
    log messages remain sortable.
    """
    chunk_fa  = workdir / f"_chunk_{chunk_idx:04d}.fasta"
    chunk_tbl = workdir / f"_chunk_{chunk_idx:04d}.tblout"
    _write_fasta(chunk_records, chunk_fa)

    if search_tool == "nhmmer":
        _run([
            tool,
            "--rna",
            "--cpu",    str(threads),
            "--tblout", str(chunk_tbl),
            "-E",       str(evalue),
            "--noali",
            str(model_file),
            str(chunk_fa.resolve()),
        ], cwd=workdir)
        hits = _parse_nhmmer_tblout(chunk_tbl, evalue)
    else:
        _run([
            tool,
            "--cpu",    str(threads),
            "--tblout", str(chunk_tbl),
            "-E",       str(evalue),
            "--noali",
            "--rfam",
            "--mxsize", str(mxsize),
            str(model_file),
            str(chunk_fa.resolve()),
        ], cwd=workdir)
        hits = _parse_cmsearch_tblout(chunk_tbl, evalue)

    chunk_fa.unlink(missing_ok=True)
    chunk_tbl.unlink(missing_ok=True)
    return hits


def run_module2_rrna(fasta: Path, kingdom: str, threads: int, evalue: float,
                     rfam_dir: Path | None, mxsize: float, search_tool: str,
                     workdir: Path, results: Path,
                     prefix: str, force: bool, chunk_size: int = 200,
                     chunk_max_bp: int = 50_000_000,
                     flag_5s_arrays: bool = True, min_5s_copies: int = 5,
                     max_5s_gap: int = 5000, cap_5s: int = 0) -> Path:
    """Run nhmmer (default) or cmsearch and write rRNA GFF3."""
    out_gff3  = results / f"mod02_rRNA_{prefix}.gff3"
    tblout    = workdir / f"{search_tool}_rRNA.tblout"

    if _checkpoint(out_gff3, f"{search_tool}-rRNA", force):
        if (flag_5s_arrays or cap_5s > 0) and out_gff3.exists() and out_gff3.stat().st_size > 0:
            _filter_5s_arrays(out_gff3, prefix, results, min_5s_copies, max_5s_gap, cap_5s, force)
        return out_gff3

    search_fa = _hard_masked(workdir)
    if not search_fa.exists():
        search_fa = fasta   # fall back if masking was skipped

    if search_tool == "nhmmer":
        hmm_file = _ensure_hmms(kingdom, rfam_dir)
        tool     = _require_tool("nhmmer")
    else:
        cm_file  = _ensure_cms(kingdom, rfam_dir)
        tool     = _require_tool("cmsearch")
    model_file = hmm_file if search_tool == "nhmmer" else cm_file

    if chunk_size > 0:
        # ── Chunked execution: stream genome, run tool per batch ──────────────
        # Chunks flush on whichever limit is hit first: sequence count
        # (chunk_size) or cumulative bp (chunk_max_bp). The bp cap matters
        # for chromosome-scale assemblies with few, very large sequences —
        # those would otherwise never reach chunk_size and the whole genome
        # would land in one chunk, defeating the point of chunking.
        bp_desc = f", max {chunk_max_bp:,} bp/chunk" if chunk_max_bp > 0 else ""
        _log(f"  Chunked mode: {chunk_size} sequences per chunk{bp_desc}")
        all_hits:     list = []
        chunk_buf:    list = []
        chunk_idx          = 0
        total_seqs         = 0
        chunk_bp           = 0

        def _chunk_full() -> bool:
            if chunk_size > 0 and len(chunk_buf) >= chunk_size:
                return True
            if chunk_max_bp > 0 and chunk_bp >= chunk_max_bp:
                return True
            return False

        for name, seq in _iter_fasta(search_fa):
            chunk_buf.append((name, seq))
            chunk_bp   += len(seq)
            total_seqs += 1
            if _chunk_full():
                chunk_idx += 1
                _log(f"  Processing chunk {chunk_idx} "
                     f"({len(chunk_buf)} seqs, {chunk_bp:,} bp, "
                     f"total so far: {total_seqs})")
                all_hits.extend(
                    _run_single_chunk(chunk_buf, chunk_idx, tool, search_tool,
                                      model_file, threads, evalue, mxsize, workdir)
                )
                chunk_buf = []
                chunk_bp  = 0

        if chunk_buf:
            chunk_idx += 1
            _log(f"  Processing chunk {chunk_idx} "
                 f"({len(chunk_buf)} seqs, {chunk_bp:,} bp, total: {total_seqs})")
            all_hits.extend(
                _run_single_chunk(chunk_buf, chunk_idx, tool, search_tool,
                                  model_file, threads, evalue, mxsize, workdir)
            )

        _log(f"  Processed {total_seqs} sequences in {chunk_idx} chunk(s)")
        hits = all_hits
    else:
        # ── Single-run execution (original behaviour) ─────────────────────────
        if search_tool == "nhmmer":
            _run([
                tool,
                "--rna",                # force RNA alphabet — avoids 'Unable to guess alphabet'
                "--cpu",    str(threads),  # on N-heavy hard-masked FASTA; RNA is compatible with hmmbuild --rna HMMs
                "--tblout", str(tblout),
                "-E",       str(evalue),
                "--noali",
                str(model_file),
                str(search_fa.resolve()),
            ], cwd=workdir)
            hits = _parse_nhmmer_tblout(tblout, evalue)
        else:
            _run([
                tool,
                "--cpu",    str(threads),
                "--tblout", str(tblout),
                "-E",       str(evalue),
                "--noali",
                "--rfam",                   # HMM pre-filter; essential for large genomes
                "--mxsize", str(mxsize),    # cap DP matrix per thread (Mb)
                str(model_file),
                str(search_fa.resolve()),
            ], cwd=workdir)
            hits = _parse_cmsearch_tblout(tblout, evalue)

    _log(f"  {search_tool} hits (E ≤ {evalue}): {len(hits)}")
    hits = _sort_hits(hits)

    with open(out_gff3, "w") as fh:
        fh.write(_GFF3_HEADER)
        for i, h in enumerate(hits, 1):
            attrs = {
                "ID":       f"rRNA_{i}",
                "Name":     h["model_name"],
                "model":    h["model_acc"],
                "score":    f"{h['score']:.2f}",
                "E-value":  f"{h['evalue']:.2e}",
            }
            fh.write(_gff3_record(
                h["seqname"], "UbboTELORNA", "rRNA",
                h["start"], h["end"], h["score"], h["strand"], ".", attrs
            ) + "\n")

    _log(f"  rRNA GFF3: {out_gff3.name}  ({len(hits)} features)")
    if (flag_5s_arrays or cap_5s > 0) and out_gff3.exists() and out_gff3.stat().st_size > 0:
        _filter_5s_arrays(out_gff3, prefix, results, min_5s_copies, max_5s_gap, cap_5s, force)
    return out_gff3


# ── 5S rRNA post-filter ───────────────────────────────────────────────────────

def _filter_5s_arrays(gff3_path: Path, prefix: str, results_dir: Path,
                      min_copies: int, max_gap: int, cap: int, force: bool) -> None:
    """Flag 5S rRNA tandem array members in GFF3 attributes and optionally cap per-sequence count."""
    array_tsv = results_dir / f"mod02_5s_arrays_{prefix}.tsv"
    if _checkpoint(array_tsv, "5S-array-flagging", force):
        return

    _log("  5S rRNA post-filter: reading GFF3 …")
    header_lines: list[str] = []
    feature_rows: list[list[str]] = []
    with open(gff3_path) as fh:
        for line in fh:
            stripped = line.rstrip("\n")
            if stripped.startswith("#") or not stripped.strip():
                header_lines.append(stripped)
                continue
            cols = stripped.split("\t")
            if len(cols) < 9:
                header_lines.append(stripped)
                continue
            # Strip any prior array attributes for idempotency on --force re-runs
            attrs = re.sub(r";?array_member=[^;]+", "", cols[8])
            attrs = re.sub(r";?array_id=[^;]+", "", attrs)
            cols[8] = attrs
            feature_rows.append(cols)

    # Collect 5S hits with their index in feature_rows
    fives_hits: list[dict] = []
    for i, cols in enumerate(feature_rows):
        name_m = re.search(r"Name=([^;]+)", cols[8])
        name = name_m.group(1).strip() if name_m else ""
        if name != "5S_rRNA":
            continue
        try:
            start = int(cols[3])
            end   = int(cols[4])
        except ValueError:
            continue
        score_m = re.search(r"score=([^;]+)", cols[8])
        score = float(score_m.group(1)) if score_m else 0.0
        fives_hits.append({
            "idx":     i,
            "seqname": cols[0],
            "start":   start,
            "end":     end,
            "score":   score,
        })

    _log(f"  5S rRNA post-filter: {len(fives_hits)} predictions found")

    # Group by sequence
    by_seq: dict[str, list] = {}
    for h in fives_hits:
        by_seq.setdefault(h["seqname"], []).append(h)

    # Detect arrays and build member map
    array_records: list[dict] = []
    member_to_array: dict[int, str] = {}
    arr_counter = 0

    for seqname in sorted(by_seq.keys(), key=_natural_key):
        seq_hits = sorted(by_seq[seqname], key=lambda x: x["start"])
        cluster: list[dict] = [seq_hits[0]]
        for h in seq_hits[1:]:
            if h["start"] - cluster[-1]["end"] <= max_gap:
                cluster.append(h)
            else:
                if len(cluster) >= min_copies:
                    arr_id = f"arr{arr_counter:04d}"
                    for m in cluster:
                        member_to_array[m["idx"]] = arr_id
                    array_records.append({
                        "array_id": arr_id,
                        "seqname":  seqname,
                        "start":    min(m["start"] for m in cluster),
                        "end":      max(m["end"] for m in cluster),
                        "n_copies": len(cluster),
                        "span_bp":  max(m["end"] for m in cluster) - min(m["start"] for m in cluster) + 1,
                    })
                    arr_counter += 1
                cluster = [h]
        if len(cluster) >= min_copies:
            arr_id = f"arr{arr_counter:04d}"
            for m in cluster:
                member_to_array[m["idx"]] = arr_id
            array_records.append({
                "array_id": arr_id,
                "seqname":  seqname,
                "start":    min(m["start"] for m in cluster),
                "end":      max(m["end"] for m in cluster),
                "n_copies": len(cluster),
                "span_bp":  max(m["end"] for m in cluster) - min(m["start"] for m in cluster) + 1,
            })
            arr_counter += 1

    _log(f"  5S arrays detected: {arr_counter} arrays, "
         f"{len(member_to_array)} array members flagged")

    # Cap: retain highest-scoring hits per sequence, mark the rest for removal
    removed_idxs: set[int] = set()
    if cap > 0:
        for seqname, seq_hits in by_seq.items():
            if len(seq_hits) > cap:
                kept = {h["idx"] for h in sorted(seq_hits, key=lambda x: -x["score"])[:cap]}
                for h in seq_hits:
                    if h["idx"] not in kept:
                        removed_idxs.add(h["idx"])
        if removed_idxs:
            _log(f"  --cap_5s {cap}: removing {len(removed_idxs)} low-scoring 5S "
                 f"predictions (keeping top {cap} per sequence)")

    # Rebuild GFF3: annotate array members, drop capped hits
    new_rows: list[list[str]] = []
    for i, cols in enumerate(feature_rows):
        if i in removed_idxs:
            continue
        if i in member_to_array:
            arr_id = member_to_array[i]
            cols[8] = cols[8].rstrip(";") + f";array_member=true;array_id={arr_id}"
        new_rows.append(cols)

    with open(gff3_path, "w") as fh:
        for h in header_lines:
            fh.write(h + "\n")
        for cols in new_rows:
            fh.write("\t".join(cols) + "\n")

    _log(f"  GFF3 updated: {len(removed_idxs)} removed, "
         f"{len(member_to_array)} array members annotated  →  {gff3_path.name}")

    # Write array summary TSV sorted by seqname, start
    array_records.sort(key=lambda a: (_natural_key(a["seqname"]), a["start"]))
    with open(array_tsv, "w") as fh:
        fh.write("array_id\tseqname\tstart\tend\tstrand\tn_copies\tspan_bp\n")
        for a in array_records:
            fh.write(
                f"{a['array_id']}\t{a['seqname']}\t{a['start']}\t"
                f"{a['end']}\t.\t{a['n_copies']}\t{a['span_bp']}\n"
            )
    _log(f"  5S array TSV: {array_tsv.name}  ({len(array_records)} arrays)")


# ── Module 3: tRNA annotation (ARAGORN) ───────────────────────────────────────

_ARAGORN_RE = re.compile(
    r"\s*\d+\s+((?:tRNA|tmRNA|mtRNA|pseudo_tRNA)-\S+)\s+"
    r"(c?)\[(\d+),(\d+)\](?:\s+\(([^)]+)\))?"
)


def _parse_aragorn(aragorn_out: Path) -> list[dict]:
    """Parse ARAGORN text output into a list of hit dicts."""
    hits    = []
    current = None
    with open(aragorn_out) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].split()[0]
                continue
            if current is None:
                continue
            m = _ARAGORN_RE.match(line)
            if m:
                trna_type  = m.group(1)
                complement = m.group(2) == "c"
                start      = int(m.group(3))
                end        = int(m.group(4))
                anticodon  = m.group(5) or ""
                if start > end:
                    start, end = end, start
                hits.append({
                    "seqname":   current,
                    "trna_type": trna_type,
                    "start":     start,
                    "end":       end,
                    "strand":    "-" if complement else "+",
                    "anticodon": anticodon,
                })
    return hits


def run_module3_trna(fasta: Path, workdir: Path, results: Path,
                     prefix: str, force: bool) -> Path:
    """Run ARAGORN and write tRNA GFF3."""
    out_gff3   = results / f"mod03_tRNA_{prefix}.gff3"
    aragorn_out = workdir / "aragorn.txt"

    if _checkpoint(out_gff3, "ARAGORN-tRNA", force):
        return out_gff3

    aragorn = _require_tool("aragorn")
    # Use the original (unmasked) FASTA — ARAGORN is a structural predictor
    # that uses the tRNA cloverleaf fold, not an alignment-based search.
    # The hard-masked FASTA (N's replacing repetitive regions) destroys the
    # structural signal and causes ARAGORN to miss all tRNAs.
    _run([
        aragorn,
        "-t",           # tRNA genes only
        "-gcstd",       # standard genetic code
        "-l",           # treat each sequence as a linear molecule
        "-w",           # write sequence headers as >name (required for parser)
        "-o", str(aragorn_out),
        str(fasta.resolve()),
    ])

    hits = _sort_hits(_parse_aragorn(aragorn_out))
    _log(f"  ARAGORN tRNA hits: {len(hits)}")

    with open(out_gff3, "w") as fh:
        fh.write(_GFF3_HEADER)
        for i, h in enumerate(hits, 1):
            attrs = {
                "ID":        f"tRNA_{i}",
                "Name":      h["trna_type"],
                "anticodon": h["anticodon"],
            }
            fh.write(_gff3_record(
                h["seqname"], "UbboTELORNA", "tRNA",
                h["start"], h["end"], ".", h["strand"], ".", attrs
            ) + "\n")

    _log(f"  tRNA GFF3: {out_gff3.name}  ({len(hits)} features)")
    return out_gff3


# ── Module 4: Integration ─────────────────────────────────────────────────────

def _collect_and_write_gff3(tel_gff: Path | None, rrna_gff: Path | None,
                             trna_gff: Path | None,
                             combined_path: Path) -> tuple[list, list, list]:
    """Read three GFF3 files, sort by position, write combined; return (tel, rrna, trna) lines."""
    all_lines: list[str] = []
    tel_lines:  list[str] = []
    rrna_lines: list[str] = []
    trna_lines: list[str] = []
    for gff, bucket in [
        (tel_gff,  tel_lines),
        (rrna_gff, rrna_lines),
        (trna_gff, trna_lines),
    ]:
        if gff is None or not gff.exists():
            continue
        lines = [l for l in open(gff) if not l.startswith("#") and l.strip()]
        bucket.extend(lines)
        all_lines.extend(lines)
    all_lines.sort(key=_gff3_sort_key)
    with open(combined_path, "w") as out:
        out.write(_GFF3_HEADER)
        for line in all_lines:
            out.write(line if line.endswith("\n") else line + "\n")
    _log(f"  Combined GFF3: {combined_path.name}")
    return tel_lines, rrna_lines, trna_lines


def run_module4_integration(tel_gff: Path | None, rrna_gff: Path | None,
                             trna_gff: Path | None,
                             results: Path, prefix: str,
                             genome_size: int = 0) -> tuple[Path, dict]:
    """Merge GFF3 files, sort by position, write summary TSV; return (combined_path, counts_dict)."""
    combined = results / f"mod04_annotation_{prefix}.gff3"
    summary  = results / f"mod04_summary_{prefix}.tsv"

    tel_lines, rrna_lines, trna_lines = _collect_and_write_gff3(
        tel_gff, rrna_gff, trna_gff, combined)

    # ── Count subtypes and accumulate lengths ──────────────────────────────────
    n_tel  = len(tel_lines)
    n_rrna = len(rrna_lines)
    n_trna = len(trna_lines)

    tel_by_end    = Counter()   # subtype  → count
    rrna_by_type  = Counter()
    trna_by_type  = Counter()
    tel_len_by_end   = Counter()  # subtype → total bp
    rrna_len_by_type = Counter()
    trna_len_by_type = Counter()
    tel_len = rrna_len = trna_len = 0

    for line in tel_lines:
        name = _gff3_name(line)
        end  = name.replace("telomere_", "") if name else "unknown"
        bp   = _gff3_length(line)
        tel_by_end[end]    += 1
        tel_len_by_end[end] += bp
        tel_len += bp

    for line in rrna_lines:
        rtype = _gff3_name(line) or "unknown"
        bp    = _gff3_length(line)
        rrna_by_type[rtype]   += 1
        rrna_len_by_type[rtype] += bp
        rrna_len += bp

    for line in trna_lines:
        ttype = _gff3_name(line) or "unknown"
        bp    = _gff3_length(line)
        trna_by_type[ttype]   += 1
        trna_len_by_type[ttype] += bp
        trna_len += bp

    def _pct(bp: int) -> str:
        if genome_size > 0:
            return f"{bp / genome_size * 100:.4f}"
        return "NA"

    # ── Write detailed summary TSV ─────────────────────────────────────────────
    with open(summary, "w") as fh:
        fh.write("feature_type\tsubtype\tcount\ttotal_length_bp\tpct_genome\n")

        # Telomeres
        fh.write(f"telomere\tTOTAL\t{n_tel}\t{tel_len}\t{_pct(tel_len)}\n")
        for end in sorted(tel_by_end):
            bp = tel_len_by_end[end]
            fh.write(f"telomere\t{end}\t{tel_by_end[end]}\t{bp}\t{_pct(bp)}\n")

        # rRNA — biological order (small → large subunit)
        fh.write(f"rRNA\tTOTAL\t{n_rrna}\t{rrna_len}\t{_pct(rrna_len)}\n")
        for rtype in _RRNA_ORDER:
            if rtype in rrna_by_type:
                bp = rrna_len_by_type[rtype]
                fh.write(f"rRNA\t{rtype}\t{rrna_by_type[rtype]}\t{bp}\t{_pct(bp)}\n")
        for rtype in sorted(rrna_by_type):
            if rtype not in _RRNA_ORDER:
                bp = rrna_len_by_type[rtype]
                fh.write(f"rRNA\t{rtype}\t{rrna_by_type[rtype]}\t{bp}\t{_pct(bp)}\n")

        # tRNA — sorted by count descending, then alphabetically
        fh.write(f"tRNA\tTOTAL\t{n_trna}\t{trna_len}\t{_pct(trna_len)}\n")
        for ttype, cnt in sorted(trna_by_type.items(),
                                  key=lambda x: (-x[1], x[0])):
            bp = trna_len_by_type[ttype]
            fh.write(f"tRNA\t{ttype}\t{cnt}\t{bp}\t{_pct(bp)}\n")

    _log(f"  Summary TSV: {summary.name}")

    # ── Log breakdown ──────────────────────────────────────────────────────────
    _log(f"  Totals  →  telomeres: {n_tel}  |  rRNA: {n_rrna}  |  tRNA: {n_trna}")
    if rrna_by_type:
        _log("  rRNA breakdown:")
        for rtype in _RRNA_ORDER:
            if rtype in rrna_by_type:
                _log(f"    {rtype:30s} {rrna_by_type[rtype]:>6}")
    if trna_by_type:
        top = sorted(trna_by_type.items(), key=lambda x: -x[1])[:5]
        _log("  tRNA top types (count): " +
             "  ".join(f"{t}:{c}" for t, c in top))

    counts = {
        "n_tel":             n_tel,
        "n_rrna":            n_rrna,
        "n_trna":            n_trna,
        "tel_len":           tel_len,
        "rrna_len":          rrna_len,
        "trna_len":          trna_len,
        "tel_by_end":        dict(tel_by_end),
        "tel_len_by_end":    dict(tel_len_by_end),
        "rrna_by_type":      dict(rrna_by_type),
        "rrna_len_by_type":  dict(rrna_len_by_type),
        "trna_by_type":      dict(trna_by_type),
        "trna_len_by_type":  dict(trna_len_by_type),
    }
    return combined, counts


# ── Module 5: Visualization ───────────────────────────────────────────────────

_C_TEL   = "#F5A623"   # amber  — telomere
_C_RRNA  = "#4C9BE8"   # blue   — rRNA
_C_TRNA  = "#E8604C"   # red    — tRNA
_C_OTHER = "#888888"   # grey   — uncharacterised / sequence bar

_RRNA_ORDER = [
    "5S_rRNA", "5_8S_rRNA",
    "SSU_rRNA_eukarya", "SSU_rRNA_bacteria", "SSU_rRNA_archaea",
    "LSU_rRNA_eukarya", "LSU_rRNA_bacteria", "LSU_rRNA_archaea",
]


def _natural_key(s: str) -> list:
    """Natural sort key: splits on digit runs so Chr2 < Chr10."""
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", s)]


def _draw_ideogram(ax: plt.Axes, top_seqs: list[tuple[str, int]],
                   tel_pos: dict, rrna_pos: dict, trna_pos: dict,
                   max_len: int, n_seqs: int, seq_lengths: dict,
                   prefix: str, sort_by: str) -> None:
    """Draw genome ideogram panel showing feature positions per sequence."""
    bar_h   = 0.65
    tel_h   = bar_h + 0.30
    min_vis = max_len * 0.002
    tel_min = max_len * 0.005

    n_rrna_total = sum(len(v) for v in rrna_pos.values())
    n_trna_total = sum(len(v) for v in trna_pos.values())
    n_tel_total  = sum(len(v) for v in tel_pos.values())

    layers = sorted([
        ("rrna", n_rrna_total, _C_RRNA, 0.40, bar_h),
        ("trna", n_trna_total, _C_TRNA, 0.60, bar_h),
        ("tel",  n_tel_total,  _C_TEL,  1.00, tel_h),
    ], key=lambda x: -x[1])

    for i, (name, slen) in enumerate(top_seqs):
        y = n_seqs - i - 1
        ax.broken_barh([(0, slen)], (y - bar_h / 2, bar_h),
                       facecolors=_C_OTHER, alpha=0.18, linewidth=0)
        for layer_name, _, color, alpha, height in layers:
            if layer_name == "rrna":
                pos_list = rrna_pos.get(name, [])
                mv = min_vis
            elif layer_name == "trna":
                pos_list = trna_pos.get(name, [])
                mv = min_vis
            else:
                pos_list = tel_pos.get(name, [])
                mv = tel_min
            segs = [(s - 1, max(e - s + 1, mv)) for s, e in pos_list]
            if segs:
                ax.broken_barh(segs, (y - height / 2, height),
                               facecolors=color, alpha=alpha, linewidth=0)

    ax.set_ylim(-0.8, n_seqs - 0.2)
    ax.set_yticks(range(n_seqs))
    ax.set_yticklabels([n for n, _ in reversed(top_seqs)], fontsize=8)
    ax.set_xlim(0, max_len * 1.01)
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x / 1e6:.0f} Mb"))
    ax.set_xlabel("Genomic position")
    sort_label = "by SeqID" if sort_by == "seqid" else "longest first"
    title_suffix = (f"top {n_seqs} of {len(seq_lengths)} sequences, {sort_label}"
                    if n_seqs < len(seq_lengths)
                    else f"{n_seqs} sequences, {sort_label}")
    ax.set_title(f"{prefix}  —  genome annotation overview  ({title_suffix})",
                 fontsize=11, pad=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(left=False)
    legend_patches = [
        mpatches.Patch(color=_C_TEL,  label=f"Telomere (n={n_tel_total:,})"),
        mpatches.Patch(color=_C_RRNA, alpha=0.40, label=f"rRNA (n={n_rrna_total:,})"),
        mpatches.Patch(color=_C_TRNA, alpha=0.60, label=f"tRNA (n={n_trna_total:,})"),
    ]
    ax.legend(handles=legend_patches, loc="lower right",
              frameon=True, framealpha=0.85, fontsize=9)


def _draw_rrna_bars(ax: plt.Axes, summ: dict) -> None:
    """Draw horizontal bar chart of rRNA subtype counts."""
    rrna_data = summ.get("rRNA", {})
    rrna_sub  = [(rt, rrna_data[rt][0]) for rt in _RRNA_ORDER if rt in rrna_data]
    if not rrna_sub:
        rrna_sub = sorted([(k, v[0]) for k, v in rrna_data.items() if k != "TOTAL"],
                          key=lambda x: -x[1])
    if rrna_sub:
        lbls, vals = zip(*rrna_sub)
        yp = list(range(len(lbls)))
        ax.barh(yp, vals, color=_C_RRNA, linewidth=0)
        ax.set_yticks(yp)
        ax.set_yticklabels(lbls, fontsize=9)
        ax.set_xlabel("Count")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        for j, v in enumerate(vals):
            ax.text(v + max(vals) * 0.01, j, f"{v:,}", va="center", fontsize=8)
        ax.set_xlim(right=max(vals) * 1.18)
    else:
        ax.text(0.5, 0.5, "No rRNA data", ha="center", va="center",
                transform=ax.transAxes, color=_C_OTHER)
    ax.set_title("rRNA subtypes")
    ax.spines[["top", "right"]].set_visible(False)


def _draw_trna_bars(ax: plt.Axes, summ: dict) -> None:
    """Draw horizontal bar chart of tRNA type counts (top 20)."""
    trna_data = summ.get("tRNA", {})
    trna_sub  = sorted([(k, v[0]) for k, v in trna_data.items() if k != "TOTAL"],
                       key=lambda x: (-x[1], x[0]))[:20]
    if trna_sub:
        lbls, vals = zip(*trna_sub)
        yp = list(range(len(lbls)))
        ax.barh(yp, vals, color=_C_TRNA, linewidth=0)
        ax.set_yticks(yp)
        ax.set_yticklabels(lbls, fontsize=9)
        ax.set_xlabel("Count")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        for j, v in enumerate(vals):
            ax.text(v + max(vals) * 0.01, j, f"{v:,}", va="center", fontsize=8)
        ax.set_xlim(right=max(vals) * 1.18)
    else:
        ax.text(0.5, 0.5, "No tRNA data", ha="center", va="center",
                transform=ax.transAxes, color=_C_OTHER)
    ax.set_title("tRNA types (top 20)")
    ax.spines[["top", "right"]].set_visible(False)


def _draw_donut(ax: plt.Axes, summ: dict, genome_size: int) -> None:
    """Draw genome composition donut chart."""
    tel_bp   = summ.get("telomere", {}).get("TOTAL", (0, 0))[1]
    rrna_bp  = summ.get("rRNA",     {}).get("TOTAL", (0, 0))[1]
    trna_bp  = summ.get("tRNA",     {}).get("TOTAL", (0, 0))[1]
    other_bp = max(0, genome_size - tel_bp - rrna_bp - trna_bp)
    slices = [(v, l, c) for v, l, c in [
        (tel_bp,   "Telomere", _C_TEL),
        (rrna_bp,  "rRNA",     _C_RRNA),
        (trna_bp,  "tRNA",     _C_TRNA),
        (other_bp, "Other",    _C_OTHER),
    ] if v > 0]
    if slices and genome_size > 0:
        vals, lbls, cols = zip(*slices)
        wedges, _ = ax.pie(
            vals, colors=cols,
            wedgeprops=dict(width=0.45, edgecolor="white", linewidth=1.5),
            startangle=90,
        )
        ax.text(0, 0, f"{genome_size / 1e6:.0f} Mb",
                ha="center", va="center", fontsize=10, fontweight="bold")
        ax.legend(
            wedges,
            [f"{l}  {v / genome_size * 100:.3f}%" for v, l in zip(vals, lbls)],
            loc="lower center", bbox_to_anchor=(0.5, -0.18),
            frameon=False, fontsize=8, ncol=2,
        )
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color=_C_OTHER)
    ax.set_title("Genome composition")


def run_module5_plot(fasta: Path,
                     tel_gff: Path | None, rrna_gff: Path | None,
                     trna_gff: Path | None, summary_tsv: Path | None,
                     genome_size: int, results: Path, prefix: str,
                     plot_formats: list, top_sequences: int,
                     sort_by: str, force: bool) -> None:
    """Generate 3-panel figure: ideogram, subtype bars, composition donut."""
    out_base  = results / f"mod05_plot_{prefix}"
    out_paths = [out_base.with_suffix(f".{fmt}") for fmt in plot_formats]
    if _checkpoint(out_paths[0], "visualization", force):
        return

    seq_lengths = _fasta_seq_lengths(fasta)
    tel_pos     = _parse_gff3_positions(tel_gff)
    rrna_pos    = _parse_gff3_positions(rrna_gff)
    trna_pos    = _parse_gff3_positions(trna_gff)
    summ        = _parse_summary_tsv(summary_tsv)

    if sort_by == "seqid":
        top_seqs = sorted(seq_lengths.items(), key=lambda x: _natural_key(x[0]))[:top_sequences]
    else:
        top_seqs = sorted(seq_lengths.items(), key=lambda x: -x[1])[:top_sequences]

    if not top_seqs:
        _log("  [Module 5] No sequences found — skipping plot")
        return

    max_len = max(slen for _, slen in top_seqs)
    n_seqs  = len(top_seqs)
    _log(f"  Plotting {n_seqs} sequences, sorted by {sort_by} "
         f"(longest: {max_len:,} bp)")

    ideo_h = max(4, min(n_seqs * 0.28, 14))
    fig    = plt.figure(figsize=(18, ideo_h + 5))
    gs     = GridSpec(2, 3, figure=fig,
                      height_ratios=[ideo_h, 4.5],
                      hspace=0.45, wspace=0.38)
    ax_ideo  = fig.add_subplot(gs[0, :])
    ax_rrna  = fig.add_subplot(gs[1, 0])
    ax_trna  = fig.add_subplot(gs[1, 1])
    ax_donut = fig.add_subplot(gs[1, 2])

    _draw_ideogram(ax_ideo, top_seqs, tel_pos, rrna_pos, trna_pos,
                   max_len, n_seqs, seq_lengths, prefix, sort_by)
    _draw_rrna_bars(ax_rrna, summ)
    _draw_trna_bars(ax_trna, summ)
    _draw_donut(ax_donut, summ, genome_size)

    for fmt in plot_formats:
        out_path = out_base.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        _log(f"  Plot saved: {out_path.name}")
    plt.close(fig)


# ── Module 6: Evolutionary analysis ──────────────────────────────────────────

_RFAM_TRNA_ACC  = "RF00005"   # universal tRNA covariance model
_TRNA_CM_SCORE_THRESHOLD = 20.0   # bits — tRNAscan-SE's own pseudogene threshold
_TRNA_LEN_MIN   = 50          # bp
_TRNA_LEN_MAX   = 150         # bp
_RRNA_ARRAY_GAP = 50_000      # bp — max gap between consecutive rRNA hits in one array
_TRNA_ARRAY_GAP = 10_000      # bp — max gap between consecutive tRNA hits in one array

def _parse_gff3_records(gff: Path | None) -> list:
    """Parse a GFF3 into full feature dicts (score, subtype, length …)."""
    records = []
    if gff is None or not gff.exists():
        return records
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            try:
                start = int(cols[3])
                end   = int(cols[4])
            except ValueError:
                continue
            attrs    = cols[8]
            id_m     = re.search(r'ID=([^;]+)',      attrs)
            name_m   = re.search(r'Name=([^;]+)',    attrs)
            score_m  = re.search(r'score=([^;]+)',   attrs)
            evalue_m = re.search(r'E-value=([^;]+)', attrs)
            records.append({
                "seqname":  cols[0],
                "start":    start,
                "end":      end,
                "strand":   cols[6],
                "feat_id":  id_m.group(1)   if id_m    else f"{cols[0]}:{start}-{end}",
                "subtype":  name_m.group(1) if name_m  else "",
                "score":    float(score_m.group(1))  if score_m  else None,
                "evalue":   float(evalue_m.group(1)) if evalue_m else None,
                "length":   end - start + 1,
                "attrs":    attrs.rstrip(),
            })
    return records


def _extract_sequences(genome: Path, records: list) -> dict:
    """Stream genome FASTA and extract feature sequences. Peak memory = one chromosome."""
    by_seq: dict = {}
    for r in records:
        by_seq.setdefault(r["seqname"], []).append(r)

    result: dict = {}
    current_name: str | None  = None
    current_parts: list | None = None

    def _flush(name, parts):
        if name not in by_seq or parts is None:
            return
        seq = "".join(parts)
        for r in by_seq[name]:
            sub = seq[r["start"] - 1 : r["end"]]
            if r.get("strand") == "-":
                sub = _revcomp(sub)
            result[r["feat_id"]] = sub

    with open(genome) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                _flush(current_name, current_parts)
                current_name  = line[1:].split()[0]
                current_parts = [] if current_name in by_seq else None
            elif current_parts is not None:
                current_parts.append(line)
    _flush(current_name, current_parts)
    return result


# ── Array detection ───────────────────────────────────────────────────────────

def _summarise_array(arr_id: int, seqname: str, cluster: list) -> dict:
    cluster  = sorted(cluster, key=lambda x: x["start"])
    starts   = [f["start"] for f in cluster]
    ends     = [f["end"]   for f in cluster]
    spacings = [cluster[i]["start"] - cluster[i - 1]["end"]
                for i in range(1, len(cluster))]
    subtypes = Counter(f["subtype"] for f in cluster)

    ssu = (subtypes.get("SSU_rRNA_eukarya", 0) +
           subtypes.get("SSU_rRNA_bacteria", 0) +
           subtypes.get("SSU_rRNA_archaea",  0))
    lsu = (subtypes.get("LSU_rRNA_eukarya", 0) +
           subtypes.get("LSU_rRNA_bacteria", 0) +
           subtypes.get("LSU_rRNA_archaea",  0))
    s58 = subtypes.get("5_8S_rRNA", 0)
    complete_units = min(ssu, s58, lsu) if (ssu and s58 and lsu) else min(ssu, lsu)

    return {
        "array_id":          arr_id,
        "seqname":           seqname,
        "array_start":       min(starts),
        "array_end":         max(ends),
        "n_copies":          len(cluster),
        "array_length_bp":   max(ends) - min(starts) + 1,
        "subtypes":          dict(subtypes),
        "spacings":          spacings,
        "mean_spacing_bp":   round(sum(spacings) / len(spacings)) if spacings else 0,
        "min_spacing_bp":    min(spacings, default=0),
        "max_spacing_bp":    max(spacings, default=0),
        "complete_rDNA_units": complete_units,
    }


def _detect_arrays(records: list, max_gap: int, min_copies: int = 2) -> list:
    """Cluster GFF3 feature records into tandem arrays by proximity."""
    by_seq: dict = {}
    for r in records:
        by_seq.setdefault(r["seqname"], []).append(r)

    arrays = []
    arr_id = 0
    for seqname in sorted(by_seq.keys(), key=_natural_key):
        feats = sorted(by_seq[seqname], key=lambda x: x["start"])
        cluster = [feats[0]]
        for f in feats[1:]:
            if f["start"] - cluster[-1]["end"] <= max_gap:
                cluster.append(f)
            else:
                if len(cluster) >= min_copies:
                    arrays.append(_summarise_array(arr_id, seqname, cluster))
                    arr_id += 1
                cluster = [f]
        if len(cluster) >= min_copies:
            arrays.append(_summarise_array(arr_id, seqname, cluster))
            arr_id += 1
    return arrays


# ── tRNA pseudogene classification ────────────────────────────────────────────

def _ensure_trna_cm(rfam_dir: Path) -> Path:
    """Download or locate the universal tRNA covariance model (RF00005)."""
    cm_path = rfam_dir / f"{_RFAM_TRNA_ACC}.cm"
    if cm_path.exists() and cm_path.stat().st_size > 0:
        _log(f"  tRNA CM found: {cm_path}")
        return cm_path
    url = _RFAM_CM_URL.format(acc=_RFAM_TRNA_ACC)
    _log(f"  Downloading {_RFAM_TRNA_ACC} (tRNA CM) …")
    try:
        urlretrieve(url, cm_path)
    except Exception as e:
        print(f"ERROR: failed to download {_RFAM_TRNA_ACC}: {e}", file=sys.stderr)
        sys.exit(1)
    return cm_path


def _score_trnas_cmsearch(trna_seqs: dict, cm_path: Path,
                           workdir: Path, threads: int, force: bool) -> dict:
    """Run cmsearch RF00005 against extracted tRNA sequences; return {feat_id: best_bit_score}."""
    cmsearch  = _require_tool("cmsearch")
    trna_fa   = workdir / "mod06_trna_seqs.fasta"
    tblout    = workdir / "mod06_trna_cmsearch.tblout"

    with open(trna_fa, "w") as fh:
        for feat_id, seq in trna_seqs.items():
            if seq:
                fh.write(f">{feat_id}\n{seq}\n")

    if not _checkpoint(tblout, "tRNA cmsearch RF00005", force):
        _run([cmsearch, "--cpu", str(threads),
              "--tblout", str(tblout), "--noali",
              "-E", "1000",          # lenient — tRNA seqs are short
              "--incE", "1000",
              str(cm_path), str(trna_fa)],
             capture_stdout=True)

    scores: dict = {}
    with open(tblout) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 16:
                continue
            try:
                score = float(cols[14])
            except ValueError:
                continue
            feat_id = cols[0]
            if feat_id not in scores or score > scores[feat_id]:
                scores[feat_id] = score
    return scores


def _classify_trnas(records: list, cm_scores: dict) -> list:
    """Classify tRNA records as functional or pseudogene candidate."""
    result = []
    for r in records:
        cm_score = cm_scores.get(r["feat_id"])
        ac_m     = re.search(r'\(([^)]+)\)', r["subtype"])
        anticodon = ac_m.group(1) if ac_m else "???"

        reasons = []
        if cm_score is None:
            reasons.append("no CM hit")
        elif cm_score < _TRNA_CM_SCORE_THRESHOLD:
            reasons.append(f"CM score {cm_score:.1f} < {_TRNA_CM_SCORE_THRESHOLD}")
        if "???" in anticodon:
            reasons.append("unrecognised anticodon")
        if r["length"] < _TRNA_LEN_MIN or r["length"] > _TRNA_LEN_MAX:
            reasons.append(f"length {r['length']} bp out of [50,150]")

        result.append({
            **r,
            "cm_score":  cm_score,
            "anticodon": anticodon,
            "category":  "pseudogene" if reasons else "functional",
            "reason":    "; ".join(reasons),
        })
    return result


# ── Module 6 plots ────────────────────────────────────────────────────────────

def _plot_evolution(rrna_records: list, trna_classified: list,
                    rrna_arrays: list, trna_arrays: list,
                    out_base: Path, plot_formats: list, prefix: str) -> None:
    """Two-row figure: rRNA score distributions / tRNA pseudogenes + array stats."""

    # Collect rRNA scores per subtype (biological order)
    rrna_by_sub: dict = {}
    for r in rrna_records:
        if r["score"] is not None:
            rrna_by_sub.setdefault(r["subtype"], []).append(r["score"])
    sub_present = [s for s in _RRNA_ORDER if s in rrna_by_sub]
    for s in rrna_by_sub:
        if s not in sub_present:
            sub_present.append(s)
    n_sub = max(len(sub_present), 3)

    fig = plt.figure(figsize=(max(18, n_sub * 4.2), 11))
    gs  = GridSpec(2, n_sub, figure=fig, hspace=0.50, wspace=0.38)

    # ── Row 1: rRNA bit score histogram per subtype ────────────────────────────
    for i, subtype in enumerate(sub_present):
        ax = fig.add_subplot(gs[0, i])
        scores = rrna_by_sub[subtype]
        ax.hist(scores, bins=50, color=_C_RRNA, edgecolor="white", linewidth=0.3)
        short = subtype.replace("_rRNA", "").replace("_", " ")
        ax.set_title(short, fontsize=9)
        ax.set_xlabel("Bit score", fontsize=8)
        if i == 0:
            ax.set_ylabel("Count", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        med = sorted(scores)[len(scores) // 2]
        ax.axvline(med, color="black", linestyle="--", linewidth=1, alpha=0.6)
        ymax = ax.get_ylim()[1]
        ax.text(med * 1.01, ymax * 0.92,
                f"med\n{med:.0f}", fontsize=7, va="top", color="black", alpha=0.8)
        ax.annotate(f"n={len(scores):,}", xy=(0.97, 0.97),
                    xycoords="axes fraction", ha="right", va="top", fontsize=8)

    # ── Row 2a: tRNA CM score — functional vs pseudogene ──────────────────────
    ax_tc = fig.add_subplot(gs[1, 0:max(1, n_sub // 2)])
    if trna_classified:
        func_sc  = [t["cm_score"] for t in trna_classified
                    if t["category"] == "functional" and t["cm_score"] is not None]
        pseudo_sc = [t["cm_score"] for t in trna_classified
                     if t["category"] == "pseudogene" and t["cm_score"] is not None]
        all_sc = func_sc + pseudo_sc
        hi = max(all_sc) if all_sc else 100
        bins = [x * 5 for x in range(int(hi // 5) + 3)]
        if func_sc:
            ax_tc.hist(func_sc, bins=bins, color=_C_RRNA, alpha=0.75,
                       label=f"Functional  n={len(func_sc):,}")
        if pseudo_sc:
            ax_tc.hist(pseudo_sc, bins=bins, color=_C_TRNA, alpha=0.75,
                       label=f"Pseudogene candidate  n={len(pseudo_sc):,}")
        ax_tc.axvline(_TRNA_CM_SCORE_THRESHOLD, color="black",
                      linestyle="--", linewidth=1, alpha=0.7)
        ax_tc.text(_TRNA_CM_SCORE_THRESHOLD + 0.5, ax_tc.get_ylim()[1] * 0.92,
                   f"{_TRNA_CM_SCORE_THRESHOLD:.0f} bit\nthreshold",
                   fontsize=7, va="top")
        ax_tc.legend(fontsize=8, frameon=False)
    else:
        ax_tc.text(0.5, 0.5, "No tRNA data", ha="center", va="center",
                   transform=ax_tc.transAxes, color=_C_OTHER)
    ax_tc.set_title("tRNA structural score (cmsearch RF00005)\nfunctional vs pseudogene candidates",
                    fontsize=9)
    ax_tc.set_xlabel("Bit score", fontsize=8)
    ax_tc.set_ylabel("Count", fontsize=8)
    ax_tc.spines[["top", "right"]].set_visible(False)

    # ── Row 2b: rRNA inter-copy spacing distribution ───────────────────────────
    ax_sp = fig.add_subplot(gs[1, max(1, n_sub // 2) : max(2, n_sub * 3 // 4)])
    all_spacings = []
    for a in rrna_arrays:
        all_spacings.extend(a["spacings"])
    if all_spacings:
        sp_kb = [s / 1_000 for s in all_spacings if s >= 0]
        ax_sp.hist(sp_kb, bins=40, color=_C_RRNA, edgecolor="white", linewidth=0.3)
        med_sp = sorted(sp_kb)[len(sp_kb) // 2]
        ax_sp.axvline(med_sp, color="black", linestyle="--", linewidth=1, alpha=0.6)
        ax_sp.text(med_sp * 1.02, ax_sp.get_ylim()[1] * 0.92,
                   f"med\n{med_sp:.1f} kb", fontsize=7, va="top")
        ax_sp.annotate(f"n={len(sp_kb):,} pairs\n{len(rrna_arrays)} arrays",
                       xy=(0.97, 0.97), xycoords="axes fraction",
                       ha="right", va="top", fontsize=8)
    else:
        ax_sp.text(0.5, 0.5, "No rRNA arrays", ha="center", va="center",
                   transform=ax_sp.transAxes, color=_C_OTHER)
    ax_sp.set_title("rRNA inter-copy spacing\n(within arrays)", fontsize=9)
    ax_sp.set_xlabel("Spacing (kb)", fontsize=8)
    ax_sp.set_ylabel("Count", fontsize=8)
    ax_sp.spines[["top", "right"]].set_visible(False)

    # ── Row 2c: Array size distribution ───────────────────────────────────────
    ax_sz = fig.add_subplot(gs[1, max(2, n_sub * 3 // 4):])
    r_sizes = [a["n_copies"] for a in rrna_arrays]
    t_sizes = [a["n_copies"] for a in trna_arrays]
    hi_sz = max((r_sizes or [0]) + (t_sizes or [0]) + [2])
    bins_sz = list(range(2, hi_sz + 2))
    if r_sizes:
        ax_sz.hist(r_sizes, bins=bins_sz, color=_C_RRNA, alpha=0.75,
                   label=f"rRNA  n={len(r_sizes)} arrays")
    if t_sizes:
        ax_sz.hist(t_sizes, bins=bins_sz, color=_C_TRNA, alpha=0.75,
                   label=f"tRNA  n={len(t_sizes)} arrays")
    ax_sz.set_title("Tandem array size\ndistribution", fontsize=9)
    ax_sz.set_xlabel("Copies per array", fontsize=8)
    ax_sz.set_ylabel("Number of arrays", fontsize=8)
    ax_sz.legend(fontsize=8, frameon=False)
    ax_sz.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"{prefix}  —  evolutionary analysis", fontsize=11, y=1.01)

    for fmt in plot_formats:
        out_path = out_base.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        _log(f"  Plot saved: {out_path.name}")
    plt.close(fig)


# ── Module 6 sub-step helpers ─────────────────────────────────────────────────

def _run_rrna_score_analysis(rrna_records: list, results: Path,
                              prefix: str, force: bool) -> None:
    """Log rRNA bit score statistics and write score TSV."""
    if not rrna_records:
        _log("  No rRNA records found — skipping 6a")
        return
    rrna_score_tsv = results / f"mod06_rrna_scores_{prefix}.tsv"
    if force or not rrna_score_tsv.exists():
        sorted_rrna = sorted(rrna_records,
                             key=lambda r: (_natural_key(r["seqname"]), r["start"]))
        with open(rrna_score_tsv, "w") as fh:
            fh.write("feature_id\tseqname\tstart\tend\tstrand\t"
                     "subtype\tlength_bp\tbit_score\tevalue\n")
            for r in sorted_rrna:
                fh.write(
                    f"{r['feat_id']}\t{r['seqname']}\t{r['start']}\t"
                    f"{r['end']}\t{r['strand']}\t{r['subtype']}\t"
                    f"{r['length']}\t"
                    f"{'NA' if r['score']  is None else r['score']}\t"
                    f"{'NA' if r['evalue'] is None else r['evalue']}\n"
                )
    scored = [r for r in rrna_records if r["score"] is not None]
    _log(f"  {len(scored):,} / {len(rrna_records):,} rRNA copies have bit scores")
    for sub in _RRNA_ORDER:
        sub_scores = [r["score"] for r in scored if r["subtype"] == sub]
        if sub_scores:
            med = sorted(sub_scores)[len(sub_scores) // 2]
            _log(f"    {sub:30s}  n={len(sub_scores):6,}  "
                 f"median score={med:.1f}  "
                 f"min={min(sub_scores):.1f}  max={max(sub_scores):.1f}")


def _run_trna_classify_analysis(trna_records: list, fasta: Path,
                                 rfam_dir: Path, workdir: Path, results: Path,
                                 prefix: str, threads: int, force: bool) -> list:
    """Score, classify tRNA records; write TSV; return classified list."""
    if not trna_records:
        _log("  No tRNA records found — skipping 6b")
        return []
    _log(f"  Extracting {len(trna_records):,} tRNA sequences …")
    trna_seqs = _extract_sequences(fasta, trna_records)
    _log(f"  Extracted {len(trna_seqs):,} sequences")

    trna_cm   = _ensure_trna_cm(rfam_dir)
    cm_scores = _score_trnas_cmsearch(trna_seqs, trna_cm, workdir, threads, force)
    _log(f"  cmsearch scored {len(cm_scores):,} tRNAs")

    trna_classified = _classify_trnas(trna_records, cm_scores)
    n_pseudo = sum(1 for t in trna_classified if t["category"] == "pseudogene")
    n_func   = len(trna_classified) - n_pseudo
    _log(f"  Functional: {n_func:,}   Pseudogene candidates: {n_pseudo:,}")

    reason_counts = Counter()
    for t in trna_classified:
        if t["category"] == "pseudogene":
            for tok in t["reason"].split(";"):
                reason_counts[tok.strip()] += 1
    for reason, cnt in reason_counts.most_common():
        _log(f"    {reason}: {cnt:,}")

    trna_class_tsv = results / f"mod06_trna_class_{prefix}.tsv"
    if force or not trna_class_tsv.exists():
        sorted_trna = sorted(trna_classified,
                             key=lambda t: (_natural_key(t["seqname"]), t["start"]))
        with open(trna_class_tsv, "w") as fh:
            fh.write("feature_id\tseqname\tstart\tend\tstrand\tsubtype\t"
                     "anticodon\tlength_bp\tcm_score\tcategory\treason\n")
            for t in sorted_trna:
                sc = f"{t['cm_score']:.1f}" if t["cm_score"] is not None else "NA"
                fh.write(
                    f"{t['feat_id']}\t{t['seqname']}\t{t['start']}\t"
                    f"{t['end']}\t{t['strand']}\t{t['subtype']}\t"
                    f"{t['anticodon']}\t{t['length']}\t"
                    f"{sc}\t{t['category']}\t{t['reason']}\n"
                )
        _log(f"  Written: {trna_class_tsv.name}")
    return trna_classified


def _run_array_analysis(rrna_records: list, trna_records: list,
                        results: Path, prefix: str,
                        force: bool) -> tuple[list, list]:
    """Detect tandem arrays, log stats, write arrays TSV; return (rrna_arrays, trna_arrays)."""
    rrna_arrays = _detect_arrays(rrna_records, max_gap=_RRNA_ARRAY_GAP) if rrna_records else []
    trna_arrays = _detect_arrays(trna_records, max_gap=_TRNA_ARRAY_GAP) if trna_records else []

    if rrna_arrays:
        total_units = sum(a["complete_rDNA_units"] for a in rrna_arrays)
        _log(f"  rRNA: {len(rrna_arrays)} arrays on "
             f"{len({a['seqname'] for a in rrna_arrays})} sequences  |  "
             f"~{total_units} complete rDNA units")
        sizes = [a["n_copies"] for a in rrna_arrays]
        _log(f"    copies/array: min={min(sizes)}  median={sorted(sizes)[len(sizes)//2]}  "
             f"max={max(sizes)}")
        all_sp = [s for a in rrna_arrays for s in a["spacings"]]
        if all_sp:
            med_sp = sorted(all_sp)[len(all_sp) // 2]
            _log(f"    inter-copy spacing: median={med_sp/1000:.1f} kb  "
                 f"min={min(all_sp)/1000:.1f} kb  max={max(all_sp)/1000:.1f} kb")
    else:
        _log("  No rRNA arrays detected")

    if trna_arrays:
        _log(f"  tRNA: {len(trna_arrays)} arrays on "
             f"{len({a['seqname'] for a in trna_arrays})} sequences")
    else:
        _log("  No tRNA arrays detected")

    arrays_tsv = results / f"mod06_arrays_{prefix}.tsv"
    all_arrays = ([{"feature_class": "rRNA", **a} for a in rrna_arrays] +
                  [{"feature_class": "tRNA", **a} for a in trna_arrays])
    all_arrays.sort(key=lambda a: (_natural_key(a["seqname"]), a["array_start"]))
    if all_arrays and (force or not arrays_tsv.exists()):
        with open(arrays_tsv, "w") as fh:
            fh.write("feature_class\tarray_id\tseqname\tarray_start\tarray_end\t"
                     "n_copies\tarray_length_bp\tmean_spacing_bp\t"
                     "min_spacing_bp\tmax_spacing_bp\t"
                     "complete_rDNA_units\tsubtype_counts\n")
            for a in all_arrays:
                sub = ";".join(f"{k}:{v}" for k, v in sorted(a["subtypes"].items()))
                fh.write(
                    f"{a['feature_class']}\t{a['array_id']}\t{a['seqname']}\t"
                    f"{a['array_start']}\t{a['array_end']}\t{a['n_copies']}\t"
                    f"{a['array_length_bp']}\t{a['mean_spacing_bp']}\t"
                    f"{a['min_spacing_bp']}\t{a['max_spacing_bp']}\t"
                    f"{a.get('complete_rDNA_units', 0)}\t{sub}\n"
                )
        _log(f"  Written: {arrays_tsv.name}")
    return rrna_arrays, trna_arrays


# ── Module 6 main ─────────────────────────────────────────────────────────────

def run_module6_evolution(fasta: Path,
                          rrna_gff: Path, trna_gff: Path,
                          rfam_dir: Path, workdir: Path, results: Path,
                          prefix: str, genome_size: int, threads: int,
                          plot_formats: list, force: bool) -> None:
    """Run evolutionary analysis: rRNA scores, tRNA pseudogene classification, tandem arrays."""
    _banner("Module 6a — rRNA Score Distribution")
    rrna_records = _parse_gff3_records(rrna_gff)
    _run_rrna_score_analysis(rrna_records, results, prefix, force)

    _banner("Module 6b — tRNA Pseudogene Classification")
    trna_records    = _parse_gff3_records(trna_gff)
    trna_classified = _run_trna_classify_analysis(
        trna_records, fasta, rfam_dir, workdir, results, prefix, threads, force)

    _banner("Module 6c — Tandem Array Detection")
    rrna_arrays, trna_arrays = _run_array_analysis(
        rrna_records, trna_records, results, prefix, force)

    _banner("Module 6d — Evolution Plots")
    out_base  = results / f"mod06_evolution_{prefix}"
    out_paths = [out_base.with_suffix(f".{fmt}") for fmt in plot_formats]
    if not _checkpoint(out_paths[0], "evolution plots", force):
        _plot_evolution(rrna_records, trna_classified,
                        rrna_arrays, trna_arrays,
                        out_base, plot_formats, prefix)


# ── Module 7: Subtelomeric tandem repeats ─────────────────────────────────────
#
# Corroborating evidence for telomere calls, and a coarse genome-completeness
# signal at chromosome ends: subtelomeric satellite/tandem-repeat arrays are
# commonly (not universally) found adjacent to true telomeres in eukaryotic
# genomes, and can still be present even where the assembly stops just short
# of a fully resolved canonical telomere array, or where that array is too
# degraded for Module 0's strict k-mer scan to call confidently. Every
# scaffold end is classified into one of three completeness tiers:
#   Tier 1  Module 0 found a confirmed telomere.
#   Tier 2  No confirmed telomere, but a genuine tandem repeat (TRF) was
#           found in the terminal window -- some evidence of proximity to a
#           chromosome end, though weaker than a direct telomere call.
#   Tier 3  Neither -- no evidence either way.

_TRF_PARAMS = ("2", "7", "7", "80", "10", "50", "2000")  # match mismatch delta PM PI minscore maxperiod
_TRF_DAT_RE = re.compile(
    r"^(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
)


def _write_terminal_windows_fasta(fasta: Path, window_bp: int,
                                  out_fasta: Path) -> dict:
    """Stream the genome once and write a small FASTA of just the terminal
    (5'/3') windows of every sequence, with synthetic headers
    '{seqname}::5prime' / '{seqname}::3prime' so TRF hits can be mapped back
    to (seqname, end, offset). Returns {synthetic_header: (seqname, end_label,
    offset_bp)} for coordinate translation. Bounded memory: only the (small)
    terminal windows are ever held/written, never the full genome."""
    offsets: dict = {}
    with open(out_fasta, "w") as out:
        for name, seq in _iter_fasta(fasta):
            length = len(seq)
            for end_label, region_seq, offset in [
                ("5prime", seq[:window_bp], 0),
                ("3prime", seq[max(0, length - window_bp):], max(0, length - window_bp)),
            ]:
                if not region_seq:
                    continue
                header = f"{name}::{end_label}"
                offsets[header] = (name, end_label, offset)
                out.write(f">{header}\n{region_seq}\n")
    return offsets


def _run_trf(windows_fasta: Path, workdir: Path) -> Path | None:
    """Run TRF on the terminal-windows FASTA; return the .dat output path, or
    None if TRF did not produce one.

    TRF is a known quirk among wrapped tools: it returns a non-zero exit
    code on ordinary *success* (the code reflects the number of repeats
    found), so it is invoked directly rather than through the shared _run()
    helper, which would otherwise treat that as a hard failure. TRF also
    names its own output file from the input filename + parameter string
    (e.g. 'input.fa.2.7.7.80.10.50.2000.dat'), which is awkward to predict
    exactly across TRF versions -- so the real output file is found by
    globbing workdir for '*.dat' after the run, with cwd=workdir controlling
    where TRF's side-effect files land (per the shared tool-wrapping
    convention used elsewhere in this pipeline)."""
    trf = _require_tool("trf")
    cmd = [trf, str(windows_fasta.resolve()), *_TRF_PARAMS, "-d", "-h"]
    _log(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                   text=True, cwd=workdir)
    dat_files = sorted(workdir.glob("*.dat"))
    if not dat_files:
        _log("  WARNING: TRF did not produce a .dat output file — "
             "skipping subtelomeric tandem repeat detection")
        return None
    return dat_files[0]


def _parse_trf_dat(dat_path: Path) -> dict:
    """Parse a TRF .dat file into {synthetic_header: [hit dicts]}."""
    hits: dict = {}
    current = None
    with open(dat_path) as fh:
        for line in fh:
            if line.startswith("Sequence:"):
                current = line.split(":", 1)[1].strip().split()[0]
                hits.setdefault(current, [])
                continue
            if current is None:
                continue
            m = _TRF_DAT_RE.match(line)
            if not m:
                continue
            start, end, period, copies_str = m.group(1), m.group(2), m.group(3), m.group(4)
            pct_matches = m.group(6)
            hits[current].append({
                "start":       int(start),
                "end":         int(end),
                "period":      int(period),
                "copies":      float(copies_str),
                "pct_matches": int(pct_matches),
            })
    return hits


def run_module7_subtelomeric(fasta: Path, tel_gff: Path | None,
                             window_bp: int, min_copies: float,
                             results: Path, workdir: Path,
                             prefix: str, force: bool) -> tuple[Path, Path, dict]:
    """Scan terminal windows for subtelomeric tandem repeats (TRF), and
    classify every scaffold end into a telomere-completeness tier.
    Returns (gff3_path, summary_tsv_path, counts)."""
    out_gff3 = results / f"mod07_subtelomeric_{prefix}.gff3"
    out_tsv  = results / f"mod07_completeness_{prefix}.tsv"
    if _checkpoint(out_gff3, "subtelomeric-scan", force) and out_tsv.exists():
        return out_gff3, out_tsv, {}

    seq_lengths = _fasta_seq_lengths(fasta)
    _log(f"  {len(seq_lengths)} sequences")

    # Which (seqname, end) pairs already have a confirmed Module 0 telomere?
    confirmed_ends: set = set()
    if tel_gff is not None and tel_gff.exists():
        with open(tel_gff) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                name_attr = _gff3_name(line)  # e.g. "telomere_5prime"
                if name_attr.startswith("telomere_"):
                    end_label = name_attr[len("telomere_"):]
                    confirmed_ends.add((cols[0], end_label))
    _log(f"  Confirmed telomeres (Module 0): {len(confirmed_ends)} scaffold end(s)")

    windows_fasta = workdir / "subtelomeric_windows.fasta"
    offsets = _write_terminal_windows_fasta(fasta, window_bp, windows_fasta)

    dat_path = _run_trf(windows_fasta, workdir)
    trf_hits_by_header = _parse_trf_dat(dat_path) if dat_path is not None else {}

    records: list = []
    summary_rows: list = []
    tel_count = 0
    tier_counts = {1: 0, 2: 0, 3: 0}

    for name, length in seq_lengths.items():
        for end_label in ("5prime", "3prime"):
            header = f"{name}::{end_label}"
            if header not in offsets:
                continue  # scaffold too short to have this end scanned
            _, _, offset = offsets[header]
            has_telomere = (name, end_label) in confirmed_ends

            hits = [h for h in trf_hits_by_header.get(header, [])
                   if h["copies"] >= min_copies]
            best_hit = max(hits, key=lambda h: h["copies"], default=None)

            if has_telomere:
                tier = 1
            elif best_hit is not None:
                tier = 2
            else:
                tier = 3
            tier_counts[tier] += 1

            for h in hits:
                tel_count += 1
                abs_s = offset + h["start"]
                abs_e = offset + h["end"]
                attrs = {
                    "ID":           f"subtel_{name}_{end_label}_{tel_count}",
                    "Name":         f"subtelomeric_repeat_{end_label}",
                    "period_size":  h["period"],
                    "copy_number":  f"{h['copies']:.1f}",
                    "pct_matches":  h["pct_matches"],
                }
                records.append(
                    _gff3_record(name, "UbboTELORNA", "subtelomeric_tandem_repeat",
                                 abs_s, abs_e, h["pct_matches"], ".", ".", attrs)
                )

            summary_rows.append({
                "seqname":            name,
                "end":                end_label,
                "seq_length_bp":      length,
                "tier":               tier,
                "confirmed_telomere": "yes" if has_telomere else "no",
                "best_period_bp":     best_hit["period"] if best_hit else "",
                "best_copy_number":   f"{best_hit['copies']:.1f}" if best_hit else "",
            })

    records.sort(key=_gff3_sort_key)
    with open(out_gff3, "w") as fh:
        fh.write(_GFF3_HEADER)
        for r in records:
            fh.write(r + "\n")

    with open(out_tsv, "w") as fh:
        fh.write("seqname\tend\tseq_length_bp\ttier\tconfirmed_telomere\t"
                 "best_period_bp\tbest_copy_number\n")
        for row in summary_rows:
            fh.write("\t".join(str(row[k]) for k in
                     ("seqname", "end", "seq_length_bp", "tier",
                      "confirmed_telomere", "best_period_bp", "best_copy_number")) + "\n")

    n_ends = sum(tier_counts.values())
    pct = {t: (100.0 * c / n_ends if n_ends else 0.0) for t, c in tier_counts.items()}
    _log(f"  Subtelomeric tandem repeats found: {tel_count}")
    _log(f"  Scaffold-end completeness — "
        f"Tier 1 (confirmed telomere): {tier_counts[1]} ({pct[1]:.1f}%), "
        f"Tier 2 (subtelomeric support only): {tier_counts[2]} ({pct[2]:.1f}%), "
        f"Tier 3 (neither): {tier_counts[3]} ({pct[3]:.1f}%)")

    counts = {
        "n_subtelomeric":  tel_count,
        "tier_counts":     {str(k): v for k, v in tier_counts.items()},
        "tier_pct":        {str(k): round(v, 1) for k, v in pct.items()},
    }
    return out_gff3, out_tsv, counts


# ── Run log and carbon tracker helpers ───────────────────────────────────────

def _open_run_log(logs_dir: Path) -> Path:
    """Open the run log, write the header, set _LOG_FH; return log path."""
    global _LOG_FH
    log_path = logs_dir / "Run_UbboTELORNA.log"
    _LOG_FH  = open(log_path, "w")
    sep = "=" * 62
    _LOG_FH.write(f"{sep}\n  UbboTELORNA {VERSION}  —  Run Log\n{sep}\n")
    _LOG_FH.write(f"Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    _LOG_FH.write(f"User      : {getpass.getuser()}\n")
    _LOG_FH.write(f"Server    : {platform.node()}\n")
    _LOG_FH.write(f"OS        : {platform.system()} {platform.release()} "
                  f"({platform.machine()})\n")
    _LOG_FH.write(f"Directory : {os.getcwd()}\n")
    _LOG_FH.write(f"Command   : {' '.join(sys.argv)}\n")
    _LOG_FH.write(f"{sep}\n\n")
    _LOG_FH.flush()
    return log_path


def _setup_co2_tracker(logs_dir: Path, prefix: str, disable: bool):
    """Start the codecarbon emissions tracker if available; return tracker or None."""
    if disable:
        _log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
        return None
    try:
        # pkg_resources is part of setuptools; in some conda environments
        # it is not on sys.path even though setuptools is installed.
        # Inject a minimal shim so codecarbon can import cleanly.
        try:
            import pkg_resources  # noqa: F401
        except ModuleNotFoundError:
            import types as _t, importlib.metadata as _m, importlib as _il
            from pathlib import Path as _P
            _shim = _t.ModuleType("pkg_resources")
            def _get_dist(name):
                try:
                    d = _m.distribution(name)
                    d.version = d.metadata["Version"]
                    return d
                except Exception:
                    return None
            def _resource_filename(pkg, resource):
                try:
                    mod = _il.import_module(pkg)
                    return str(_P(mod.__file__).parent / resource)
                except Exception:
                    return resource
            _shim.get_distribution     = _get_dist
            _shim.resource_filename    = _resource_filename
            _shim.DistributionNotFound = Exception
            sys.modules["pkg_resources"] = _shim

        from codecarbon import EmissionsTracker
        tracker = EmissionsTracker(
            output_dir=str(logs_dir),
            output_file=f"{prefix}.emissions.csv",
            project_name="UbboTELORNA",
            log_level="warning",
        )
        tracker.start()
        _log("  codecarbon tracker started")
        return tracker
    except ImportError as e:
        _log(f"  codecarbon not installed — carbon tracking skipped ({e})")
    except Exception as e:
        _log(f"  codecarbon failed to start — carbon tracking skipped ({e})")
    return None


# ── Argument parser ───────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="UbboTELORNA",
        description=(
            "Annotate telomeres, rRNA, and tRNA in a genome assembly.\n\n"
            "Modules:\n"
            "  0  Telomere identification  (k-mer scan at contig ends)\n"
            "  1  Low-complexity masking   (tantan)\n"
            "  2  rRNA annotation          (nhmmer [default] or cmsearch + Rfam profiles)\n"
            "  3  tRNA annotation          (ARAGORN)\n"
            "  4  Integration              (merged GFF3 + summary table)\n"
            "  5  Visualization            (ideogram, subtype bars, composition donut)\n"
            "  6  Evolutionary analysis    (rRNA scores, tRNA pseudogenes, tandem arrays)\n"
            "  7  Subtelomeric tandem repeats (TRF; telomere-completeness tiering)\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    req = ap.add_argument_group("Required")
    req.add_argument("--fasta",   required=True, type=Path,
                     help="Input genome assembly FASTA")
    req.add_argument("--output",  required=True, type=Path,
                     help="Output directory (created if absent)")

    tel = ap.add_argument_group("Module 0 — Telomere")
    tel.add_argument("--telomere_repeat", default=None,
                     help="Known repeat unit (default: auto-detect). "
                          "Example: TTTAGGG (plants), TTAGGG (vertebrates)")
    tel.add_argument("--telomere_window", type=int, default=10_000,
                     help="bp to scan at each contig end (default: 10000)")
    tel.add_argument("--telomere_density", type=float, default=0.5,
                     help="Minimum repeat density to call a telomere (default: 0.5)")
    tel.add_argument("--telomere_min_len", type=int, default=100,
                     help="Minimum telomere length to report (default: 100 bp)")

    subtel = ap.add_argument_group("Module 7 — Subtelomeric tandem repeats")
    subtel.add_argument("--subtelomeric_window_bp", type=int, default=20_000,
                        help="bp to scan at each scaffold end for subtelomeric "
                             "tandem repeats via TRF (default: 20000). Larger "
                             "than --telomere_window since satellite arrays can "
                             "sit further from the true terminus than the "
                             "telomere repeat itself.")
    subtel.add_argument("--subtelomeric_min_copies", type=float, default=3.0,
                        help="Minimum tandem copy number (as reported by TRF) "
                             "for a repeat to be reported and counted toward "
                             "Tier 2 completeness (default: 3.0)")

    rrna = ap.add_argument_group("Module 2 — rRNA")
    rrna.add_argument("--kingdom", choices=["euka", "bacteria", "archaea"],
                      default="euka",
                      help="Organism kingdom — selects Rfam model set (default: euka)")
    rrna.add_argument("--evalue", type=float, default=1e-5,
                      help="E-value threshold for cmsearch (default: 1e-5)")
    rrna.add_argument("--rfam_dir", type=Path, default=None,
                      help=f"Directory with pre-downloaded Rfam .cm files "
                           f"(default: auto-download to {_DEFAULT_CACHE})")
    rrna.add_argument("--search_tool", choices=["nhmmer", "cmsearch"],
                      default="nhmmer",
                      help="rRNA search tool: nhmmer (default, faster, lower "
                           "memory, uses Rfam HMM profiles) or cmsearch "
                           "(higher sensitivity, uses Rfam covariance models)")
    rrna.add_argument("--cmsearch_mxsize", type=float, default=512.0,
                      help="Max DP matrix size per cmsearch thread in Mb — "
                           "only used with --search_tool cmsearch "
                           "(default: 512)")
    rrna.add_argument("--flag_5s_arrays", dest="flag_5s_arrays",
                      action="store_true", default=True,
                      help="Flag 5S rRNA tandem array members in GFF3 attributes "
                           "(default: on; use --no_flag_5s_arrays to disable)")
    rrna.add_argument("--no_flag_5s_arrays", dest="flag_5s_arrays",
                      action="store_false",
                      help="Disable 5S tandem array flagging")
    rrna.add_argument("--5s_array_min_copies", dest="s5_array_min_copies",
                      type=int, default=5,
                      help="Minimum consecutive copies to define a 5S array "
                           "(default: 5)")
    rrna.add_argument("--5s_array_max_gap", dest="s5_array_max_gap",
                      type=int, default=5000,
                      help="Maximum inter-copy gap to define a 5S array, bp "
                           "(default: 5000)")
    rrna.add_argument("--cap_5s", type=int, default=0,
                      help="If >0, cap 5S predictions per sequence to this number, "
                           "keeping the highest-scoring hits (default: 0 = no cap)")

    gen = ap.add_argument_group("General")
    gen.add_argument("--threads", type=int, default=4,
                     help="CPU threads for cmsearch (default: 4)")
    gen.add_argument("--skip_module", default="",
                     help="Comma-separated list of module numbers to skip "
                          "(e.g. --skip_module 0,1,2). "
                          "0=telomere  1=masking  2=rRNA  3=tRNA  "
                          "4=integration  5=visualization  6=evolution  "
                          "7=subtelomeric tandem repeats. "
                          "Skipped modules are not rerun, but their outputs "
                          "from previous runs are picked up automatically.")
    gen.add_argument("--format", default="pdf",
                     help="Plot format(s): pdf, png, svg — comma-separated "
                          "(default: pdf)")
    gen.add_argument("--top_sequences", type=int, default=50,
                     help="Number of sequences to show in the ideogram (default: 50)")
    gen.add_argument("--sort_sequences", choices=["length", "seqid"], default="length",
                     help="Order sequences in the ideogram: 'length' (longest first) "
                          "or 'seqid' (natural alphabetic/numeric sort); default: length")
    gen.add_argument("--chunk_size", type=int, default=200,
                     help="Sequences per chunk fed to nhmmer/cmsearch "
                          "(default: 200; set to 0 to disable chunking and "
                          "pass the full genome to the tool in one run)")
    gen.add_argument("--chunk_max_bp", type=int, default=50_000_000,
                     help="Also flush a chunk once its cumulative sequence "
                          "length reaches this many bp (default: 50,000,000), "
                          "regardless of --chunk_size. Chromosome-scale "
                          "assemblies with few, very large sequences never "
                          "reach the sequence-count threshold, so the whole "
                          "genome would otherwise land in a single chunk; "
                          "this caps peak nhmmer/cmsearch memory on those "
                          "genomes. Set to 0 to disable the bp cap.")
    gen.add_argument("--force", action="store_true",
                     help="Rerun all steps even if outputs already exist")
    gen.add_argument("--dry_run", action="store_true",
                     help="Validate inputs and print steps, then exit")
    gen.add_argument("--disable_co2_tracking", action="store_true",
                     help="Disable carbon footprint tracking (codecarbon)")
    gen.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    return ap


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning, message=".*pynvml.*")
    _print_quote()

    ap   = _build_parser()
    args = ap.parse_args()

    # Resolve paths
    args.fasta  = args.fasta.resolve()
    run_dir     = args.output.resolve()

    _validate_inputs([("--fasta", args.fasta)])

    # Create output layout
    results  = run_dir / "results"
    workdir  = run_dir / "workdir"
    logs_dir = run_dir / "logs"
    for d in (results, workdir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    prefix = run_dir.name

    log_path = _open_run_log(logs_dir)

    _banner(f"UbboTELORNA {VERSION}")
    _log(f"  Input FASTA : {args.fasta}")
    _log(f"  Output dir  : {run_dir}")
    _log(f"  Kingdom     : {args.kingdom}")
    _log(f"  rRNA tool   : {args.search_tool}")
    chunk_desc = args.chunk_size if args.chunk_size > 0 else "disabled"
    if args.chunk_max_bp > 0:
        chunk_desc = f"{chunk_desc} (max {args.chunk_max_bp:,} bp/chunk)"
    _log(f"  Chunk size  : {chunk_desc}")
    _log(f"  Threads     : {args.threads}")
    genome_size = _fasta_total_length(args.fasta)
    _log(f"  Genome size : {genome_size:,} bp")
    rfam_dir = args.rfam_dir or _DEFAULT_CACHE
    rfam_dir.mkdir(parents=True, exist_ok=True)

    # Parse --skip_module into a set of ints
    skip_modules: set = set()
    if args.skip_module.strip():
        for tok in args.skip_module.split(","):
            tok = tok.strip()
            if tok:
                try:
                    skip_modules.add(int(tok))
                except ValueError:
                    print(f"ERROR: --skip_module value '{tok}' is not a module number",
                          file=sys.stderr)
                    sys.exit(1)
    if skip_modules:
        _log(f"  Skipping modules: {sorted(skip_modules)}")

    if args.force:
        _log("  --force set: all steps will rerun regardless of existing outputs")
    elif workdir.exists() and any(workdir.iterdir()):
        _log("  Existing workdir found — resuming from checkpoints "
             "(use --force to rerun all steps from scratch)")

    # Dry run
    if args.dry_run:
        _banner("Dry run — no steps will be executed")
        _log(f"  FASTA       : {args.fasta}")
        _log(f"  Output      : {run_dir}/")
        _log("  Steps that would run:")
        if 0 not in skip_modules:
            _log("    [0] Telomere identification  →  results/mod00_telomeres_*.gff3")
        if 1 not in skip_modules:
            _log("    [1] Low-complexity masking   →  workdir/masked_soft.fasta")
        if 2 not in skip_modules:
            chunk_desc = (f"chunk_size={args.chunk_size}"
                          if args.chunk_size > 0 else "no chunking")
            _log(f"    [2] rRNA annotation ({args.search_tool}, {chunk_desc})"
                 f"  →  results/mod02_rRNA_*.gff3")
        if 3 not in skip_modules:
            _log("    [3] tRNA annotation          →  results/mod03_tRNA_*.gff3")
        if 4 not in skip_modules:
            _log("    [4] Integration              →  results/mod04_annotation_*.gff3")
        if 5 not in skip_modules:
            fmt0 = args.format.split(",")[0].strip()
            _log(f"    [5] Visualization ({args.sort_sequences})  →  results/mod05_plot_*.{fmt0}")
        if 6 not in skip_modules:
            _log("    [6] Evolutionary analysis       →  results/mod06_*.tsv + mod06_evolution_*")
        if 7 not in skip_modules:
            _log("    [7] Subtelomeric tandem repeats (TRF)  →  "
                 "results/mod07_subtelomeric_*.gff3 + mod07_completeness_*.tsv")
        _log("  Exiting (--dry_run).")
        if _LOG_FH:
            _LOG_FH.close()
        sys.exit(0)

    # Carbon tracking
    t_start  = time.monotonic()
    _tracker = _setup_co2_tracker(logs_dir, prefix, args.disable_co2_tracking)

    def _skip(n: int, label: str, path: Path | None = None) -> bool:
        if n in skip_modules:
            if path and path.exists() and path.stat().st_size > 0:
                _log(f"  [Module {n}] skipped — picking up existing {path.name}")
            else:
                _log(f"  [Module {n}] skipped ({label})")
            return True
        return False

    # ── Module 0: Telomere identification ─────────────────────────────────────
    tel_gff     = results / f"mod00_telomeres_{prefix}.gff3"
    repeat_used = None
    if not _skip(0, "telomere identification", tel_gff):
        _banner("Module 0 — Telomere Identification")
        tel_gff, repeat_used = run_module0_telomeres(
            fasta       = args.fasta,
            repeat_unit = args.telomere_repeat,
            tel_window  = args.telomere_window,
            tel_density = args.telomere_density,
            tel_min_len = args.telomere_min_len,
            results     = results,
            workdir     = workdir,
            prefix      = prefix,
            force       = args.force,
        )

    # ── Module 7: Subtelomeric tandem repeats ─────────────────────────────────
    # Runs right after Module 0, the only module it depends on -- it uses
    # tel_gff to know which scaffold ends already have a confirmed telomere.
    # Kept separate from Module 4 (Integration)'s merged tel/rRNA/tRNA GFF3;
    # not folded in there since that module's existing contract wasn't part
    # of this addition's scope.
    subtel_gff = results / f"mod07_subtelomeric_{prefix}.gff3"
    subtel_counts: dict = {}
    if not _skip(7, "subtelomeric tandem repeats", subtel_gff):
        _banner("Module 7 — Subtelomeric Tandem Repeats")
        subtel_gff, _subtel_tsv, subtel_counts = run_module7_subtelomeric(
            fasta       = args.fasta,
            tel_gff     = tel_gff,
            window_bp   = args.subtelomeric_window_bp,
            min_copies  = args.subtelomeric_min_copies,
            results     = results,
            workdir     = workdir,
            prefix      = prefix,
            force       = args.force,
        )

    # ── Module 1: Low-complexity masking ──────────────────────────────────────
    masked_soft = workdir / "masked_soft.fasta"
    if not _skip(1, "low-complexity masking", masked_soft):
        _banner("Module 1 — Low-Complexity Masking")
        run_module1_masking(
            fasta   = args.fasta,
            workdir = workdir,
            force   = args.force,
        )

    # ── Module 2: rRNA annotation ─────────────────────────────────────────────
    rrna_gff = results / f"mod02_rRNA_{prefix}.gff3"
    if not _skip(2, "rRNA annotation", rrna_gff):
        _banner("Module 2 — rRNA Annotation")
        rrna_gff = run_module2_rrna(
            fasta            = args.fasta,
            kingdom          = args.kingdom,
            threads          = args.threads,
            evalue           = args.evalue,
            rfam_dir         = args.rfam_dir,
            mxsize           = args.cmsearch_mxsize,
            search_tool      = args.search_tool,
            workdir          = workdir,
            results          = results,
            prefix           = prefix,
            force            = args.force,
            chunk_size       = args.chunk_size,
            chunk_max_bp     = args.chunk_max_bp,
            flag_5s_arrays   = args.flag_5s_arrays,
            min_5s_copies    = args.s5_array_min_copies,
            max_5s_gap       = args.s5_array_max_gap,
            cap_5s           = args.cap_5s,
        )

    # ── Module 3: tRNA annotation ─────────────────────────────────────────────
    trna_gff = results / f"mod03_tRNA_{prefix}.gff3"
    if not _skip(3, "tRNA annotation", trna_gff):
        _banner("Module 3 — tRNA Annotation")
        trna_gff = run_module3_trna(
            fasta   = args.fasta,
            workdir = workdir,
            results = results,
            prefix  = prefix,
            force   = args.force,
        )

    # ── Module 4: Integration ─────────────────────────────────────────────────
    counts: dict = {}
    summary_tsv = results / f"mod04_summary_{prefix}.tsv"
    if not _skip(4, "integration", summary_tsv):
        _banner("Module 4 — Integration")
        _, counts = run_module4_integration(
            tel_gff, rrna_gff, trna_gff, results, prefix,
            genome_size=genome_size)

    # ── Module 5: Visualization ───────────────────────────────────────────────
    plot_formats = [f.strip().lstrip(".") for f in args.format.split(",")]
    out_plot = results / f"mod05_plot_{prefix}.{plot_formats[0]}"
    if not _skip(5, "visualization", out_plot):
        _banner("Module 5 — Visualization")
        run_module5_plot(
            fasta         = args.fasta,
            tel_gff       = tel_gff,
            rrna_gff      = rrna_gff,
            trna_gff      = trna_gff,
            summary_tsv   = summary_tsv,
            genome_size   = genome_size,
            results       = results,
            prefix        = prefix,
            plot_formats  = plot_formats,
            top_sequences = args.top_sequences,
            sort_by       = args.sort_sequences,
            force         = args.force,
        )

    # ── Module 6: Evolutionary analysis ───────────────────────────────────────
    evo_out = results / f"mod06_rrna_scores_{prefix}.tsv"
    if not _skip(6, "evolutionary analysis", evo_out):
        run_module6_evolution(
            fasta        = args.fasta,
            rrna_gff     = rrna_gff,
            trna_gff     = trna_gff,
            rfam_dir     = rfam_dir,
            workdir      = workdir,
            results      = results,
            prefix       = prefix,
            genome_size  = genome_size,
            threads      = args.threads,
            plot_formats = plot_formats,
            force        = args.force,
        )

    # ── Resource usage & run summary ──────────────────────────────────────────
    elapsed_s   = time.monotonic() - t_start
    ru          = resource.getrusage(resource.RUSAGE_SELF)
    peak_mem_mb = (ru.ru_maxrss / (1024 * 1024)
                   if platform.system() == "Darwin"
                   else ru.ru_maxrss / 1024)

    emissions_kg = None
    if _tracker is not None:
        try:
            emissions_kg = _tracker.stop()
        except Exception:
            pass

    summary = {
        "date":    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version": VERSION,
        "input_fasta":    str(args.fasta),
        "genome_size_bp": genome_size,
        "kingdom":        args.kingdom,
        "telomere_repeat": repeat_used,
        "parameters": {
            "telomere_window":  args.telomere_window,
            "telomere_density": args.telomere_density,
            "telomere_min_len": args.telomere_min_len,
            "evalue":           args.evalue,
            "search_tool":      args.search_tool,
            "threads":          args.threads,
        },
        "feature_counts": {
            "telomere": {
                "total":          counts.get("n_tel",  0),
                "total_length_bp": counts.get("tel_len", 0),
                "by_end":         counts.get("tel_by_end",     {}),
                "len_by_end":     counts.get("tel_len_by_end", {}),
            },
            "rRNA": {
                "total":           counts.get("n_rrna", 0),
                "total_length_bp": counts.get("rrna_len", 0),
                "by_type":         counts.get("rrna_by_type",     {}),
                "len_by_type":     counts.get("rrna_len_by_type", {}),
            },
            "tRNA": {
                "total":           counts.get("n_trna", 0),
                "total_length_bp": counts.get("trna_len", 0),
                "by_type":         counts.get("trna_by_type",     {}),
                "len_by_type":     counts.get("trna_len_by_type", {}),
            },
        } if counts else {},
        "telomere_completeness": {
            "total_subtelomeric_repeats": subtel_counts.get("n_subtelomeric", 0),
            "scaffold_ends": {
                "tier1_confirmed_telomere":       subtel_counts.get("tier_counts", {}).get("1", 0),
                "tier2_subtelomeric_support_only": subtel_counts.get("tier_counts", {}).get("2", 0),
                "tier3_neither":                   subtel_counts.get("tier_counts", {}).get("3", 0),
            },
            "scaffold_ends_pct": {
                "tier1_confirmed_telomere":       subtel_counts.get("tier_pct", {}).get("1", 0.0),
                "tier2_subtelomeric_support_only": subtel_counts.get("tier_pct", {}).get("2", 0.0),
                "tier3_neither":                   subtel_counts.get("tier_pct", {}).get("3", 0.0),
            },
        } if subtel_counts else {},
        "resource_usage": {
            "wall_clock_s":       round(elapsed_s, 1),
            "peak_mem_mb":        round(peak_mem_mb, 1),
            "emissions_kg_CO2eq": emissions_kg,
        },
    }
    summary_path = results / f"{prefix}.run_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    _banner("Run complete")
    _log(f"  Wall clock  : {elapsed_s:.1f} s")
    _log(f"  Peak memory : {peak_mem_mb:.1f} MB")
    _log(f"  Results     : {results}/")
    _log(f"  Log         : {log_path}")

    if _LOG_FH is not None:
        _LOG_FH.close()


if __name__ == "__main__":
    main()
