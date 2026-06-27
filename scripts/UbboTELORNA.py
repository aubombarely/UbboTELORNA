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

VERSION = "v0.1.0"

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
    """Extract the Name= attribute value from a GFF3 data line.
    Strips any trailing parenthesised anticodon, e.g. tRNA-Phe(GAA) → tRNA-Phe.
    Returns empty string if Name is absent.
    """
    attrs = line.split("\t")[8] if line.count("\t") >= 8 else ""
    for field in attrs.split(";"):
        if field.startswith("Name="):
            name = field[5:].strip()
            return re.sub(r"\([^)]+\)$", "", name)
    return ""


def _gff3_sort_key(line: str) -> tuple:
    """Natural-sort key for a GFF3 data line: (seqid parts, start).

    Splits the SeqID on digit/non-digit boundaries so that e.g.
    Chr2 sorts before Chr10 rather than after it.
    """
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
                 start: int, end: int, score, strand: str,
                 frame: str, attrs: dict) -> str:
    score_str = f"{score:.2f}" if isinstance(score, (int, float)) else "."
    attr_str  = ";".join(f"{k}={v}" for k, v in attrs.items())
    return f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score_str}\t{strand}\t{frame}\t{attr_str}"


# ── Telomere detection helpers ────────────────────────────────────────────────

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _canonical(kmer: str) -> str:
    """Lexicographically smaller of kmer and its reverse complement."""
    rc = _revcomp(kmer.upper())
    return min(kmer.upper(), rc)


def _all_rotations(kmer: str) -> set:
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


def _detect_repeat_unit(seqs: list[tuple[str, str]],
                        window: int, k_range=(5, 6, 7, 8)) -> str | None:
    """
    Auto-detect the dominant telomere repeat unit from terminal windows.
    Returns the canonical k-mer with highest normalized terminal count, or None.
    """
    terminal_counter: Counter = Counter()
    total_bases = 0

    for _name, seq in seqs:
        ends = []
        w    = min(window, len(seq))
        ends.append(seq[:w])
        if len(seq) > w:
            ends.append(seq[-w:])
        for region in ends:
            total_bases += len(region)
            ru = region.upper()
            for k in k_range:
                for i in range(len(ru) - k + 1):
                    kmer = ru[i:i + k]
                    if set(kmer) <= set("ACGT"):
                        terminal_counter[_canonical(kmer)] += 1

    if not terminal_counter:
        return None

    # Normalise by number of sequences and pick the most frequent canonical k-mer
    top_kmer, top_count = terminal_counter.most_common(1)[0]
    # Require it to be meaningfully enriched (> 0.5 % of terminal bases)
    if top_count / max(total_bases, 1) < 0.005:
        return None
    return top_kmer


# ── Module 0: Telomere identification ─────────────────────────────────────────

def run_module0_telomeres(fasta: Path, repeat_unit: str | None,
                          tel_window: int, tel_density: float,
                          tel_min_len: int,
                          results: Path, workdir: Path,
                          prefix: str, force: bool) -> tuple[Path, str | None]:
    """
    Scan contig ends for telomeric repeats.
    Returns (gff3_path, repeat_unit_used).
    """
    out_gff3 = results / f"mod00_telomeres_{prefix}.gff3"
    if _checkpoint(out_gff3, "telomere-scan", force):
        # Try to recover the repeat unit from existing GFF3
        with open(out_gff3) as fh:
            for line in fh:
                if line.startswith("##repeat-unit"):
                    return out_gff3, line.split()[-1]
        return out_gff3, repeat_unit

    seqs = _read_fasta(fasta)
    _log(f"  Loaded {len(seqs)} sequences")

    if repeat_unit is None:
        _log("  Auto-detecting telomere repeat unit …")
        repeat_unit = _detect_repeat_unit(seqs, tel_window)
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

    for name, seq in seqs:
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
            in_tel = False
            tel_start_rel = 0

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
    """
    Soft-mask (lowercase) low-complexity regions with tantan.
    Returns path to masked FASTA.
    """
    masked = workdir / "masked_soft.fasta"
    if _checkpoint(masked, "tantan-masking", force):
        return masked

    tantan = _require_tool("tantan")
    result = _run([tantan, str(fasta.resolve())], capture_stdout=True)
    masked.write_text(result.stdout)
    _log(f"  Soft-masked FASTA: {masked.name}")

    # Create hard-masked version (lowercase → N) for search tools
    hard = workdir / "masked_hard.fasta"
    _log(f"  Creating hard-masked FASTA for search tools: {hard.name}")
    recs = _read_fasta(masked)
    hard_recs = [(name, re.sub(r"[acgt]", "N", seq)) for name, seq in recs]
    _write_fasta(hard_recs, hard)

    return masked


def _hard_masked(workdir: Path) -> Path:
    return workdir / "masked_hard.fasta"


# ── Module 2: rRNA annotation (Infernal + Rfam) ───────────────────────────────

def _ensure_cms(kingdom: str, rfam_dir: Path | None) -> Path:
    """
    Download (if needed) and return path to a concatenated .cm file
    containing all models for the given kingdom.
    """
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
    """
    Build (if needed) and return path to a concatenated HMMER3 .hmm file for
    the given kingdom.  Strategy: download the Rfam seed Stockholm alignment
    for each family and build an HMM with hmmbuild.  The Rfam REST API does
    not expose standalone .hmm files, but the Stockholm endpoint is stable.
    """
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
    """Parse nhmmer --tblout output into a list of hit dicts.

    nhmmer tblout columns (whitespace-separated):
      0  target name   1  target acc   2  query name   3  query acc
      4  hmmfrom       5  hmmto        6  alifrom      7  alito
      8  envfrom       9  envto        10 sq len       11 strand
      12 E-value       13 score        14 bias         15+ description
    """
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


def run_module2_rrna(fasta: Path, kingdom: str, threads: int, evalue: float,
                     rfam_dir: Path | None, mxsize: float, search_tool: str,
                     workdir: Path, results: Path,
                     prefix: str, force: bool) -> Path:
    """Run nhmmer (default) or cmsearch and write rRNA GFF3."""
    out_gff3  = results / f"mod02_rRNA_{prefix}.gff3"
    tblout    = workdir / f"{search_tool}_rRNA.tblout"

    if _checkpoint(out_gff3, f"{search_tool}-rRNA", force):
        return out_gff3

    search_fa = _hard_masked(workdir)
    if not search_fa.exists():
        search_fa = fasta   # fall back if masking was skipped

    if search_tool == "nhmmer":
        hmm_file = _ensure_hmms(kingdom, rfam_dir)
        tool     = _require_tool("nhmmer")
        _run([
            tool,
            "--cpu",    str(threads),
            "--tblout", str(tblout),
            "-E",       str(evalue),
            "--noali",
            str(hmm_file),
            str(search_fa.resolve()),
        ], cwd=workdir)
        hits = _parse_nhmmer_tblout(tblout, evalue)
    else:
        cm_file  = _ensure_cms(kingdom, rfam_dir)
        tool     = _require_tool("cmsearch")
        _run([
            tool,
            "--cpu",    str(threads),
            "--tblout", str(tblout),
            "-E",       str(evalue),
            "--noali",
            "--rfam",                   # HMM pre-filter; essential for large genomes
            "--mxsize", str(mxsize),    # cap DP matrix per thread (Mb)
            str(cm_file),
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
    return out_gff3


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

    raw_lines = sum(1 for l in open(aragorn_out) if not l.startswith("#") and l.strip())
    _log(f"  ARAGORN output lines: {raw_lines}")
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

def run_module4_integration(tel_gff: Path | None, rrna_gff: Path | None,
                             trna_gff: Path | None,
                             results: Path, prefix: str,
                             genome_size: int = 0) -> tuple[Path, dict]:
    """Merge GFF3 files, sort by position, and write a detailed summary TSV.

    Returns (combined_gff3_path, counts_dict) where counts_dict has keys
    n_tel, n_rrna, n_trna, rrna_by_type {name: count}, trna_by_type {name: count},
    tel_by_end {5prime: n, 3prime: n}.
    """
    from collections import Counter

    combined = results / f"mod04_annotation_{prefix}.gff3"
    summary  = results / f"mod04_summary_{prefix}.tsv"

    # ── Collect and sort all feature lines ────────────────────────────────────
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

    with open(combined, "w") as out:
        out.write(_GFF3_HEADER)
        for line in all_lines:
            out.write(line if line.endswith("\n") else line + "\n")

    _log(f"  Combined GFF3: {combined.name}")

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
    _RRNA_ORDER = [
        "5S_rRNA", "5_8S_rRNA",
        "SSU_rRNA_eukarya", "SSU_rRNA_bacteria", "SSU_rRNA_archaea",
        "LSU_rRNA_eukarya", "LSU_rRNA_bacteria", "LSU_rRNA_archaea",
    ]

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

    # ── Load data ──────────────────────────────────────────────────────────────
    seq_lengths = _fasta_seq_lengths(fasta)
    tel_pos     = _parse_gff3_positions(tel_gff)
    rrna_pos    = _parse_gff3_positions(rrna_gff)
    trna_pos    = _parse_gff3_positions(trna_gff)
    summ        = _parse_summary_tsv(summary_tsv)

    # Select and sort top N sequences
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

    # ── Figure layout ──────────────────────────────────────────────────────────
    ideo_h = max(4, min(n_seqs * 0.28, 14))   # scale height with seq count
    fig    = plt.figure(figsize=(18, ideo_h + 5))
    gs     = GridSpec(2, 3, figure=fig,
                      height_ratios=[ideo_h, 4.5],
                      hspace=0.45, wspace=0.38)
    ax_ideo  = fig.add_subplot(gs[0, :])
    ax_rrna  = fig.add_subplot(gs[1, 0])
    ax_trna  = fig.add_subplot(gs[1, 1])
    ax_donut = fig.add_subplot(gs[1, 2])

    # ── Panel 1: Genome ideogram ───────────────────────────────────────────────
    bar_h   = 0.65
    tel_h   = bar_h + 0.30          # telomeres drawn taller so they stand out
    min_vis = max_len * 0.002       # 0.2 % of longest seq — minimum visible width
    tel_min = max_len * 0.005       # 0.5 % — telomeres get a wider minimum

    # Count total features per type to decide draw order (most abundant = first/bottom)
    n_rrna_total = sum(len(v) for v in rrna_pos.values())
    n_trna_total = sum(len(v) for v in trna_pos.values())
    n_tel_total  = sum(len(v) for v in tel_pos.values())

    # Sort layers: most abundant drawn first (bottom), rarest drawn last (top)
    layers = sorted([
        ("rrna", n_rrna_total, _C_RRNA, 0.40, bar_h),
        ("trna", n_trna_total, _C_TRNA, 0.60, bar_h),
        ("tel",  n_tel_total,  _C_TEL,  1.00, tel_h),
    ], key=lambda x: -x[1])    # descending count → bottom to top

    for i, (name, slen) in enumerate(top_seqs):
        y = n_seqs - i - 1

        # Sequence bar (grey background)
        ax_ideo.broken_barh([(0, slen)], (y - bar_h / 2, bar_h),
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
                ax_ideo.broken_barh(segs, (y - height / 2, height),
                                    facecolors=color, alpha=alpha, linewidth=0)

    ax_ideo.set_ylim(-0.8, n_seqs - 0.2)
    ax_ideo.set_yticks(range(n_seqs))
    ax_ideo.set_yticklabels([n for n, _ in reversed(top_seqs)], fontsize=8)
    ax_ideo.set_xlim(0, max_len * 1.01)
    ax_ideo.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, _: f"{x / 1e6:.0f} Mb"))
    ax_ideo.set_xlabel("Genomic position")
    sort_label = "by SeqID" if sort_by == "seqid" else "longest first"
    title_suffix = (f"top {n_seqs} of {len(seq_lengths)} sequences, {sort_label}"
                    if n_seqs < len(seq_lengths)
                    else f"{n_seqs} sequences, {sort_label}")
    ax_ideo.set_title(f"{prefix}  —  genome annotation overview  ({title_suffix})",
                      fontsize=11, pad=8)
    ax_ideo.spines[["top", "right", "left"]].set_visible(False)
    ax_ideo.tick_params(left=False)

    legend_patches = [
        mpatches.Patch(color=_C_TEL,  label=f"Telomere (n={n_tel_total:,})"),
        mpatches.Patch(color=_C_RRNA, alpha=0.40, label=f"rRNA (n={n_rrna_total:,})"),
        mpatches.Patch(color=_C_TRNA, alpha=0.60, label=f"tRNA (n={n_trna_total:,})"),
    ]
    ax_ideo.legend(handles=legend_patches, loc="lower right",
                   frameon=True, framealpha=0.85, fontsize=9)

    # ── Panel 2: rRNA subtype bars ─────────────────────────────────────────────
    rrna_data = summ.get("rRNA", {})
    rrna_sub  = [(rt, rrna_data[rt][0]) for rt in _RRNA_ORDER if rt in rrna_data]
    if not rrna_sub:   # fallback: any subtypes not in the standard list
        rrna_sub = sorted([(k, v[0]) for k, v in rrna_data.items() if k != "TOTAL"],
                           key=lambda x: -x[1])

    if rrna_sub:
        lbls, vals = zip(*rrna_sub)
        yp = list(range(len(lbls)))
        ax_rrna.barh(yp, vals, color=_C_RRNA, linewidth=0)
        ax_rrna.set_yticks(yp)
        ax_rrna.set_yticklabels(lbls, fontsize=9)
        ax_rrna.set_xlabel("Count")
        ax_rrna.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        for j, v in enumerate(vals):
            ax_rrna.text(v + max(vals) * 0.01, j, f"{v:,}",
                         va="center", fontsize=8)
        ax_rrna.set_xlim(right=max(vals) * 1.18)
    else:
        ax_rrna.text(0.5, 0.5, "No rRNA data", ha="center", va="center",
                     transform=ax_rrna.transAxes, color=_C_OTHER)
    ax_rrna.set_title("rRNA subtypes")
    ax_rrna.spines[["top", "right"]].set_visible(False)

    # ── Panel 3: tRNA type bars (top 20) ──────────────────────────────────────
    trna_data = summ.get("tRNA", {})
    trna_sub  = sorted([(k, v[0]) for k, v in trna_data.items() if k != "TOTAL"],
                        key=lambda x: (-x[1], x[0]))[:20]

    if trna_sub:
        lbls, vals = zip(*trna_sub)
        yp = list(range(len(lbls)))
        ax_trna.barh(yp, vals, color=_C_TRNA, linewidth=0)
        ax_trna.set_yticks(yp)
        ax_trna.set_yticklabels(lbls, fontsize=9)
        ax_trna.set_xlabel("Count")
        ax_trna.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        for j, v in enumerate(vals):
            ax_trna.text(v + max(vals) * 0.01, j, f"{v:,}",
                         va="center", fontsize=8)
        ax_trna.set_xlim(right=max(vals) * 1.18)
    else:
        ax_trna.text(0.5, 0.5, "No tRNA data", ha="center", va="center",
                     transform=ax_trna.transAxes, color=_C_OTHER)
    ax_trna.set_title("tRNA types (top 20)")
    ax_trna.spines[["top", "right"]].set_visible(False)

    # ── Panel 4: Genome composition donut ─────────────────────────────────────
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
        wedges, _ = ax_donut.pie(
            vals, colors=cols,
            wedgeprops=dict(width=0.45, edgecolor="white", linewidth=1.5),
            startangle=90,
        )
        ax_donut.text(0, 0, f"{genome_size / 1e6:.0f} Mb",
                      ha="center", va="center", fontsize=10, fontweight="bold")
        ax_donut.legend(
            wedges,
            [f"{l}  {v / genome_size * 100:.3f}%" for v, l in zip(vals, lbls)],
            loc="lower center", bbox_to_anchor=(0.5, -0.18),
            frameon=False, fontsize=8, ncol=2,
        )
    else:
        ax_donut.text(0.5, 0.5, "No data", ha="center", va="center",
                      transform=ax_donut.transAxes, color=_C_OTHER)
    ax_donut.set_title("Genome composition")

    # ── Save ───────────────────────────────────────────────────────────────────
    for fmt in plot_formats:
        out_path = out_base.with_suffix(f".{fmt}")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        _log(f"  Plot saved: {out_path.name}")
    plt.close(fig)


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

    gen = ap.add_argument_group("General")
    gen.add_argument("--threads", type=int, default=4,
                     help="CPU threads for cmsearch (default: 4)")
    gen.add_argument("--skip_module0", action="store_true",
                     help="Skip Module 0 — telomere identification")
    gen.add_argument("--skip_module1", action="store_true",
                     help="Skip Module 1 — low-complexity masking")
    gen.add_argument("--skip_module2", action="store_true",
                     help="Skip Module 2 — rRNA annotation")
    gen.add_argument("--skip_module3", action="store_true",
                     help="Skip Module 3 — tRNA annotation")
    gen.add_argument("--skip_integration", action="store_true",
                     help="Skip Module 4 — do not produce merged GFF3")
    gen.add_argument("--skip_module5", action="store_true",
                     help="Skip Module 5 — do not produce visualization figure")
    gen.add_argument("--format", default="pdf",
                     help="Plot format(s): pdf, png, svg — comma-separated "
                          "(default: pdf)")
    gen.add_argument("--top_sequences", type=int, default=50,
                     help="Number of sequences to show in the ideogram (default: 50)")
    gen.add_argument("--sort_sequences", choices=["length", "seqid"], default="length",
                     help="Order sequences in the ideogram: 'length' (longest first) "
                          "or 'seqid' (natural alphabetic/numeric sort); default: length")
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

    # Open log
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

    _banner(f"UbboTELORNA {VERSION}")
    _log(f"  Input FASTA : {args.fasta}")
    _log(f"  Output dir  : {run_dir}")
    _log(f"  Kingdom     : {args.kingdom}")
    _log(f"  rRNA tool   : {args.search_tool}")
    _log(f"  Threads     : {args.threads}")
    genome_size = _fasta_total_length(args.fasta)
    _log(f"  Genome size : {genome_size:,} bp")

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
        if not args.skip_module0:
            _log("    [0] Telomere identification  →  results/mod00_telomeres_*.gff3")
        if not args.skip_module1:
            _log("    [1] Low-complexity masking   →  workdir/masked_soft.fasta")
        if not args.skip_module2:
            _log(f"    [2] rRNA annotation ({args.search_tool})  →  results/mod02_rRNA_*.gff3")
        if not args.skip_module3:
            _log("    [3] tRNA annotation          →  results/mod03_tRNA_*.gff3")
        if not args.skip_integration:
            _log("    [4] Integration              →  results/mod04_annotation_*.gff3")
        if not args.skip_module5:
            fmt0 = args.format.split(",")[0].strip()
            _log(f"    [5] Visualization ({args.sort_sequences})  →  results/mod05_plot_*.{fmt0}")
        _log("  Exiting (--dry_run).")
        if _LOG_FH:
            _LOG_FH.close()
        sys.exit(0)

    # Carbon tracking
    t_start  = time.monotonic()
    _tracker = None
    if args.disable_co2_tracking:
        _log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
    else:
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
                _shim.get_distribution    = _get_dist
                _shim.resource_filename   = _resource_filename
                _shim.DistributionNotFound = Exception
                sys.modules["pkg_resources"] = _shim

            from codecarbon import EmissionsTracker
            _tracker = EmissionsTracker(
                output_dir=str(logs_dir),
                output_file=f"{prefix}.emissions.csv",
                project_name="UbboTELORNA",
                log_level="warning",
            )
            _tracker.start()
            _log("  codecarbon tracker started")
        except ImportError as e:
            _log(f"  codecarbon not installed — carbon tracking skipped ({e})")
        except Exception as e:
            _log(f"  codecarbon failed to start — carbon tracking skipped ({e})")

    # ── Module 0: Telomere identification ─────────────────────────────────────
    tel_gff     = results / f"mod00_telomeres_{prefix}.gff3"
    repeat_used = None
    if not args.skip_module0:
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

    # ── Module 1: Low-complexity masking ──────────────────────────────────────
    if not args.skip_module1:
        _banner("Module 1 — Low-Complexity Masking")
        run_module1_masking(
            fasta   = args.fasta,
            workdir = workdir,
            force   = args.force,
        )

    # ── Module 2: rRNA annotation ─────────────────────────────────────────────
    rrna_gff = results / f"mod02_rRNA_{prefix}.gff3"
    if not args.skip_module2:
        _banner("Module 2 — rRNA Annotation")
        rrna_gff = run_module2_rrna(
            fasta       = args.fasta,
            kingdom     = args.kingdom,
            threads     = args.threads,
            evalue      = args.evalue,
            rfam_dir    = args.rfam_dir,
            mxsize      = args.cmsearch_mxsize,
            search_tool = args.search_tool,
            workdir     = workdir,
            results     = results,
            prefix      = prefix,
            force       = args.force,
        )

    # ── Module 3: tRNA annotation ─────────────────────────────────────────────
    trna_gff = results / f"mod03_tRNA_{prefix}.gff3"
    if not args.skip_module3:
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
    if not args.skip_integration:
        _banner("Module 4 — Integration")
        _, counts = run_module4_integration(
            tel_gff, rrna_gff, trna_gff, results, prefix,
            genome_size=genome_size)

    # ── Module 5: Visualization ───────────────────────────────────────────────
    plot_formats = [f.strip().lstrip(".") for f in args.format.split(",")]
    if not args.skip_module5:
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
