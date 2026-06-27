#!/usr/bin/env python3
"""
UbboTELORNA.py  —  Telomere, rRNA, and tRNA annotation for genome assemblies.

Annotates three classes of ancient, conserved genomic elements:

  Module 0  Telomere identification  — k-mer density scan at contig ends
  Module 1  Low-complexity masking   — tantan soft-mask (prevents search failures
                                       on telomeric and repetitive sequences)
  Module 2  rRNA annotation          — Infernal cmsearch + Rfam covariance models
  Module 3  tRNA annotation          — ARAGORN
  Module 4  Integration              — merged GFF3 + summary table

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
    _run([tantan, "-m", "0", str(fasta)],
         capture_stdout=True, cwd=workdir)

    # tantan writes to stdout; capture and save
    result = subprocess.run(
        [tantan, "-m", "0", str(fasta)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        print(f"ERROR: tantan failed:\n{result.stderr[-2000:]}", file=sys.stderr)
        sys.exit(1)
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
                     rfam_dir: Path | None, workdir: Path, results: Path,
                     prefix: str, force: bool) -> Path:
    """Run cmsearch and write rRNA GFF3."""
    out_gff3  = results / f"mod02_rRNA_{prefix}.gff3"
    tblout    = workdir / "cmsearch_rRNA.tblout"
    cmsearch_out = workdir / "cmsearch_rRNA.out"

    if _checkpoint(out_gff3, "cmsearch-rRNA", force):
        return out_gff3

    cm_file   = _ensure_cms(kingdom, rfam_dir)
    search_fa = _hard_masked(workdir)
    if not search_fa.exists():
        search_fa = fasta   # fall back if masking was skipped

    cmsearch  = _require_tool("cmsearch")
    _run([
        cmsearch,
        "--cpu",    str(threads),
        "--tblout", str(tblout),
        "-E",       str(evalue),
        "--noali",
        str(cm_file),
        str(search_fa),
    ], cwd=workdir)

    hits = _parse_cmsearch_tblout(tblout, evalue)
    _log(f"  cmsearch hits (E ≤ {evalue}): {len(hits)}")

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
    r"\s+\d+\s+((?:tRNA|tmRNA|mtRNA|pseudo_tRNA)-\S+)\s+"
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

    aragorn     = _require_tool("aragorn")
    search_fa   = _hard_masked(workdir)
    if not search_fa.exists():
        search_fa = fasta

    _run([
        aragorn,
        "-t",           # tRNA genes only
        "-gcstd",       # standard genetic code
        "-l",           # treat each sequence as a linear molecule
        "-w",           # show sequence information
        "-o", str(aragorn_out),
        str(search_fa),
    ])

    hits = _parse_aragorn(aragorn_out)
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
                             results: Path, prefix: str) -> Path:
    """Merge GFF3 files and write a summary TSV."""
    combined = results / f"mod04_annotation_{prefix}.gff3"
    summary  = results / f"mod04_summary_{prefix}.tsv"

    n_tel  = 0
    n_rrna = 0
    n_trna = 0

    with open(combined, "w") as out:
        out.write(_GFF3_HEADER)
        for gff, label, counter_name in [
            (tel_gff,  "telomere", "n_tel"),
            (rrna_gff, "rRNA",     "n_rrna"),
            (trna_gff, "tRNA",     "n_trna"),
        ]:
            if gff is None or not gff.exists():
                continue
            count = 0
            with open(gff) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    out.write(line)
                    count += 1
            if label == "telomere":
                n_tel = count
            elif label == "rRNA":
                n_rrna = count
            else:
                n_trna = count

    _log(f"  Combined GFF3: {combined.name}")

    with open(summary, "w") as fh:
        fh.write("feature_type\tcount\tsource_file\n")
        for label, count, src in [
            ("telomere", n_tel,  str(tel_gff)  if tel_gff  else "skipped"),
            ("rRNA",     n_rrna, str(rrna_gff) if rrna_gff else "skipped"),
            ("tRNA",     n_trna, str(trna_gff) if trna_gff else "skipped"),
        ]:
            fh.write(f"{label}\t{count}\t{src}\n")

    _log(f"  Summary TSV: {summary.name}")
    _log(f"  Totals  →  telomeres: {n_tel}  |  rRNA: {n_rrna}  |  tRNA: {n_trna}")
    return combined


# ── Argument parser ───────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="UbboTELORNA",
        description=(
            "Annotate telomeres, rRNA, and tRNA in a genome assembly.\n\n"
            "Modules:\n"
            "  0  Telomere identification  (k-mer scan at contig ends)\n"
            "  1  Low-complexity masking   (tantan)\n"
            "  2  rRNA annotation          (Infernal cmsearch + Rfam CMs)\n"
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
    warnings.filterwarnings("ignore", category=FutureWarning, module="pynvml")
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
    _log(f"  Threads     : {args.threads}")

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
            _log("    [2] rRNA annotation          →  results/mod02_rRNA_*.gff3")
        if not args.skip_module3:
            _log("    [3] tRNA annotation          →  results/mod03_tRNA_*.gff3")
        if not args.skip_integration:
            _log("    [4] Integration              →  results/mod04_annotation_*.gff3")
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
    tel_gff    = None
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
    rrna_gff = None
    if not args.skip_module2:
        _banner("Module 2 — rRNA Annotation")
        rrna_gff = run_module2_rrna(
            fasta    = args.fasta,
            kingdom  = args.kingdom,
            threads  = args.threads,
            evalue   = args.evalue,
            rfam_dir = args.rfam_dir,
            workdir  = workdir,
            results  = results,
            prefix   = prefix,
            force    = args.force,
        )

    # ── Module 3: tRNA annotation ─────────────────────────────────────────────
    trna_gff = None
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
    if not args.skip_integration:
        _banner("Module 4 — Integration")
        run_module4_integration(tel_gff, rrna_gff, trna_gff, results, prefix)

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
        "kingdom":        args.kingdom,
        "telomere_repeat": repeat_used,
        "parameters": {
            "telomere_window":  args.telomere_window,
            "telomere_density": args.telomere_density,
            "telomere_min_len": args.telomere_min_len,
            "evalue":           args.evalue,
            "threads":          args.threads,
        },
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
