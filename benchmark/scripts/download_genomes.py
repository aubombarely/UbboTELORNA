#!/usr/bin/env python3
"""Download a genome assembly (FASTA + GFF3) from NCBI RefSeq.

Uses the NCBI datasets CLI (ncbi-datasets-cli conda package) instead of
ncbi-genome-download, which breaks whenever NCBI adds a column to their
assembly summary TSV.

Usage (called by Snakemake rules):
    python3 download_genomes.py \\
        --accession GCF_000005845.2 \\
        --outdir    results/genomes/ecoli_k12 \\
        --genome    ecoli_k12

The --group argument is accepted but ignored (datasets CLI resolves groups
automatically from the accession, so it is kept only for Snakemake
compatibility with the old interface).
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

_MAX_ATTEMPTS = 3
_RETRY_DELAY  = 60  # seconds between retries


def _find_single(directory: Path, pattern: str) -> Path | None:
    matches = list(directory.glob(pattern))
    return matches[0] if matches else None


def _require_datasets() -> str:
    tool = shutil.which("datasets")
    if tool is None:
        print(
            "ERROR: 'datasets' not found in PATH.\n"
            "       Install with: conda install -c conda-forge ncbi-datasets-cli",
            file=sys.stderr,
        )
        sys.exit(1)
    return tool


def download(accession: str, outdir: Path, genome: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_out = outdir / f"{genome}.fasta"
    gff3_out  = outdir / f"{genome}_ref.gff3"

    if (fasta_out.exists() and gff3_out.exists()
            and fasta_out.stat().st_size > 0
            and gff3_out.stat().st_size > 0):
        print(f"[download] {genome}: already downloaded, skipping")
        return

    datasets = _require_datasets()
    print(f"[download] {genome}: fetching {accession} via NCBI datasets CLI …")

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = Path(tmpdir) / "download.zip"

        cmd = [
            datasets, "download", "genome", "accession", accession,
            "--include", "genome,gff3",
            "--filename", str(zip_path),
            "--no-progressbar",
        ]
        last_err = ""
        for attempt in range(1, _MAX_ATTEMPTS + 1):
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                break
            last_err = result.stderr[-3000:]
            if attempt < _MAX_ATTEMPTS:
                print(
                    f"[download] {genome}: attempt {attempt}/{_MAX_ATTEMPTS} failed; "
                    f"retrying in {_RETRY_DELAY}s …",
                    file=sys.stderr,
                )
                time.sleep(_RETRY_DELAY)
        else:
            print(
                f"ERROR: datasets download failed for {accession} "
                f"after {_MAX_ATTEMPTS} attempts:\n{last_err}",
                file=sys.stderr,
            )
            sys.exit(1)

        if not zip_path.exists() or zip_path.stat().st_size == 0:
            print(
                f"ERROR: datasets returned success but wrote no zip file for {accession}.\n"
                f"       The assembly record may exist but have no downloadable data.\n"
                f"       Check: datasets summary genome accession {accession}",
                file=sys.stderr,
            )
            sys.exit(1)

        extract_dir = Path(tmpdir) / "extracted"
        extract_dir.mkdir()
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)

        # NCBI datasets writes files under ncbi_dataset/data/<accession>/
        data_root = extract_dir / "ncbi_dataset" / "data"
        data_dir  = data_root / accession
        if not data_dir.exists():
            # Accession folder may have version suffix stripped; try glob
            candidates = sorted(data_root.glob("GC*"))
            if not candidates:
                print(
                    f"ERROR: datasets zip for {accession} contained no genome data.\n"
                    f"       The assembly may be suppressed or not packaged for download.\n"
                    f"       Find a replacement with:\n"
                    f"         datasets summary genome taxon '<organism>' --reference",
                    file=sys.stderr,
                )
                sys.exit(1)
            data_dir = candidates[0]

        fasta_src = _find_single(data_dir, "*.fna")
        gff3_src  = _find_single(data_dir, "*.gff")

        if fasta_src is None:
            print(f"ERROR: no .fna file in {data_dir}", file=sys.stderr)
            sys.exit(1)
        if gff3_src is None:
            print(f"ERROR: no .gff file in {data_dir}", file=sys.stderr)
            sys.exit(1)

        shutil.copy2(fasta_src, fasta_out)
        shutil.copy2(gff3_src,  gff3_out)

    print(
        f"[download] {genome}: done — "
        f"{fasta_out.stat().st_size // 1_000_000} MB FASTA, "
        f"{gff3_out.stat().st_size // 1_000} KB GFF3"
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Download NCBI genome FASTA + GFF3 via the datasets CLI"
    )
    ap.add_argument("--accession", required=True,
                    help="NCBI RefSeq accession (e.g. GCF_000005845.2)")
    ap.add_argument("--group",     required=False, default=None,
                    help="Ignored — kept for Snakemake compatibility")
    ap.add_argument("--outdir",    required=True, type=Path,
                    help="Output directory; FASTA and GFF3 written here")
    ap.add_argument("--genome",    required=True,
                    help="Short genome key used for output filenames")
    args = ap.parse_args()
    download(args.accession, args.outdir, args.genome)


if __name__ == "__main__":
    main()
