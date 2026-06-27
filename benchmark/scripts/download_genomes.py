#!/usr/bin/env python3
"""Download a genome assembly (FASTA + GFF3) from NCBI RefSeq.

Uses ncbi-genome-download so the download logic is handled by a maintained
library rather than manual FTP URL construction.

Usage (called by Snakemake rules):
    python3 download_genomes.py \\
        --accession GCF_000005845.2 \\
        --group     bacteria \\
        --outdir    results/genomes/ecoli_k12 \\
        --genome    ecoli_k12
"""

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def _decompress_gz(src: Path, dst: Path) -> None:
    with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def _find_single(directory: Path, pattern: str) -> Path:
    matches = list(directory.glob(pattern))
    if not matches:
        print(f"ERROR: no file matching '{pattern}' in {directory}", file=sys.stderr)
        sys.exit(1)
    return matches[0]


def download(accession: str, group: str, outdir: Path, genome: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    fasta_out = outdir / f"{genome}.fasta"
    gff3_out  = outdir / f"{genome}_ref.gff3"

    if fasta_out.exists() and gff3_out.exists() \
            and fasta_out.stat().st_size > 0 and gff3_out.stat().st_size > 0:
        print(f"[download] {genome}: already downloaded, skipping")
        return

    dl_dir = outdir / "ncbi_dl"
    dl_dir.mkdir(exist_ok=True)

    print(f"[download] {genome}: fetching {accession} ({group}) …")
    cmd = [
        "ncbi-genome-download",
        "--formats", "fasta,gff",
        "--assembly-accessions", accession,
        "--output-folder", str(dl_dir),
        "--flat-output",
        "--no-cache",
        group,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: ncbi-genome-download failed:\n{result.stderr[-3000:]}",
              file=sys.stderr)
        sys.exit(1)

    # ncbi-genome-download writes compressed files into dl_dir with flat layout
    fasta_gz = _find_single(dl_dir, "*_genomic.fna.gz")
    gff_gz   = _find_single(dl_dir, "*_genomic.gff.gz")

    print(f"[download] {genome}: decompressing FASTA …")
    _decompress_gz(fasta_gz, fasta_out)

    print(f"[download] {genome}: decompressing GFF3 …")
    _decompress_gz(gff_gz, gff3_out)

    # clean up the bulky download directory
    shutil.rmtree(dl_dir)
    print(f"[download] {genome}: done — {fasta_out.stat().st_size // 1_000_000} MB FASTA, "
          f"{gff3_out.stat().st_size // 1_000} KB GFF3")


def main() -> None:
    ap = argparse.ArgumentParser(description="Download NCBI genome FASTA + GFF3")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--group",     required=True,
                    help="ncbi-genome-download group (bacteria, plant, …)")
    ap.add_argument("--outdir",    required=True, type=Path)
    ap.add_argument("--genome",    required=True,
                    help="Short genome key used for output filenames")
    args = ap.parse_args()
    download(args.accession, args.group, args.outdir, args.genome)


if __name__ == "__main__":
    main()
