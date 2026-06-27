#!/usr/bin/env python3
"""Create adversarial and fragmented genome FASTAs for robustness testing.

Adversarial mode:
    Prepends N copies of a telomere repeat unit to the start and end of
    every sequence, then reverses the repeat for the end.  This mimics
    an assembly where chromosomes begin and end with telomeric arrays —
    the exact condition that causes barrnap to fail.

Fragmented mode:
    Splits every sequence into non-overlapping chunks of a fixed size,
    producing a highly fragmented assembly that tests tool behaviour on
    short scaffolds.
"""

import argparse
import sys
from pathlib import Path


def _read_fasta(path: Path):
    """Yield (header, seq) pairs from a FASTA file."""
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line, []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def _write_seq(fh, header: str, seq: str, width: int = 80) -> None:
    fh.write(header + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def _rev(seq: str) -> str:
    return seq[::-1]


def make_adversarial(fasta: Path, repeat: str, n_copies: int,
                     output: Path) -> None:
    """Prepend and append a telomere repeat array to every sequence."""
    prefix = repeat * n_copies
    suffix = _rev(repeat) * n_copies          # reverse complement not needed
    n_seq = 0
    with open(output, "w") as fout:
        for header, seq in _read_fasta(fasta):
            modified = prefix + seq + suffix
            _write_seq(fout, header + " [adversarial_telomere]", modified)
            n_seq += 1
    print(f"[adversarial] {n_seq} sequences modified: "
          f"+{len(prefix)} bp prefix, +{len(suffix)} bp suffix "
          f"({n_copies} × '{repeat}')", file=sys.stderr)


def make_fragmented(fasta: Path, fragment_size: int, output: Path) -> None:
    """Split every sequence into non-overlapping fragments."""
    n_seq = n_frag = 0
    with open(output, "w") as fout:
        for header, seq in _read_fasta(fasta):
            name = header.lstrip(">").split()[0]
            n_seq += 1
            for i, start in enumerate(range(0, len(seq), fragment_size)):
                chunk = seq[start:start + fragment_size]
                if not chunk:
                    continue
                _write_seq(fout,
                           f">{name}_frag{i + 1}_{start + 1}_{start + len(chunk)}",
                           chunk)
                n_frag += 1
    print(f"[fragmented] {n_seq} sequences → {n_frag} fragments "
          f"({fragment_size} bp chunks)", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create adversarial / fragmented FASTA for robustness tests")
    ap.add_argument("--fasta",         required=True, type=Path)
    ap.add_argument("--output",        required=True, type=Path)
    ap.add_argument("--repeat",        default=None,
                    help="Telomere repeat unit for adversarial mode")
    ap.add_argument("--n_copies",      type=int, default=1000,
                    help="Copies of repeat to prepend/append (default 1000)")
    ap.add_argument("--fragment_size", type=int, default=None,
                    help="Fragment size in bp for fragmented mode")
    args = ap.parse_args()

    if args.repeat and args.fragment_size:
        print("ERROR: specify --repeat OR --fragment_size, not both",
              file=sys.stderr)
        sys.exit(1)
    if not args.repeat and not args.fragment_size:
        print("ERROR: specify one of --repeat or --fragment_size", file=sys.stderr)
        sys.exit(1)

    if args.repeat:
        make_adversarial(args.fasta, args.repeat, args.n_copies, args.output)
    else:
        make_fragmented(args.fasta, args.fragment_size, args.output)


if __name__ == "__main__":
    main()
