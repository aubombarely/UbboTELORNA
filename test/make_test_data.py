#!/usr/bin/env python3
"""
make_test_data.py  —  Generate synthetic test genome for UbboTELORNA.

Produces test/test_genome.fasta with four sequences:

  seq1_telomeric    Plant telomere repeats (TTTAGGG) at both ends
  seq2_rRNA         Eukaryotic 5S rRNA fragment (RF00001-detectable)
  seq3_tRNA         tRNA-Phe(GAA) (ARAGORN-detectable)
  seq4_lowcomplex   AC-repeat low-complexity region (tests masking module)

Usage
-----
    python3 test/make_test_data.py

Outputs
-------
    test/test_genome.fasta   — 4 sequences, ~6 kb total
"""

import random
from pathlib import Path

random.seed(42)

# ── Embedded biological sequences ─────────────────────────────────────────────

# Human 5S rRNA (121 bp) — matches Rfam RF00001 covariance model
_5S_RRNA = (
    "GCCTACGGCCATAACCCTGACCCTGCAGCCAATCTGCGTAAACGAATGGAG"
    "AGTTTGAGTCTGGCCGATCTGAGACCGAAGCTTCCGGGGTCAGGCGGAGCC"
    "TGTCTGAGCGATCT"
)

# Canonical eukaryotic tRNA-Phe(GAA) (73 bp) — ARAGORN-detectable cloverleaf
_TRNA_PHE = (
    "GCGGATTTAGCTCAGTGGTAGAGCACTTGCATGGCATGCAAGAGGACGGGG"
    "TTCGAATCCCGTAATCCGCCA"
)

_PLANT_TEL = "TTTAGGG"


def rand_seq(n: int) -> str:
    return "".join(random.choices("ACGT", k=n))


def write_fasta(path: Path, records: dict) -> None:
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ── Build sequences ───────────────────────────────────────────────────────────

records = {
    # ~2.6 kb: 30 telomere repeats (210 bp) + 2000 bp body + 30 repeats
    "seq1_telomeric": (
        _PLANT_TEL * 30 + rand_seq(2000) + _PLANT_TEL * 30
    ),
    # ~1.1 kb: random flanks + embedded 5S rRNA
    "seq2_rRNA": (
        rand_seq(500) + _5S_RRNA + rand_seq(500)
    ),
    # ~0.7 kb: random flanks + embedded tRNA-Phe
    "seq3_tRNA": (
        rand_seq(300) + _TRNA_PHE + rand_seq(300)
    ),
    # ~1.8 kb: AC-repeat at both ends (tests masking robustness)
    "seq4_lowcomplex": (
        "AC" * 200 + rand_seq(1000) + "AC" * 200
    ),
}

out = Path(__file__).parent / "test_genome.fasta"
write_fasta(out, records)

print(f"Written: {out}")
for name, seq in records.items():
    print(f"  {name}: {len(seq):,} bp")
print(f"  Total: {sum(len(s) for s in records.values()):,} bp")
