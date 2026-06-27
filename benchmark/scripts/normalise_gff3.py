#!/usr/bin/env python3
"""Normalise tool-specific GFF3 output to a common format.

Each external tool uses slightly different feature types, Name attribute
values, and attribute vocabularies.  This script converts all of them to
the internal naming used by UbboTELORNA so that compare_annotations.py
can work against a single reference.

Supported tools: barrnap, trnascan
UbboTELORNA output is already in the canonical format and does not need
normalisation.

Canonical rRNA Name values (matching Rfam accession names):
    5S_rRNA  5_8S_rRNA  SSU_rRNA_eukarya  LSU_rRNA_eukarya
    SSU_rRNA_bacteria  LSU_rRNA_bacteria
    SSU_rRNA_archaea   LSU_rRNA_archaea

Canonical tRNA Name values:  tRNA-{AminoAcid}  (e.g. tRNA-Phe, tRNA-Ser)
"""

import argparse
import re
import sys
from pathlib import Path


# barrnap Name → canonical Name
_BARRNAP_RRNA = {
    "5S_rRNA":   "5S_rRNA",
    "5.8S_rRNA": "5_8S_rRNA",
    "18S_rRNA":  "SSU_rRNA_eukarya",
    "28S_rRNA":  "LSU_rRNA_eukarya",
    "16S_rRNA":  "SSU_rRNA_bacteria",
    "23S_rRNA":  "LSU_rRNA_bacteria",
    "12S_rRNA":  "SSU_rRNA_mitochondria",
    "16S_rRNA_mito": "LSU_rRNA_mitochondria",
}

# RefSeq GFF3 product strings → canonical Name (used to parse reference GFF3)
_REFSEQ_RRNA = {
    "5S ribosomal RNA":              "5S_rRNA",
    "5.8S ribosomal RNA":            "5_8S_rRNA",
    "18S ribosomal RNA":             "SSU_rRNA_eukarya",
    "28S ribosomal RNA":             "LSU_rRNA_eukarya",
    "16S ribosomal RNA":             "SSU_rRNA_bacteria",
    "23S ribosomal RNA":             "LSU_rRNA_bacteria",
    "small subunit ribosomal RNA":   "SSU_rRNA_eukarya",
    "large subunit ribosomal RNA":   "LSU_rRNA_eukarya",
}


def _parse_attrs(attr_str: str) -> dict:
    """Parse GFF3 attribute column into a dict."""
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _fmt_attrs(attrs: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


def _normalise_barrnap(cols: list) -> list | None:
    """Return normalised columns or None to skip the row."""
    if len(cols) < 9 or cols[2] != "rRNA":
        return None
    attrs = _parse_attrs(cols[8])
    raw_name = attrs.get("Name", "")
    canonical = _BARRNAP_RRNA.get(raw_name)
    if canonical is None:
        # Try prefix match (barrnap sometimes appends extra text)
        for k, v in _BARRNAP_RRNA.items():
            if raw_name.startswith(k):
                canonical = v
                break
    if canonical is None:
        canonical = raw_name          # pass through unknown names unchanged
    attrs["Name"]  = canonical
    attrs["tool"]  = "barrnap"
    cols[8] = _fmt_attrs(attrs)
    return cols


def _normalise_trnascan(cols: list) -> list | None:
    """Return normalised columns or None to skip the row."""
    if len(cols) < 9 or cols[2] not in ("tRNA", "pseudogene"):
        return None
    attrs = _parse_attrs(cols[8])
    # tRNAscan-SE writes Name=tRNA-Phe  or  Name=Pseudo
    name = attrs.get("Name", "")
    if name.startswith("tRNA-"):
        pass                          # already canonical
    elif "amino_acid" in attrs:
        aa = attrs["amino_acid"]
        name = f"tRNA-{aa}"
    # strip anticodon from Name if present: tRNA-Phe(GAA) → tRNA-Phe
    name = re.sub(r"\([A-Za-z]{3}\)$", "", name)
    attrs["Name"] = name
    attrs["tool"] = "trnascan"
    cols[2] = "tRNA"                  # normalise "pseudogene" → "tRNA"
    cols[8] = _fmt_attrs(attrs)
    return cols


_NORMALISERS = {
    "barrnap":  _normalise_barrnap,
    "trnascan": _normalise_trnascan,
}


def normalise(tool: str, input_path: Path, output_path: Path) -> None:
    fn = _NORMALISERS.get(tool)
    if fn is None:
        print(f"ERROR: unknown tool '{tool}'", file=sys.stderr)
        sys.exit(1)

    written = 0
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            result = fn(cols)
            if result is not None:
                fout.write("\t".join(result) + "\n")
                written += 1

    print(f"[normalise] {tool}: {written} features written to {output_path}",
          file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description="Normalise GFF3 to UbboTELORNA format")
    ap.add_argument("--tool",   required=True, choices=list(_NORMALISERS))
    ap.add_argument("--input",  required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    normalise(args.tool, args.input, args.output)


if __name__ == "__main__":
    main()
