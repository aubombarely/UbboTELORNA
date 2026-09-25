#!/usr/bin/env python3
"""Compute sensitivity, precision, and F1 by comparing two GFF3 files.

Standard mode (--feature rRNA or tRNA):
    Predicted GFF3 vs. reference GFF3 using reciprocal overlap.
    Outputs one TSV row per feature subtype plus an "ALL" summary row.

Telomere mode (--feature telomere):
    Compares UbboTELORNA GFF3 telomere features against tidk TSV output.
    No reference ground truth is available for most genomes, so agreement
    between the two tools is reported instead of TP/FP/FN.

Robustness mode (--robustness_mode):
    Compares native vs. adversarial vs. fragmented GFF3 counts for both
    UbboTELORNA and barrnap side by side.  Produces a summary TSV.
"""

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from pathlib import Path


# ── GFF3 parsing ─────────────────────────────────────────────────────────────

def _parse_attrs(s: str) -> dict:
    attrs = {}
    for field in s.strip().split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _load_gff3(path: Path, feature_type: str) -> dict:
    """Return {seqname: [(start, end, strand, subtype, array_id), ...]} sorted by start.

    array_id is the predicted feature's UbboTELORNA --flag_5s_arrays tag
    (array_member=true;array_id=arrXXXX) when present, else None. Array
    members are NOT collapsed here — each is kept as its own feature so it
    can be matched individually against the reference (a real array with
    N genuine copies should earn up to N real TPs when the reference
    actually annotates that many). The over-counting problem this
    replaced (a correctly-detected array being penalized as N-1 separate
    false positives against a sparse reference) is instead handled in
    _match_features, which caps *unmatched* array members at one FP per
    array rather than one per member. See _match_features docstring.
    """
    records = defaultdict(list)
    feature_type_lc = feature_type.lower()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2].lower()
            if not (ftype == feature_type_lc
                    or (feature_type_lc == "telomere" and ftype == "telomere_region")
                    or (feature_type_lc == "rrna" and ftype == "rrna")
                    or (feature_type_lc == "trna" and ftype in ("trna", "ncrna"))):
                continue
            seqname = cols[0]
            start   = int(cols[3])
            end     = int(cols[4])
            strand  = cols[6]
            attrs   = _parse_attrs(cols[8])
            # resolve subtype
            if feature_type_lc == "rrna":
                subtype = (attrs.get("Name")
                           or attrs.get("gene")
                           or attrs.get("product", "rRNA_unknown"))
            elif feature_type_lc == "trna":
                name = attrs.get("Name", "")
                import re
                name = re.sub(r"\([A-Za-z]{3}\)$", "", name)
                subtype = name or "tRNA_unknown"
            else:
                subtype = feature_type_lc

            array_id = attrs.get("array_id") if feature_type_lc == "rrna" else None
            records[seqname].append((start, end, strand, subtype, array_id))

    for seqname in records:
        records[seqname].sort(key=lambda x: x[0])
    return dict(records)


# ── overlap logic ─────────────────────────────────────────────────────────────

def _reciprocal_overlap(a_start, a_end, b_start, b_end, threshold: float) -> bool:
    """True if overlap ≥ threshold of the shorter feature."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap <= 0:
        return False
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)
    return overlap / shorter >= threshold


def _match_features(predicted: dict, reference: dict,
                    threshold: float) -> tuple[int, int, int, dict, dict, dict]:
    """Return (TP, FP, FN, tp_by_sub, fp_by_sub, fn_by_sub).

    Each reference feature can be claimed by at most one predicted feature
    (ref_matched gates re-use), so TP never exceeds the true reference
    count even when many predicted array members overlap the same region.

    Array members (predicted features sharing an array_id — see
    _load_gff3) are matched individually, so a real array with N genuine
    copies earns up to N real TPs when the reference annotates that many.
    Only *unmatched* members are treated specially: instead of counting
    each one as its own FP (which would punish a correctly-detected array
    as N-1 separate errors whenever the reference is sparser than the
    real biology), all of an array's unmatched members together count as
    at most one FP — "this array has more copies than the reference
    lists" is one finding, not N-1 of them. Non-array predictions are
    unaffected: each unmatched one is still its own FP as before.
    """
    tp_by_sub: dict[str, int] = defaultdict(int)
    fp_by_sub: dict[str, int] = defaultdict(int)
    fn_by_sub: dict[str, int] = defaultdict(int)

    all_seqs = set(predicted) | set(reference)

    for seq in all_seqs:
        pred_list = predicted.get(seq, [])
        ref_list  = reference.get(seq, [])

        ref_starts   = [r[0] for r in ref_list]
        ref_matched  = [False] * len(ref_list)

        unmatched_array_sub: dict[str, str] = {}  # array_id -> subtype, for any array with >=1 unmatched member

        for (p_start, p_end, p_strand, p_sub, p_array_id) in pred_list:
            # binary search: find first ref feature that could overlap
            lo = bisect.bisect_left(ref_starts, p_start - (p_end - p_start + 1))
            matched = False
            for i in range(lo, len(ref_list)):
                r_start, r_end, r_strand, r_sub, _r_array_id = ref_list[i]
                if r_start > p_end:
                    break
                if not ref_matched[i] and _reciprocal_overlap(p_start, p_end, r_start, r_end, threshold):
                    matched = True
                    ref_matched[i] = True
                    tp_by_sub[p_sub] += 1
                    break
            if not matched:
                if p_array_id:
                    unmatched_array_sub.setdefault(p_array_id, p_sub)
                else:
                    fp_by_sub[p_sub] += 1

        for sub in unmatched_array_sub.values():
            fp_by_sub[sub] += 1

        for matched, (r_start, r_end, r_strand, r_sub, _r_array_id) in zip(ref_matched, ref_list):
            if not matched:
                fn_by_sub[r_sub] += 1

    tp = sum(tp_by_sub.values())
    fp = sum(fp_by_sub.values())
    fn = sum(fn_by_sub.values())
    return tp, fp, fn, dict(tp_by_sub), dict(fp_by_sub), dict(fn_by_sub)


def _metrics(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    sensitivity = tp / (tp + fn)  if (tp + fn) > 0 else 0.0
    precision   = tp / (tp + fp)  if (tp + fp) > 0 else 0.0
    f1 = (2 * sensitivity * precision / (sensitivity + precision)
          if (sensitivity + precision) > 0 else 0.0)
    return round(sensitivity, 4), round(precision, 4), round(f1, 4)


# ── standard correctness mode ─────────────────────────────────────────────────

def run_standard(args) -> None:
    predicted = _load_gff3(Path(args.predicted), args.feature)
    reference = _load_gff3(Path(args.reference), args.feature)

    tp, fp, fn, tp_sub, fp_sub, fn_sub = _match_features(
        predicted, reference, args.threshold)

    all_subtypes = sorted(set(list(tp_sub) + list(fp_sub) + list(fn_sub)))

    rows = []
    for sub in all_subtypes:
        s_tp = tp_sub.get(sub, 0)
        s_fp = fp_sub.get(sub, 0)
        s_fn = fn_sub.get(sub, 0)
        sens, prec, f1 = _metrics(s_tp, s_fp, s_fn)
        rows.append({
            "genome":    args.genome,
            "tool":      args.tool,
            "feature":   args.feature,
            "subtype":   sub,
            "TP":        s_tp,
            "FP":        s_fp,
            "FN":        s_fn,
            "n_ref":     s_tp + s_fn,
            "n_pred":    s_tp + s_fp,
            "sensitivity": sens,
            "precision":   prec,
            "F1":          f1,
        })

    # overall summary row
    sens_all, prec_all, f1_all = _metrics(tp, fp, fn)
    rows.append({
        "genome": args.genome, "tool": args.tool,
        "feature": args.feature, "subtype": "ALL",
        "TP": tp, "FP": fp, "FN": fn,
        "n_ref": tp + fn, "n_pred": tp + fp,
        "sensitivity": sens_all, "precision": prec_all, "F1": f1_all,
    })

    fields = ["genome","tool","feature","subtype",
              "TP","FP","FN","n_ref","n_pred",
              "sensitivity","precision","F1"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[compare] {args.genome} {args.tool} {args.feature}: "
          f"TP={tp} FP={fp} FN={fn} F1={f1_all}", file=sys.stderr)


# ── telomere comparison (tool-vs-tool) ────────────────────────────────────────

def _load_tidk(path: Path) -> dict:
    """Load tidk TSV (chromosome, forward, reverse counts).
    Returns {seqname: (forward_count, reverse_count)}."""
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("id") or row.get("chromosome") or list(row.values())[0]
            try:
                fwd = int(row.get("forward_repeat_number", 0) or 0)
                rev = int(row.get("reverse_repeat_number", 0) or 0)
            except ValueError:
                continue
            if fwd > 0 or rev > 0:
                result[seq] = (fwd, rev)
    return result


def run_telomere(args) -> None:
    ubbo_tel  = _load_gff3(Path(args.predicted), "telomere")
    tidk_tel  = _load_tidk(Path(args.tidk_tsv))

    ubbo_seqs = set(ubbo_tel)
    tidk_seqs = set(tidk_tel)

    both   = ubbo_seqs & tidk_seqs
    ubbo_only = ubbo_seqs - tidk_seqs
    tidk_only = tidk_seqs - ubbo_seqs
    neither_n = max(0, 0)  # unknown without reference

    rows = [
        {"genome": args.genome, "category": "both_tools",
         "n_sequences": len(both)},
        {"genome": args.genome, "category": "ubbotelorna_only",
         "n_sequences": len(ubbo_only)},
        {"genome": args.genome, "category": "tidk_only",
         "n_sequences": len(tidk_only)},
    ]

    fields = ["genome", "category", "n_sequences"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[compare_tel] {args.genome}: "
          f"both={len(both)} ubbo_only={len(ubbo_only)} tidk_only={len(tidk_only)}",
          file=sys.stderr)


# ── robustness summary mode ───────────────────────────────────────────────────

def _count_features(gff3_path: Path | str) -> dict:
    """Count rRNA and tRNA features in a GFF3."""
    counts = defaultdict(int)
    p = Path(gff3_path)
    if not p.exists() or p.stat().st_size == 0:
        return counts
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 3:
                counts[cols[2].lower()] += 1
    return dict(counts)


def run_robustness(args) -> None:
    scenarios = [
        ("native",     "ubbotelorna", args.native_ubbo),
        ("native",     "barrnap",     args.native_barrnap),
        ("adversarial","ubbotelorna", args.advers_ubbo),
        ("adversarial","barrnap",     args.advers_barrnap),
        ("fragmented", "ubbotelorna", args.fragment_ubbo),
        ("fragmented", "barrnap",     args.fragment_barrnap),
    ]

    rows = []
    for scenario, tool, path in scenarios:
        counts = _count_features(path)
        rows.append({
            "genome":   args.genome,
            "scenario": scenario,
            "tool":     tool,
            "n_rrna":   counts.get("rrna", 0),
            "n_trna":   counts.get("trna", 0),
            "n_tel":    counts.get("telomere_region", 0) + counts.get("telomere", 0),
        })

    fields = ["genome","scenario","tool","n_rrna","n_trna","n_tel"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[robustness] {args.genome}: summary written to {args.output}",
          file=sys.stderr)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compare annotations: sensitivity, precision, F1")

    ap.add_argument("--robustness_mode",    action="store_true")

    # standard / telomere mode
    ap.add_argument("--predicted")
    ap.add_argument("--reference")
    ap.add_argument("--tidk_tsv")
    ap.add_argument("--feature",   default="rRNA",
                    choices=["rRNA","tRNA","telomere"])
    ap.add_argument("--tool",      default="unknown")
    ap.add_argument("--genome",    default="unknown")
    ap.add_argument("--threshold", type=float, default=0.80)
    ap.add_argument("--output",    required=True)

    # robustness mode
    ap.add_argument("--native_ubbo")
    ap.add_argument("--native_barrnap")
    ap.add_argument("--advers_ubbo")
    ap.add_argument("--advers_barrnap")
    ap.add_argument("--fragment_ubbo")
    ap.add_argument("--fragment_barrnap")

    args = ap.parse_args()

    if args.robustness_mode:
        run_robustness(args)
    elif args.feature == "telomere":
        if not args.tidk_tsv:
            print("ERROR: --tidk_tsv required for telomere comparison", file=sys.stderr)
            sys.exit(1)
        run_telomere(args)
    else:
        if not args.predicted or not args.reference:
            print("ERROR: --predicted and --reference required", file=sys.stderr)
            sys.exit(1)
        run_standard(args)


if __name__ == "__main__":
    main()
