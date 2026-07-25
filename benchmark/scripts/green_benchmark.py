#!/usr/bin/env python3
"""Fast (6-genome) codecarbon benchmark: UbboTELORNA vs. barrnap +
tRNAscan-SE + tidk, measuring real emissions per tool rather than just
wall-clock time.

Every tool invocation is wrapped in its own fresh codecarbon
EmissionsTracker, so all four tools are measured with the same
methodology (UbboTELORNA's own built-in tracker is disabled here via
--disable_co2_tracking, to avoid mixing self-reported and externally
observed numbers).

UbboTELORNA's CORE pipeline (Modules 0, 2, 3, 4, 5 in the current
numbering: telomere, masking, rRNA, tRNA, integration) is run once per
genome and compared against the SUM of barrnap + tRNAscan-SE + tidk run
separately -- the same "full pipeline" framing used in the main
benchmark's README. Module 1 (Subtelomeric tandem repeats, added
v0.4.0) and Module 6/7 (Visualization/Evolutionary analysis) are
excluded from this head-to-head: Module 1 has no equivalent on the
comparator side (no standalone tool does subtelomeric completeness
tiering), so including it would shrink UbboTELORNA's apparent
speed/emissions advantage without a fair addition to the other side of
the comparison -- not a like-for-like measurement. Instead, Module 1's
own marginal cost is measured SEPARATELY (core pipeline with vs.
without Module 1 included) and reported in its own table, since that's
the honest way to report a capability with no comparator: state its
real cost plainly, don't fold it into a comparison it was never part of.

Assumes the 6 genomes below are already downloaded (see `download_all` in
the main Snakefile, or run `download_genomes.py` directly for just these).

Usage:
    conda activate ubbotelorna_bench   # barrnap, tRNAscan-SE, tidk, snakemake, ...
    python3 green_benchmark.py --outdir results/green/ --threads 8

    # UbboTELORNA itself runs in its own env, activated automatically
    # if `ubbotelorna_env` in config.yaml names a conda env with `conda run`;
    # otherwise it must be importable/runnable from the current env too.
"""

import argparse
import csv
import shutil
import subprocess
import sys
import time
import yaml
from pathlib import Path

_basedir = Path(__file__).resolve().parent.parent  # benchmark/

# Six genomes spanning ~3 orders of magnitude in size, all already defined
# in config.yaml (subset of performance_genomes, plus csinensis) so results
# are directly comparable to the existing 40-genome correctness benchmark.
# csinensis (Citrus sinensis, GCF_022201045.2) added specifically to
# exercise Module 1 (Subtelomeric tandem repeats): a real, publicly citable
# draft-quality plant assembly with `telomere_repeat: "auto"`, the exact
# scenario Module 1's completeness tiering is designed for.
GENOMES = ["ecoli_k12", "scerevisiae", "athaliana", "celegans", "osativa",
          "csinensis"]


def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(f"ERROR: '{name}' not found in PATH. Activate the benchmark "
              f"environment first: conda activate ubbotelorna_bench",
              file=sys.stderr)
        sys.exit(1)
    return tool


def _load_config() -> dict:
    with open(_basedir / "config.yaml") as fh:
        return yaml.safe_load(fh)


def _genome_fasta(config: dict, genome: str) -> Path:
    fasta = Path(config["results_dir"]) / "genomes" / genome / f"{genome}.fasta"
    if not fasta.exists():
        print(f"ERROR: {fasta} not found. Download it first, e.g.:\n"
              f"  python3 scripts/download_genomes.py --accession "
              f"{config['genomes'][genome]['accession']} --outdir "
              f"{fasta.parent} --genome {genome}", file=sys.stderr)
        sys.exit(1)
    return fasta


def run_tracked(cmd: list, label: str, genome: str, outdir: Path,
                run_fn=subprocess.run) -> dict:
    """Run `cmd` (or call `run_fn` with no args, for multi-step tools like
    tidk) inside a fresh codecarbon EmissionsTracker. Returns a result dict;
    never raises on tool failure — failures are recorded, not fatal, so one
    bad genome/tool doesn't abort the whole (fast, small) benchmark."""
    from codecarbon import EmissionsTracker

    tracker = EmissionsTracker(
        project_name=f"{label}_{genome}",
        output_dir=str(outdir), output_file="emissions.csv",
        log_level="warning", save_to_file=True,
    )
    tracker.start()
    t0 = time.monotonic()
    ok, err = True, ""
    try:
        if cmd is not None:
            result = subprocess.run(cmd, capture_output=True, text=True)
            ok = result.returncode == 0
            err = result.stderr[-2000:] if not ok else ""
        else:
            run_fn()
    except Exception as e:
        ok, err = False, str(e)
    elapsed_s = time.monotonic() - t0
    emissions_kg = tracker.stop()

    # Flatten embedded newlines/tabs in captured stderr so every row stays
    # on exactly one line -- otherwise plain-text tools (cut/awk/grep)
    # can't parse the file even when the delimiter itself is correct.
    err_flat = err.replace("\n", " | ").replace("\t", " ") if err else err
    row = {"genome": genome, "tool": label, "ok": ok,
          "elapsed_s": round(elapsed_s, 2),
          "emissions_kg_co2eq": emissions_kg, "error": err_flat}
    status = "OK" if ok else "FAILED"
    print(f"[{status}] {label:12s} {genome:14s} "
          f"{elapsed_s:7.1f}s  {emissions_kg or 0:.6f} kg CO2eq")
    return row


def _ubbotelorna_cmd(config: dict, genome: str, fasta: Path, threads: int,
                     run_dir: Path, skip_module: str) -> list:
    """UbboTELORNA's own dependencies (tantan, infernal/cmsearch, aragorn,
    nhmmer, trf) live in its own conda env (`ubbotelorna`), separate from
    this benchmark's `ubbotelorna_bench` env -- mirrors the main Snakefile's
    `conda: "../envs/UbboTELORNA.yaml"` directive on its UbboTELORNA rules.
    Running it with whatever env happens to be active otherwise fails
    partway through with a missing dependency."""
    env_name = config.get("ubbotelorna_conda_env", "ubbotelorna")
    return [
        "conda", "run", "-n", env_name, "--no-capture-output",
        "python3", str((_basedir / config["ubbotelorna_script"]).resolve()),
        "--fasta", str(fasta), "--output", str(run_dir),
        "--kingdom", config["genomes"][genome]["kingdom"],
        "--threads", str(threads),
        "--disable_co2_tracking",   # avoid double-tracking; see module docstring
        "--skip_module", skip_module,
        "--force",
    ]


def run_ubbotelorna(config: dict, genome: str, fasta: Path, threads: int,
                    outdir: Path) -> dict:
    """CORE pipeline only: Modules 0, 2, 3, 4, 5 (telomere, masking, rRNA,
    tRNA, integration) -- the head-to-head comparison against
    barrnap+tRNAscan-SE+tidk. Module 1 (subtelomeric, no comparator) and
    Modules 6/7 (visualization/evolutionary analysis, not part of the core
    annotation task being compared) are skipped here; see
    run_ubbotelorna_with_subtelomeric() for Module 1's own cost, measured
    separately rather than folded into this comparison."""
    run_dir = outdir / "runs" / f"ubbotelorna_{genome}"
    cmd = _ubbotelorna_cmd(config, genome, fasta, threads, run_dir,
                           skip_module="1,6,7")
    return run_tracked(cmd, "ubbotelorna", genome, outdir / "emissions")


def run_ubbotelorna_with_subtelomeric(config: dict, genome: str, fasta: Path,
                                      threads: int, outdir: Path) -> dict:
    """Same CORE pipeline as run_ubbotelorna(), plus Module 1 (subtelomeric
    tandem repeats) included. Reported in its own table, not mixed into the
    vs.-comparator numbers above -- there is no equivalent capability on
    the barrnap+tRNAscan-SE+tidk side to compare it against, so the honest
    way to report it is its own real marginal cost, not an artificially
    shrunk "advantage" over tools that don't do this at all."""
    run_dir = outdir / "runs" / f"ubbotelorna_subtel_{genome}"
    cmd = _ubbotelorna_cmd(config, genome, fasta, threads, run_dir,
                           skip_module="6,7")
    return run_tracked(cmd, "ubbotelorna_with_subtelomeric", genome,
                       outdir / "emissions")


def run_barrnap(config: dict, genome: str, fasta: Path, threads: int,
                outdir: Path) -> dict:
    kingdom = config["genomes"][genome]["barrnap_kingdom"]
    gff3_out = outdir / "runs" / f"barrnap_{genome}.gff3"
    gff3_out.parent.mkdir(parents=True, exist_ok=True)

    def _run():
        with open(gff3_out, "w") as fh:
            subprocess.run(["barrnap", "--kingdom", kingdom, "--threads",
                           str(threads), str(fasta)], stdout=fh,
                          stderr=subprocess.PIPE)

    return run_tracked(None, "barrnap", genome, outdir / "emissions", run_fn=_run)


def run_trnascan(config: dict, genome: str, fasta: Path, threads: int,
                 outdir: Path) -> dict:
    """tRNAscan-SE prompts interactively to overwrite --gff/-m if they
    already exist and hangs forever waiting for stdin that never comes on
    a reruns -- same issue already fixed for the perf sweep in the main
    Snakefile (see git history). Delete stale outputs first so a rerun of
    this script is always safe."""
    mode = config["genomes"][genome]["trnascan_mode"]
    run_dir = outdir / "runs" / f"trnascan_{genome}"
    run_dir.mkdir(parents=True, exist_ok=True)
    gff3_out = run_dir / f"{genome}.gff3"
    stats_out = run_dir / f"{genome}_stats.txt"
    gff3_out.unlink(missing_ok=True)
    stats_out.unlink(missing_ok=True)
    cmd = ["tRNAscan-SE", mode, "--thread", str(threads),
          "--gff", str(gff3_out), "-m", str(stats_out), str(fasta)]
    return run_tracked(cmd, "trnascan", genome, outdir / "emissions")


def run_tidk(config: dict, genome: str, fasta: Path, outdir: Path) -> dict:
    """Mirrors the main Snakefile's run_tidk logic exactly (explore-then-
    search when the repeat unit isn't already known), wrapped as one
    tracked unit since that's the equivalent single task UbboTELORNA's
    Module 0 performs internally in one step."""
    repeat = config["genomes"][genome]["telomere_repeat"]
    run_dir = outdir / "runs" / f"tidk_{genome}"
    run_dir.mkdir(parents=True, exist_ok=True)
    # Precautionary: same reruns-hang risk as tRNAscan-SE if tidk also
    # prompts on an existing output (unconfirmed, but cheap to guard against).
    (run_dir / f"{genome}_telomeric_repeat_windows.tsv").unlink(missing_ok=True)

    def _run():
        rep = repeat
        if rep == "auto":
            r = subprocess.run(["tidk", "explore", "--minimum", "5",
                               "--maximum", "12", str(fasta)],
                              capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"tidk explore failed:\n{r.stderr}")
            rep = None
            for line in r.stdout.splitlines():
                cols = line.strip().split("\t")
                if cols and cols[0] and all(c in "ACGTacgt" for c in cols[0]):
                    rep = cols[0].upper()
                    break
            if rep is None:
                raise RuntimeError(f"tidk explore found no telomeric repeat "
                                  f"in {fasta}")
        r2 = subprocess.run(["tidk", "search", "--string", rep,
                            "--output", genome, "--dir", str(run_dir),
                            str(fasta)], capture_output=True, text=True)
        if r2.returncode != 0:
            raise RuntimeError(f"tidk search failed:\n{r2.stderr}")

    return run_tracked(None, "tidk", genome, outdir / "emissions", run_fn=_run)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=Path("results/green"),
                    help="Output directory (default: results/green)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--genomes", nargs="+", default=GENOMES,
                    help=f"Genome keys to benchmark (default: {' '.join(GENOMES)})")
    args = ap.parse_args(argv)

    for tool in ("barrnap", "tRNAscan-SE", "tidk"):
        _require_tool(tool)
    try:
        import codecarbon  # noqa: F401
    except ImportError:
        print("ERROR: codecarbon not installed. "
              "conda install -c conda-forge codecarbon", file=sys.stderr)
        sys.exit(1)

    config = _load_config()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "emissions").mkdir(parents=True, exist_ok=True)

    rows = []
    for genome in args.genomes:
        fasta = _genome_fasta(config, genome)
        print(f"\n=== {genome} ({config['genomes'][genome]['organism']}, "
              f"{config['genomes'][genome]['size_mb']} Mb) ===")
        rows.append(run_ubbotelorna(config, genome, fasta, args.threads, args.outdir))
        rows.append(run_ubbotelorna_with_subtelomeric(config, genome, fasta,
                                                       args.threads, args.outdir))
        rows.append(run_barrnap(config, genome, fasta, args.threads, args.outdir))
        rows.append(run_trnascan(config, genome, fasta, args.threads, args.outdir))

        # Genomes with telomere_repeat: null (e.g. circular bacterial
        # chromosomes like ecoli_k12) have no linear telomere to search for
        # at all -- mirrors the main Snakefile's own exclusion
        # (`if config["genomes"][g].get("telomere_repeat")`), which this
        # script missed on first pass and crashed on (tidk given a literal
        # None as --string).
        if config["genomes"][genome].get("telomere_repeat"):
            rows.append(run_tidk(config, genome, fasta, args.outdir))
        else:
            print(f"[SKIP]   tidk         {genome:14s} "
                  f"no telomere_repeat configured (circular genome)")
            rows.append({"genome": genome, "tool": "tidk", "ok": None,
                        "elapsed_s": None, "emissions_kg_co2eq": None,
                        "error": "skipped: telomere_repeat is null in config.yaml"})

    summary_path = args.outdir / "green_benchmark_summary.tsv"
    with open(summary_path, "w", newline="") as fh:
        # Actually tab-delimited (previously defaulted to comma despite the
        # .tsv name), and the "error" field is already newline-flattened in
        # run_tracked() -- both needed for cut/awk/grep to work on this file,
        # since the error column can otherwise contain embedded newlines
        # and commas that break naive line-based parsing regardless of
        # delimiter.
        writer = csv.DictWriter(fh, fieldnames=["genome", "tool", "ok",
                                               "elapsed_s", "emissions_kg_co2eq",
                                               "error"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWritten: {summary_path}")

    # Quick "core pipeline" comparison: UbboTELORNA (Modules 0,2,3,4,5) vs.
    # barrnap+trnascan+tidk summed. Module 1 (subtelomeric) deliberately
    # excluded from both sides -- see module docstring.
    comparator_tools = {"barrnap", "trnascan", "tidk"}
    print("\n=== Core-pipeline comparison (per genome) ===")
    print(f"{'Genome':14s} {'UbboTELORNA':>14s} {'barrnap+tRNAscan+tidk':>24s} {'Delta':>10s}")
    for genome in args.genomes:
        g_rows = [r for r in rows if r["genome"] == genome and r["ok"]]
        ubbo = sum(r["emissions_kg_co2eq"] or 0 for r in g_rows if r["tool"] == "ubbotelorna")
        others = sum(r["emissions_kg_co2eq"] or 0 for r in g_rows if r["tool"] in comparator_tools)
        delta = ubbo - others
        print(f"{genome:14s} {ubbo:14.6f} {others:24.6f} {delta:+10.6f}")

    # Module 1 (Subtelomeric tandem repeats) marginal cost, reported
    # separately since there's no comparator to measure it against: real
    # cost of the new completeness-tiering capability, stated plainly
    # rather than folded into (and distorting) the comparison above.
    print("\n=== Module 1 (Subtelomeric tandem repeats) marginal cost (per genome) ===")
    print(f"{'Genome':14s} {'Core-only (kg)':>16s} {'Core+Module1 (kg)':>20s} {'Module1 cost (kg)':>20s}")
    for genome in args.genomes:
        g_rows = [r for r in rows if r["genome"] == genome and r["ok"]]
        core = sum(r["emissions_kg_co2eq"] or 0 for r in g_rows if r["tool"] == "ubbotelorna")
        with_subtel = sum(r["emissions_kg_co2eq"] or 0 for r in g_rows
                          if r["tool"] == "ubbotelorna_with_subtelomeric")
        module1_cost = with_subtel - core
        print(f"{genome:14s} {core:16.6f} {with_subtel:20.6f} {module1_cost:20.6f}")


if __name__ == "__main__":
    main()
