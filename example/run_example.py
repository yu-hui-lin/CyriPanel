#!/usr/bin/env python3
"""Run the seven-BAM example and compare with its frozen observed outputs."""
import argparse
import csv
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

EXAMPLE = Path(__file__).resolve().parent
CODE = EXAMPLE.parent
ROLES = dict.fromkeys(
    ["NA18522", "HG00280", "HG00133", "NA19043", "HG01530"], "reference"
)
ROLES.update(HG00639="test", HG02129="test")
TESTS = [sid for sid, role in ROLES.items() if role == "test"]
REGIONS = ["CYP2D6plusREP6_hapcn2", "CYP2D7", "D7spacer_hapcn1"]


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checksum_list(root, wanted, verify):
    entries = {}
    for line in (root / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        require(re.fullmatch(r"[0-9a-f]{64}", digest) and name not in entries,
                "Invalid or duplicate checksum entry")
        entries[name] = digest
    require(set(entries) == set(wanted), "Unexpected checksum paths: " + str(root))
    for name, digest in entries.items():
        require((root / name).is_file(), "Missing file: " + name)
        if verify:
            require(sha256(root / name) == digest, "Checksum mismatch: " + name)
    return entries


def preflight():
    with (EXAMPLE / "manifest.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["sample_id", "role", "bam", "index"],
                "Unexpected manifest columns")
        rows = list(reader)
    require(len(rows) == 7 and {r["sample_id"]: r["role"] for r in rows} == ROLES,
            "Expected five references and two tests")
    identities = set()
    for row in rows:
        sid = row["sample_id"]
        for field, suffix in (("bam", ".bam"), ("index", ".bam.bai")):
            require(row[field] == "data/bams/" + sid + suffix,
                    "Unexpected manifest path: " + sid)
            # Preserve lexical paths and the explicitly declared index.
            path = EXAMPLE / row[field]
            require(path.is_file() and path.stat().st_size > 0, "Missing file: " + str(path))
            row[field] = path
        stat = row["bam"].stat()
        identities.add((stat.st_dev, stat.st_ino))
    require(len(identities) == 7, "A BAM is reused under multiple sample IDs")
    for item in ("star_caller.py", "caller", "depth_calling", "data"):
        require((CODE / item).exists(), "Missing runtime resource: " + item)
    require(sha256(CODE / "depth_calling/run_CNVPanelizer.R") ==
            "8a864517cb066fcd45af6c21ae2c9a187fc57a8cb270ebbe1a097f2a27a308b4",
            "R wrapper differs from the tested skip=0 version")
    bed = CODE / "data/PGxProbe_region_hg38.bed"
    require(list((CODE / "data").glob("*.bed")) == [bed], "Expected exactly one panel BED")
    require(sha256(bed) == "e7912fc4e9066a534189449afdea5697c1e5f2c98da5c33e29ec663581803f61",
            "Panel BED differs from the tested version")
    expected_files = ["cyripanel.tsv", "cyripanel.json"] + [
        "cnv_reports/" + sid + "_CNV_exon_level_report.csv" for sid in TESTS
    ] + ["validation/" + name for name in (
        "json-comparison.json", "panel-count-summary.tsv", "primary-comparison.tsv"
    )]
    checksum_list(EXAMPLE / "expected", expected_files, verify=True)
    data_files = ["bams/" + sid + suffix for sid in ROLES
                  for suffix in (".bam", ".bam.bai")]
    hashes = checksum_list(EXAMPLE / "data", data_files, verify=False)
    for executable in ("Rscript", "samtools"):
        require(shutil.which(executable), "Missing executable: " + executable)
    return rows, hashes


def read_table(path, delimiter, key=None):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require(reader.fieldnames, "Empty table: " + str(path))
        key = key or reader.fieldnames[0]
        rows = list(reader)
    require(rows and len({r[key] for r in rows}) == len(rows),
            "Empty table or duplicate IDs: " + str(path))
    return {r[key]: r for r in rows}


def compare(out):
    expected = EXAMPLE / "expected"
    actual = out / "results"
    calls = [read_table(p / "cyripanel.tsv", "\t", "Sample") for p in (expected, actual)]
    details = [json.loads((p / "cyripanel.json").read_text()) for p in (expected, actual)]
    require(all(set(table) == set(TESTS) for table in calls + details),
            "Unexpected result sample IDs")
    checks = []
    for sid in TESTS:
        for field in ("Genotype", "Filter"):
            a, b = calls[0][sid][field], calls[1][sid][field]
            require(all(str(data[sid].get(field)) == table[sid][field]
                        for data, table in zip(details, calls)), "TSV/JSON inconsistency")
            checks.append([sid, field, a, b, a == b])
        tables = [read_table(p / (sid + "_CNV_exon_level_report.csv"), ",")
                  for p in (expected / "cnv_reports", actual / "cnv_reports")]
        for region in REGIONS:
            a, b = [Decimal(table[region]["MeanRatio"]) for table in tables]
            checks.append([sid, region, str(a), str(b), a.is_finite() and b.is_finite() and a == b])
    with (out / "comparison.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Metric", "Expected", "Observed", "Match"])
        writer.writerows(checks)
    equal = {sid: details[0][sid] == details[1][sid] for sid in TESTS}
    (out / "json-comparison.json").write_text(json.dumps(equal, indent=2) + "\n")
    require(all(row[-1] for row in checks) and all(equal.values()),
            "Outputs differ from frozen results; inspect comparison reports")
    for sid in TESTS:
        print(sid, calls[1][sid]["Genotype"], calls[1][sid]["Filter"], flush=True)


def execute(out, rows, hashes):
    require(not out.exists() and not out.is_symlink(), "Output path already exists")
    protected = [CODE / name for name in ("caller", "depth_calling", "data")]
    protected += [EXAMPLE / "data", EXAMPLE / "expected"]
    require(not any(p.resolve() in out.resolve().parents for p in protected),
            "Output directory must be outside runtime resources and example inputs")
    out.mkdir(parents=True)
    work, results = out / "work", out / "results"
    work.mkdir()
    results.mkdir()
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME", "R_HOME", "R_LIBS", "R_LIBS_SITE"):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"
    for key in ("R_ENVIRON_USER", "R_PROFILE_USER", "R_ENVIRON", "R_PROFILE"):
        env[key] = os.devnull
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    library = out / "empty-r-user-library"
    library.mkdir()
    env["R_LIBS_USER"] = str(library)

    def run(command, name):
        with (out / (name + ".stdout.log")).open("w") as stdout, \
                (out / (name + ".stderr.log")).open("w") as stderr:
            subprocess.run(command, cwd=work, env=env, stdout=stdout, stderr=stderr, check=True)

    print("Verifying 14 BAM/index checksums...", flush=True)
    for name, digest in hashes.items():
        require(sha256(EXAMPLE / "data" / name) == digest, "Checksum mismatch: " + name)
    (out / "input-SHA256SUMS").write_bytes((EXAMPLE / "data/SHA256SUMS").read_bytes())
    shutil.copy2(CODE / "star_caller.py", work / "star_caller.py")
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc", "*.bam", "*.bai", "*.csi",
                                    "*.cram", "cnv_panelizer_results", "results")
    for name in ("caller", "depth_calling", "data"):
        shutil.copytree(CODE / name, work / name, ignore=ignore)
    sources = {p.relative_to(work).as_posix(): sha256(p)
               for p in work.rglob("*") if p.is_file()}
    (out / "source-sha256.json").write_text(json.dumps(sources, indent=2, sort_keys=True) + "\n")
    for name in ("ref_dir", "input_tests"):
        (work / name).mkdir()
    tests = []
    for row in rows:
        folder = "ref_dir" if row["role"] == "reference" else "input_tests"
        bam = work / folder / (row["sample_id"] + ".bam")
        bam.symlink_to(row["bam"])
        Path(str(bam) + ".bai").symlink_to(row["index"])
        if row["role"] == "test":
            tests.append(bam)
    manifest = work / "test.manifest"
    manifest.write_text("".join(str(p) + "\n" for p in tests))
    run(["samtools", "quickcheck", "-v"] + [str(r["bam"]) for r in rows], "bam-check")
    run(["samtools", "--version"], "samtools-version")
    run([sys.executable, "-I", "-c",
         "import sys,numpy,scipy,pysam,pandas,statsmodels; print(sys.version); "
         "[print(p.__name__,p.__version__) for p in (numpy,scipy,pysam,pandas,statsmodels)]"],
        "python-versions")
    run(["Rscript", "--vanilla", "-e",
         'suppressPackageStartupMessages(library(CNVPanelizer));'
         'stopifnot(as.character(packageVersion("CNVPanelizer"))=="1.30.0");'
         'print(.libPaths()); print(find.package("CNVPanelizer")); sessionInfo()'], "R-sessionInfo")
    print("Running CyriPanel for HG00639 and HG02129...", flush=True)
    run([sys.executable, "-u", str(work / "star_caller.py"), "--manifest", str(manifest),
         "--genome", "38", "--outDir", str(results), "--prefix", "cyripanel", "--threads", "1"],
        "caller")
    logs = "\n".join((out / ("caller." + stream + ".log")).read_text(errors="replace")
                     for stream in ("stdout", "stderr"))
    require(not re.search(r"pysam fallback|Fallback CN result", logs, re.I),
            "CNVPanelizer depth fallback detected")
    (results / "cnv_reports").mkdir()
    for sid in TESTS:
        name = sid + "_CNV_exon_level_report.csv"
        shutil.copy2(work / "depth_calling/cnv_panelizer_results" / name, results / "cnv_reports" / name)
    compare(out)
    (out / "EXAMPLE_COMPLETE").write_text("Frozen observed outputs matched; not a biological truth assessment.\n")
    print("EXAMPLE_OUTPUTS_MATCH: 10 metrics and 2 complete parsed JSON records", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Lightweight preflight; no BAM hashing or analysis")
    parser.add_argument("--out", type=Path, help="New output directory")
    args = parser.parse_args()
    if not args.check and args.out is None:
        parser.error("--out is required unless --check is used")
    rows, hashes = preflight()
    if args.check:
        print("EXAMPLE_PREFLIGHT_OK: 5 references, 2 tests, 7 frozen outputs verified")
        print("BAM checksums and analysis will run in the compute job.")
        return
    execute(Path(os.path.abspath(args.out.expanduser())), rows, hashes)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("FAILED: " + str(exc), file=sys.stderr)
        sys.exit(1)
