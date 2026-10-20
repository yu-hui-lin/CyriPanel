#!/usr/bin/env python3
"""Install the verified example data archive beside this script."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import tarfile
import tempfile

EXAMPLE = Path(__file__).resolve().parent
ROLES = dict.fromkeys(
    ["NA18522", "HG00280", "HG00133", "NA19043", "HG01530"], "reference"
)
ROLES.update(HG00639="test", HG02129="test")
ARCHIVE_NAME = "cyripanel-example-data.tar.gz"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def hashes(text):
    result = {}
    for line in text.splitlines():
        digest, name = line.split("  ", 1)
        require(re.fullmatch(r"[0-9a-f]{64}", digest) and name not in result,
                "Invalid or duplicate checksum entry")
        result[name] = digest
    return result


def transfer(source, destination=None):
    digest = hashlib.sha256()
    for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
        digest.update(block)
        if destination is not None:
            destination.write(block)
    return digest.hexdigest()


def metadata():
    data = EXAMPLE / "data"
    require(data.is_dir() and not data.is_symlink(), "Expected a real data directory")
    names = ["SHA256SUMS", "bam-validation.json", "archive-members.sha256",
             "archive-report.json", ARCHIVE_NAME + ".sha256"]
    small = {}
    for name in names:
        path = data / name
        require(path.is_file() and not path.is_symlink(), "Missing regular metadata: " + name)
        small[name] = path.read_bytes()
    report = json.loads(small["archive-report.json"])
    require(report["archive"] == ARCHIVE_NAME and report["sample_roles"] == ROLES,
            "Unexpected archive identity or sample roles")
    require(report["regular_file_count"] == 24 and report["bam_count"] == 7
            and report["index_count"] == 7 and report["header_count"] == 7
            and report["symlink_member_count"] == 0
            and report["member_checksums_passed"] is True
            and report["gzip_integrity_passed"] is True,
            "Unexpected archive verification report")
    require(hashes(small[ARCHIVE_NAME + ".sha256"].decode()) ==
            {ARCHIVE_NAME: report["archive_sha256"]}, "Archive checksum records disagree")
    inputs = hashes(small["SHA256SUMS"].decode())
    expected_inputs = {"bams/" + sid + suffix for sid in ROLES
                       for suffix in (".bam", ".bam.bai")}
    require(set(inputs) == expected_inputs, "Expected 14 BAM/index checksums")
    members = hashes(small["archive-members.sha256"].decode())
    expected_members = {"data/" + name for name in expected_inputs}
    expected_members.update("data/headers/" + sid + ".sam" for sid in ROLES)
    expected_members.update("data/" + name for name in
                            ("SHA256SUMS", "bam-validation.json", "header-preparation.json"))
    require(set(members) == expected_members, "Unexpected archive member paths")
    require(all(members["data/" + name] == digest for name, digest in inputs.items()),
            "BAM/index checksum lists disagree")
    for name in ("SHA256SUMS", "bam-validation.json"):
        require(hashlib.sha256(small[name]).hexdigest() == members["data/" + name],
                "Tracked metadata differs from archive: " + name)
    return report, members, small


def require_absent(path):
    require(not path.exists() and not path.is_symlink(),
            "Refusing to replace existing data: " + str(path))


def install(archive_path, report, members, small):
    data = EXAMPLE / "data"
    additions = ("bams", "headers", "header-preparation.json")
    for name in additions:
        require_absent(data / name)
    print("Verifying archive SHA256...", flush=True)
    with archive_path.open("rb") as raw:
        require(transfer(raw) == report["archive_sha256"], "Archive SHA256 mismatch")
        raw.seek(0)
        with tempfile.TemporaryDirectory(prefix=".install-data-", dir=EXAMPLE) as temporary:
            staging = Path(temporary)
            seen = set()
            print("Extracting and verifying 24 regular files...", flush=True)
            with gzip.GzipFile(fileobj=raw, mode="rb") as uncompressed:
                with tarfile.open(fileobj=uncompressed, mode="r|") as archive:
                    for member in archive:
                        require(member.name in members and member.name not in seen
                                and member.isfile() and not member.sparse,
                                "Unexpected, duplicate, or non-regular archive member")
                        target = staging / member.name
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with archive.extractfile(member) as source, target.open("xb") as output:
                            require(transfer(source, output) == members[member.name],
                                    "Member checksum mismatch: " + member.name)
                        require(target.stat().st_size == member.size,
                                "Extracted size mismatch: " + member.name)
                        seen.add(member.name)
                while uncompressed.read(4 * 1024 * 1024):
                    pass
            require(seen == set(members), "Archive members are missing")
            for name in ("SHA256SUMS", "bam-validation.json"):
                require((staging / "data" / name).read_bytes() == small[name],
                        "Extracted metadata differs: " + name)
            require(metadata() == (report, members, small), "Metadata changed during installation")
            for name in additions:
                require_absent(data / name)
            for name in additions:
                require_absent(data / name)
                (staging / "data" / name).rename(data / name)
    print("EXAMPLE_DATA_INSTALLED: 7 BAMs, 7 indexes, 7 headers; all 24 file hashes matched")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path, help="Downloaded archive path")
    parser.add_argument("--check", action="store_true",
                        help="Check small metadata and archive size only; no hashing or extraction")
    args = parser.parse_args()
    report, members, small = metadata()
    require(args.archive.is_file() and args.archive.stat().st_size == report["archive_bytes"],
            "Archive missing or size differs from the verification report")
    if args.check:
        print("INSTALLER_METADATA_OK: 24 members, 14 BAM/index checksums, archive size matches")
        print("No archive hashing or extraction performed; destination readiness was not checked.")
        return
    install(args.archive, report, members, small)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, EOFError, tarfile.TarError) as error:
        sys.exit("STOP: " + str(error))
