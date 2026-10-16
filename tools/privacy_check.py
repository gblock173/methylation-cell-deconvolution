"""Check allowlisted Git content and history without printing matched secrets."""
import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 200000
RULES = {
    "email": re.compile(r"[\w.+%-]+@[\w.-]+\.[A-Za-z]{2,}"),
    "ipv4": re.compile(r"(?<![\w.])(?:\d{1,3}\.){3}\d{1,3}(?![\w.])"),
    "windows_path": re.compile(r"(?i)\b[A-Z]:[\\/]"),
    "unc_path": re.compile(r"\\\\[A-Za-z0-9_.-]+\\"),
    "unix_private_path": re.compile(r"(?<![\w:])/(?!usr/bin/env\b|bin/sh\b)(?:[A-Za-z_][\w.-]+/){2,}", re.I),
    "ssh_target": re.compile(r"\b[\w.-]+@[\w.-]+:"),
    "private_key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
    "access_token": re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|AKIA[A-Z0-9]{16})\b"),
    "cohort_identifier": re.compile(r"\b(?:[A-Z]{2,6}\d{4,}|[A-Fa-f0-9]{12}|[A-Fa-f0-9]{8}-(?:[A-Fa-f0-9]{4}-){3}[A-Fa-f0-9]{12})\b"),
    "sequencing_identifier": re.compile(r"\b(?:GSM|SRR|ERR|DRR)\d{4,}\b"),
    "credential_literal": re.compile(r"""(?i)\b(?:password|token|secret|api_key)\s*[:=]\s*['"][A-Za-z0-9_/-]{12,}['"]"""),
}

def git(*args, check=True):
    result = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True)
    if check and result.returncode:
        raise RuntimeError("Cannot inspect Git state")
    return result.stdout

def allowed_names():
    return {line for line in (ROOT / ".release-files").read_text().splitlines()
            if line and not line.startswith("#")}

def scan(label, payload, failures):
    if len(payload) > MAX_BYTES:
        failures.append((label, "oversized_file"))
        return
    try:
        content = payload.decode("utf-8")
    except UnicodeDecodeError:
        failures.append((label, "binary_content"))
        return
    if "\0" in content:
        failures.append((label, "binary_content"))
    for category, pattern in RULES.items():
        if pattern.search(content):
            failures.append((label, category))

def check_entry(name, mode, payload, allowed, failures):
    if name not in allowed:
        failures.append(("unlisted file", "not_allowlisted"))
        return
    if mode not in ("100644", "100755"):
        failures.append((name, "symlink_or_submodule"))
    scan(name + " (name)", name.encode(), failures)
    scan(name, payload, failures)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--working-tree", action="store_true")
    args = parser.parse_args()
    failures, seen = [], set()
    allowed = allowed_names()
    if args.working_tree:
        for name in sorted(allowed):
            path = ROOT / name
            if not path.is_file() or path.is_symlink():
                failures.append((name, "missing_or_symlink"))
                continue
            check_entry(name, "100644", path.read_bytes(), allowed, failures)
    else:
        entries = git("ls-files", "--stage", "-z").decode().split("\0")
        for entry in filter(None, entries):
            header, name = entry.split("\t", 1)
            mode, oid, stage = header.split()
            if stage != "0":
                failures.append((name, "unmerged"))
            check_entry(name, mode, git("cat-file", "blob", oid), allowed, failures)
            seen.add(oid)
        commits = git("rev-list", "--all", check=False).decode().splitlines()
        for commit in commits:
            scan("commit metadata", git("cat-file", "commit", commit), failures)
            identities = git("show", "-s", "--format=%an%n%ae%n%cn%n%ce", commit).decode().splitlines()
            if identities != ["Yujin Jeon", "", "Yujin Jeon", ""]:
                failures.append(("commit metadata", "unexpected_identity"))
            tree = git("ls-tree", "-r", "-z", commit).decode().split("\0")
            for entry in filter(None, tree):
                header, name = entry.split("\t", 1)
                mode, kind, oid = header.split()
                if kind != "blob":
                    failures.append(("historical entry", "non_blob"))
                    continue
                if oid in seen and name in allowed and mode in ("100644", "100755"):
                    continue
                check_entry(name, mode, git("cat-file", "blob", oid), allowed, failures)
                seen.add(oid)
        for entry in git("for-each-ref", "--format=%(objecttype) %(objectname) %(refname)").decode().splitlines():
            kind, oid, name = entry.split(" ", 2)
            scan("reference name", name.encode(), failures)
            if kind == "tag":
                scan("tag metadata", git("cat-file", "tag", oid), failures)
    if failures:
        for label, category in sorted(set(failures)):
            print("FAIL:", label, "-", category)
        return 1
    print("PASS: release allowlist, content and applicable history checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
