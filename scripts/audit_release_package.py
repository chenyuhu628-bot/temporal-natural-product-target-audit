#!/usr/bin/env python3
"""Run the eight fail-closed public-release gates and write audit evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


TEXT_SUFFIXES = {".md", ".txt", ".tsv", ".csv", ".json", ".yml", ".yaml", ".cff", ".py", ".gitignore"}
FORBIDDEN_FILE_PARTS = {"raw", "interim", "processed", "scoring_inputs", "evaluation_inputs", "__pycache__", ".git"}
FORBIDDEN_NAMES = {
    "historical_pairs.tsv.gz",
    "candidate_sequences.fasta",
    "evaluation_pairs.tsv.gz",
    "corrective_prediction_ranks.tsv.gz",
    "query_compounds.tsv.gz",
}
ENTITY_PATTERNS = {
    "InChIKey": re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"),
    "NPASS identifier": re.compile(r"\b(?:NPC|NPT|NPASS)[_-]?[0-9]{3,}\b", re.I),
    "PMID value": re.compile(r"\bPMID\s*[:=_-]?\s*[0-9]{7,9}\b", re.I),
    "UniProt accession": re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b"),
}
SECRET_PATTERNS = {
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "GitHub token": re.compile(r"\b(?:ghp|github_pat)_[A-Za-z0-9_]{20,}\b"),
    "Zenodo token": re.compile(r"(?i)\bzenodo[_-]?(?:access[_-]?)?token\s*[:=]\s*['\"]?[A-Za-z0-9._-]{20,}"),
    "assigned secret": re.compile(r"(?i)\b(?:api[_-]?key|access[_-]?token|password|client[_-]?secret)\s*[:=]\s*['\"][^'\"]{8,}['\"]"),
}
ABSOLUTE_PATH = re.compile(r"(?<![A-Za-z0-9+.-])(?:[A-Za-z]:[\\/](?:Users|NPASS|home)[\\/]|/(?:home|Users)/)", re.I)
FORBIDDEN_EXACT_HEADERS = {
    "smiles",
    "inchi",
    "inchikey",
    "uniprot_accession",
    "pmid",
    "np_id",
    "query_id",
    "compound_id",
    "target_id",
    "relation_id",
    "sequence",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_excluded(path: Path, root: Path) -> bool:
    relative = path.relative_to(root)
    parts = {part.casefold() for part in relative.parts}
    return bool(parts & {".git", "__pycache__"} or path.suffix.casefold() in {".pyc", ".pyo", ".tmp"})


def files(root: Path) -> list[Path]:
    return sorted(
        (path for path in root.rglob("*") if path.is_file() and not is_excluded(path, root)),
        key=lambda path: path.relative_to(root).as_posix(),
    )


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig", errors="replace")


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def manifest_check(root: Path, manifest_relative: str, checksum_relative: str) -> dict[str, object]:
    manifest_path = root / manifest_relative
    checksum_path = root / checksum_relative
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    recorded = {row["path"]: row for row in manifest_rows}
    excluded = {manifest_relative, checksum_relative}
    actual = {
        path.relative_to(root).as_posix(): path
        for path in files(root)
        if path.relative_to(root).as_posix() not in excluded
    }
    manifest_missing = sorted(set(actual) - set(recorded))
    manifest_extra = sorted(set(recorded) - set(actual))
    manifest_bad = sorted(
        name
        for name in set(actual) & set(recorded)
        if recorded[name]["sha256"] != sha256(actual[name]) or int(recorded[name]["bytes"]) != actual[name].stat().st_size
    )

    checksums: dict[str, str] = {}
    malformed: list[int] = []
    for number, line in enumerate(read_text(checksum_path).splitlines(), 1):
        fields = line.split("  ", 1)
        if len(fields) != 2 or not re.fullmatch(r"[0-9a-f]{64}", fields[0]):
            malformed.append(number)
        else:
            checksums[fields[1]] = fields[0]
    checksum_actual = {
        path.relative_to(root).as_posix(): sha256(path)
        for path in files(root)
        if path.relative_to(root).as_posix() != checksum_relative
    }
    checksum_missing = sorted(set(checksum_actual) - set(checksums))
    checksum_extra = sorted(set(checksums) - set(checksum_actual))
    checksum_bad = sorted(name for name in set(checksum_actual) & set(checksums) if checksum_actual[name] != checksums[name])
    passed = not any([manifest_missing, manifest_extra, manifest_bad, malformed, checksum_missing, checksum_extra, checksum_bad])
    return {
        "package": root.name,
        "status": "PASS" if passed else "FAIL",
        "manifest_entries": len(recorded),
        "manifest_missing": manifest_missing,
        "manifest_extra": manifest_extra,
        "manifest_hash_or_size_mismatch": manifest_bad,
        "checksum_malformed_lines": malformed,
        "checksum_missing": checksum_missing,
        "checksum_extra": checksum_extra,
        "checksum_mismatch": checksum_bad,
    }


def code_licence_gate(repository: Path, software: Path) -> dict[str, object]:
    license_text = read_text(repository / "LICENSE")
    metadata = json.loads(read_text(repository / ".zenodo.json"))
    passed = (
        "MIT License" in license_text
        and metadata.get("license") == "MIT"
        and (software / "LICENSE").is_file()
        and sha256(repository / "LICENSE") == sha256(software / "LICENSE")
    )
    return {"status": "PASS" if passed else "FAIL", "repository_license": metadata.get("license"), "software_license_copy_exact": (software / "LICENSE").is_file() and sha256(repository / "LICENSE") == sha256(software / "LICENSE")}


def copyright_gate(repository: Path) -> dict[str, object]:
    texts = "\n".join(read_text(repository / name) for name in ["LICENSE", "NOTICE", "CITATION.cff"])
    passed = "Chenyu Hu" in texts and "2026" in texts
    return {"status": "PASS" if passed else "FAIL", "owner": "Chenyu Hu" if "Chenyu Hu" in texts else None, "year": 2026 if "2026" in texts else None}


def scan_release(repository: Path, dataset: Path, software: Path) -> tuple[dict[str, object], dict[str, object]]:
    forbidden_files: list[str] = []
    identifiers: list[dict[str, object]] = []
    absolute_paths: list[dict[str, object]] = []
    secrets: list[dict[str, object]] = []
    packages = [("repository", repository), ("dataset", dataset), ("software", software)]
    for package_name, root in packages:
        for path in files(root):
            rel = path.relative_to(root).as_posix()
            parts = {part.casefold() for part in path.relative_to(root).parts}
            if parts & FORBIDDEN_FILE_PARTS or path.name.casefold() in FORBIDDEN_NAMES or "complete_rank" in rel.casefold():
                forbidden_files.append(f"{package_name}:{rel}")
            if path.suffix.casefold() not in TEXT_SUFFIXES and path.name not in {"LICENSE", "NOTICE"}:
                continue
            text = read_text(path)
            for label, pattern in SECRET_PATTERNS.items():
                if pattern.search(text) and path.name != "audit_release_package.py":
                    secrets.append({"package": package_name, "path": rel, "finding": label})
            if path.name != "audit_release_package.py" and "ziyu" in text.casefold():
                absolute_paths.append({"package": package_name, "path": rel, "finding": "local username"})
            for number, line in enumerate(text.splitlines(), 1):
                if ABSOLUTE_PATH.search(line) and "re.compile" not in line and path.name != "audit_release_package.py":
                    absolute_paths.append({"package": package_name, "path": rel, "line": number, "finding": "absolute local path"})
            if package_name == "dataset" and path.suffix.casefold() in {".tsv", ".csv", ".md", ".json"}:
                for label, pattern in ENTITY_PATTERNS.items():
                    match = pattern.search(text)
                    if match:
                        identifiers.append({"path": rel, "finding": label, "match": match.group(0)})

    archive = software / "temporal-natural-product-target-audit-1.0.0.zip"
    archive_bad: list[str] = []
    with ZipFile(archive) as handle:
        for name in handle.namelist():
            parts = {part.casefold() for part in Path(name).parts}
            if parts & FORBIDDEN_FILE_PARTS or Path(name).name.casefold() in FORBIDDEN_NAMES or name.endswith((".pyc", ".pyo")):
                archive_bad.append(name)
    restricted_pass = not forbidden_files and not identifiers and not absolute_paths and not archive_bad
    credentials_pass = not secrets
    restricted = {
        "status": "PASS" if restricted_pass else "FAIL",
        "forbidden_files": forbidden_files,
        "identifier_findings": identifiers,
        "absolute_path_findings": absolute_paths,
        "archive_forbidden_entries": archive_bad,
        "npass_decision": "PASS_AS_LINK_ONLY_NO_REDISTRIBUTION",
    }
    credentials = {"status": "PASS" if credentials_pass else "FAIL", "findings": secrets}
    return restricted, credentials


def aggregate_gate(dataset: Path) -> dict[str, object]:
    review_path = dataset / "AGGREGATE_RELEASE_REVIEW.tsv"
    with review_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_path = {row["path"]: row for row in rows}
    data_paths = {
        path.relative_to(dataset).as_posix()
        for folder in ["figure_source_data", "table_source_data", "reproduction"]
        for path in (dataset / folder).rglob("*")
        if path.is_file()
    }
    missing_review = sorted(data_paths - set(by_path))
    extra_review = sorted(set(by_path) - data_paths)
    bad_review = sorted(
        path
        for path, row in by_path.items()
        if row["contains_identifiers"] != "no"
        or row["reverse_reconstruction_possible"] != "no"
        or row["manual_review"] != "PASS"
        or row["release_decision"] != "INCLUDE"
        or row["license"] != "CC BY 4.0"
    )
    forbidden_headers: list[dict[str, object]] = []
    for path in sorted((dataset / "figure_source_data").glob("*.tsv")) + sorted((dataset / "table_source_data").glob("*.tsv")):
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            headers = {value.casefold() for value in next(reader)}
        bad = sorted(headers & FORBIDDEN_EXACT_HEADERS)
        if bad:
            forbidden_headers.append({"path": path.relative_to(dataset).as_posix(), "headers": bad})
    passed = not missing_review and not extra_review and not bad_review and not forbidden_headers
    return {
        "status": "PASS" if passed else "FAIL",
        "reviewed_files": len(rows),
        "missing_review": missing_review,
        "extra_review": extra_review,
        "failed_review_rows": bad_review,
        "forbidden_identifier_headers": forbidden_headers,
    }


def reproduction_gate(release_root: Path) -> dict[str, object]:
    receipt_path = release_root / "audit" / "clean" / "CLEAN_ENVIRONMENT_REPRODUCTION_RECEIPT.json"
    receipt = json.loads(read_text(receipt_path))
    passed = (
        receipt.get("status") == "COMPLETED"
        and receipt.get("clean_environment_author_side_reproduction") == "COMPLETED"
        and receipt.get("source_to_frozen", {}).get("locked_inputs_exact_sha256") == 16
        and receipt.get("source_to_frozen", {}).get("required_inputs") == 16
        and receipt.get("manuscript_aggregate_validation", {}).get("checks") == 17
        and len(receipt.get("revision_and_v4_validations", [])) == 7
        and all(item.get("status") == "PASS" for item in receipt.get("revision_and_v4_validations", []))
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "author_side": receipt.get("clean_environment_author_side_reproduction"),
        "independent_third_party": receipt.get("independent_third_party_reproduction"),
        "locked_inputs": receipt.get("source_to_frozen", {}).get("locked_inputs_exact_sha256"),
        "manuscript_checks": receipt.get("manuscript_aggregate_validation", {}).get("checks"),
    }


def readme_gate(repository: Path, python: str) -> dict[str, object]:
    missing: list[str] = []
    for name in ["README.md", "README_zh-CN.md"]:
        text = read_text(repository / name)
        for match in re.finditer(r"python\s+([^\s`]+\.py)", text):
            if not (repository / match.group(1)).is_file():
                missing.append(f"{name}:{match.group(1)}")
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    commands = [
        [python, "scripts/run_reproduction.py", "--mode", "verify-chain"],
        [python, "scripts/run_reproduction.py", "--mode", "smoke"],
        [python, "-m", "unittest", "discover", "-s", "tests", "-v"],
    ]
    executions = []
    for command in commands:
        completed = subprocess.run(command, cwd=repository, env=env, capture_output=True, text=True, timeout=120)
        executions.append({"command": " ".join(command[1:]), "returncode": completed.returncode})
    required_text = all(
        value in read_text(repository / "README.md")
        for value in ["LINK_ONLY_NO_REDISTRIBUTION", "COMPLETED", "NOT YET PERFORMED", "materialize_execution_layout.py"]
    )
    passed = not missing and required_text and all(item["returncode"] == 0 for item in executions)
    return {"status": "PASS" if passed else "FAIL", "missing_command_paths": missing, "required_status_text": required_text, "executions": executions}


def included_files(repository: Path, dataset: Path, software: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for package, root, licence in [("github_repository", repository, "MIT"), ("zenodo_dataset", dataset, "CC BY 4.0"), ("zenodo_software", software, "MIT")]:
        for path in files(root):
            rows.append(
                {
                    "package": package,
                    "path": path.relative_to(root).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "licence_scope": licence,
                    "decision": "INCLUDE",
                }
            )
    return rows


def excluded_rows(repository: Path) -> list[dict[str, object]]:
    rows = [
        {"scope": "all_public_packages", "path_or_class": "NPASS rows, identifiers, structures, activities, mappings, and relations", "decision": "EXCLUDE_LINK_ONLY", "reason": "No redistribution authorization"},
        {"scope": "all_public_packages", "path_or_class": "ChEMBL, PubMed, UniProt, MMseqs2, RDKit, package binaries, and source snapshots", "decision": "EXCLUDE", "reason": "Third-party content; provide links and reconstruction instructions only"},
        {"scope": "all_public_packages", "path_or_class": "Identifier-bearing ledgers, endpoints, ranks, per-query outputs, mappings, structures, and sequences", "decision": "EXCLUDE", "reason": "Restricted source-derived row-level content"},
        {"scope": "all_public_packages", "path_or_class": "Credentials, local configuration, caches, logs, and private paths", "decision": "EXCLUDE", "reason": "Security and privacy boundary"},
        {"scope": "manuscript_history", "path_or_class": "Historical Table S6 and Table S12", "decision": "REPLACED_AND_EXCLUDED", "reason": "Stale release-state metadata replaced by sanitized public versions"},
    ]
    for path in repository.rglob("*"):
        if path.is_file() and is_excluded(path, repository):
            rows.append({"scope": "github_worktree", "path_or_class": path.relative_to(repository).as_posix(), "decision": "EXCLUDE", "reason": "Ignored local bytecode/cache artifact"})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", required=True, type=Path)
    args = parser.parse_args()
    release_root = args.release_root.resolve()
    repository = release_root / "temporal-natural-product-target-audit"
    dataset = release_root / "zenodo_dataset"
    software = release_root / "zenodo_software"
    audit = release_root / "audit"
    if not audit.is_dir():
        raise FileNotFoundError(audit)

    restricted, credentials = scan_release(repository, dataset, software)
    manifest_results = [
        manifest_check(repository, "manifests/PACKAGE_MANIFEST.tsv", "manifests/CHECKSUMS.sha256"),
        manifest_check(dataset, "MANIFEST.tsv", "CHECKSUMS.sha256"),
        manifest_check(software, "MANIFEST.tsv", "CHECKSUMS.sha256"),
    ]
    gates = {
        "code_license": code_licence_gate(repository, software),
        "code_copyright": copyright_gate(repository),
        "restricted_content_scan": restricted,
        "credentials_scan": credentials,
        "aggregate_non_reconstruction": aggregate_gate(dataset),
        "manifest_and_checksums": {"status": "PASS" if all(item["status"] == "PASS" for item in manifest_results) else "FAIL", "packages": manifest_results},
        "author_side_clean_reproduction": reproduction_gate(release_root),
        "readme_and_commands": readme_gate(repository, sys.executable),
    }
    all_pass = len(gates) == 8 and all(gate["status"] == "PASS" for gate in gates.values())
    created = datetime.now(timezone.utc).isoformat()
    report = {
        "schema_version": "jcheminform_prepublication_audit_v1",
        "created_at_utc": created,
        "overall_status": "PASS" if all_pass else "FAIL",
        "publication_mode": "public_cleared_artifacts_only",
        "gate_count": 8,
        "gates": gates,
        "npass_status": "PASS_AS_LINK_ONLY_NO_REDISTRIBUTION" if restricted["status"] == "PASS" else "FAIL",
        "independent_third_party_reproduction": "PENDING_POST_RELEASE_VALIDATION",
    }
    write_json(audit / "package_audit_report.json", report)
    write_tsv(audit / "included_files.tsv", included_files(repository, dataset, software))
    write_tsv(audit / "excluded_files.tsv", excluded_rows(repository))

    rights_source = dataset / "table_source_data" / "Table_S12_public_release_rights_and_exclusion_matrix.tsv"
    with rights_source.open("r", encoding="utf-8", newline="") as handle:
        rights = list(csv.DictReader(handle, delimiter="\t"))
    write_tsv(audit / "rights_decision_matrix.tsv", rights)

    gate_lines = "\n".join(f"- `{name}`: **{value['status']}**" for name, value in gates.items())
    markdown = f"""# Pre-publication check report

Overall status: **{'PASS' if all_pass else 'FAIL'}**  
Publication mode: `public_cleared_artifacts_only`  
Generated: {created}

## Eight fail-closed gates

{gate_lines}

NPASS: **{'PASS_AS_LINK_ONLY_NO_REDISTRIBUTION' if restricted['status'] == 'PASS' else 'FAIL'}**.  
Clean-environment author-side reproduction: **COMPLETED**.  
Independent third-party reproduction: **PENDING_POST_RELEASE_VALIDATION**.

External publication may be attempted only if this report is PASS and the
relevant authenticated GitHub/Zenodo sessions are already configured. No URL,
release, record, or DOI is inferred from local preparation.
"""
    (audit / "PRE_PUBLICATION_CHECK_REPORT.md").write_text(markdown, encoding="utf-8", newline="\n")
    (audit / "package_audit_report.md").write_text(markdown, encoding="utf-8", newline="\n")
    decision = {
        "schema_version": "publication_decision_v1",
        "created_at_utc": created,
        "decision": "LOCAL_GATES_PASS_EXTERNAL_AUTH_CHECK_ALLOWED" if all_pass else "DO_NOT_PUBLISH",
        "all_eight_gates_pass": all_pass,
        "github_publication_attempt_allowed": all_pass,
        "zenodo_publication_attempt_allowed": all_pass,
        "external_publication_completed": False,
        "reason": "All local hard gates passed; authenticated external sessions must be checked next." if all_pass else "One or more local hard gates failed.",
    }
    write_json(audit / "PUBLICATION_DECISION.json", decision)
    print(json.dumps({"status": report["overall_status"], "gates": {name: value["status"] for name, value in gates.items()}}, indent=2))
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
