"""Compare manuscript-reported headline values with local aggregate source TSVs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def document_text(path: Path) -> str:
    with ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))
    return "\n".join(
        "".join((node.text or "") for node in paragraph.iter() if node.tag in (W + "t", W + "delText"))
        for paragraph in root.iter(W + "p")
    )


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric(rows: list[dict[str, str]], scope: str, baseline: str, scenario: str) -> float:
    matches = [
        row for row in rows
        if row["provenance_scope"] == scope
        and row["baseline"] == baseline
        and row["scenario_or_subset"] == scenario
        and row["metric"] == "Recall@50"
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one metric row, found {len(matches)} for {scope}/{baseline}/{scenario}")
    return float(matches[0]["estimate"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--manuscript", required=True, type=Path)
    args = parser.parse_args()
    manuscript = document_text(args.manuscript)
    figure_root = args.project_root / "manuscript_molecular_diversity_v4_20260729/figures"
    table_root = args.project_root / "manuscript_molecular_diversity_v4_20260729/tables"
    figure_1 = tsv(figure_root / "Figure_1_source_data.tsv")
    figure_2 = tsv(figure_root / "Figure_2_source_data.tsv")
    rights = tsv(table_root / "Table_S12_rights_and_controlled_access.tsv")
    values = {
        "morgan_temporal_primary": metric(figure_2, "temporal_strict_ab", "weighted_morgan_transfer", "day_only_conservative"),
        "morgan_temporal_interval": metric(figure_2, "temporal_strict_ab", "weighted_morgan_transfer", "interval_certain_pre_cutoff"),
        "pair_temporal_primary": metric(figure_2, "temporal_strict_ab", "structure_sequence_pair_neighbor", "day_only_conservative"),
        "pair_temporal_interval": metric(figure_2, "temporal_strict_ab", "structure_sequence_pair_neighbor", "interval_certain_pre_cutoff"),
        "morgan_scaffold_primary": metric(figure_2, "scaffold_cold_strict_ab", "weighted_morgan_transfer", "day_only_conservative"),
        "morgan_scaffold_interval": metric(figure_2, "scaffold_cold_strict_ab", "weighted_morgan_transfer", "interval_certain_pre_cutoff"),
        "pair_scaffold_primary": metric(figure_2, "scaffold_cold_strict_ab", "structure_sequence_pair_neighbor", "day_only_conservative"),
        "pair_scaffold_interval": metric(figure_2, "scaffold_cold_strict_ab", "structure_sequence_pair_neighbor", "interval_certain_pre_cutoff"),
    }
    expected_rounded = {
        "morgan_temporal_primary": "0.2442",
        "morgan_temporal_interval": "0.2529",
        "pair_temporal_primary": "0.2480",
        "pair_temporal_interval": "0.2510",
        "morgan_scaffold_primary": "0.2898",
        "morgan_scaffold_interval": "0.3106",
        "pair_scaffold_primary": "0.2784",
        "pair_scaffold_interval": "0.2879",
    }
    checks: list[dict[str, object]] = []
    for name, expected in expected_rounded.items():
        observed = f"{values[name]:.4f}"
        checks.append({"check": name, "expected": expected, "observed": observed, "manuscript_contains": expected in manuscript, "status": "PASS" if observed == expected and expected in manuscript else "FAIL"})
    headline = {row["item"]: row["value"] for row in figure_1 if row["panel"] == "A"}
    for item, expected in {
        "historical_strict_v2_source_rows": "20647",
        "day_only_selected_rows": "13885",
        "interval_certain_selected_rows": "20455",
        "tier_B_to_A_upgrades": "141",
    }.items():
        observed = headline.get(item)
        manuscript_form = f"{int(expected):,}"
        checks.append({"check": item, "expected": expected, "observed": observed, "manuscript_contains": manuscript_form in manuscript, "status": "PASS" if observed == expected and manuscript_form in manuscript else "FAIL"})
    checks.extend(
        [
            {"check": "methods_present", "status": "PASS" if "2. Methods" in manuscript else "FAIL"},
            {"check": "data_availability_present", "status": "PASS" if "Data availability" in manuscript else "FAIL"},
            {"check": "code_availability_present", "status": "PASS" if "Code availability" in manuscript else "FAIL"},
            {"check": "supplementary_table_s12_legend_present", "status": "PASS" if ("Supplementary Table S12 | Rights and controlled-access matrix" in manuscript or "Supplementary Table S12 | Public-release rights and exclusion matrix" in manuscript) else "FAIL"},
            {"check": "rights_matrix_has_ten_families", "observed": len(rights), "status": "PASS" if len(rights) == 10 else "FAIL"},
        ]
    )
    result = {"status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL", "checks": checks, "note": "Consistency check only. Cleared non-reconstructive aggregates are authorized; NPASS and row-level source-derived data remain excluded."}
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
