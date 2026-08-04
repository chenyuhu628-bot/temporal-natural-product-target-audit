"""Aggregate-only PMID/source-document concentration and component audit."""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from audit_common import (
    PROTOCOL_ID,
    choose_field,
    distribution,
    finalize_manifest,
    gini_nonnegative,
    input_descriptor,
    open_dict_reader,
    require_new_output_dir,
    require_protocol_lock,
    write_json_new,
    write_tsv_new,
)


AUDIT_ID = "pmid_source_concentration_v1"
LABEL_PATTERN = re.compile(r"^[a-z][a-z0-9_]{0,47}$")
LOCKED_ROW_STATUSES = {
    "eligible_pre_cutoff",
    "excluded_not_numeric_pmid",
    "excluded_pubmed_not_found",
    "excluded_non_day_precision",
    "excluded_after_cutoff",
}


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.size: dict[str, int] = {}

    def add(self, item: str) -> None:
        if item not in self.parent:
            self.parent[item] = item
            self.size[item] = 1

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.size[root_left] < self.size[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        self.size[root_left] += self.size[root_right]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument(
        "--cohort",
        required=True,
        action="append",
        help="Repeat LABEL=EVIDENCE_ROWS.tsv.gz; rows need pair, compound, ref_id_type, ref_id",
    )
    result.add_argument(
        "--pair-filter",
        action="append",
        default=[],
        help="Optional repeat LABEL=PAIR_LEDGER.tsv.gz defining the exact relation/query denominator",
    )
    result.add_argument(
        "--row-status-ledger",
        action="append",
        default=[],
        help="Optional repeat LABEL=ROW_STATUS_LEDGER.tsv.gz for date/precision attrition",
    )
    result.add_argument("--document-type", default="PMID")
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def labeled_paths(values: Iterable[str], role: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{role} must use LABEL=PATH syntax")
        label, raw_path = value.split("=", 1)
        if not LABEL_PATTERN.fullmatch(label):
            raise ValueError(f"Unsafe or invalid cohort label: {label!r}")
        if label in result:
            raise ValueError(f"Duplicate {role} label: {label}")
        result[label] = Path(raw_path)
    return result


def load_pair_filter(path: Path) -> tuple[set[str], dict[str, str]]:
    pairs: set[str] = set()
    pair_to_query: dict[str, str] = {}
    with open_dict_reader(path) as reader:
        fields = reader.fieldnames or []
        pair_field = choose_field(fields, ["canonical_pair_key", "pair_key"], "pair filter")
        query_field = choose_field(
            fields,
            ["inchikey_full", "query_compound_inchikey_full"],
            "pair filter",
        )
        for row in reader:
            pair = row[pair_field].strip()
            query = row[query_field].strip()
            if not pair or not query:
                raise ValueError("Pair filter contains an empty relation or query identifier")
            if pair in pair_to_query and pair_to_query[pair] != query:
                raise ValueError("A pair filter relation maps to multiple query compounds")
            pairs.add(pair)
            pair_to_query[pair] = query
    if not pairs:
        raise ValueError("Pair filter is empty")
    return pairs, pair_to_query


def component_summary(
    left_values: set[str], right_values: set[str], edges: set[tuple[str, str]], left_prefix: str
) -> dict[str, Any]:
    union = UnionFind()
    left_nodes = {value: f"{left_prefix}:{value}" for value in left_values}
    right_nodes = {value: f"document:{value}" for value in right_values}
    for node in left_nodes.values():
        union.add(node)
    for node in right_nodes.values():
        union.add(node)
    for left, right in edges:
        union.union(left_nodes[left], right_nodes[right])

    left_counts: Counter[str] = Counter()
    right_counts: Counter[str] = Counter()
    for node in left_nodes.values():
        left_counts[union.find(node)] += 1
    for node in right_nodes.values():
        right_counts[union.find(node)] += 1
    roots = set(left_counts) | set(right_counts)
    largest_left = max(left_counts.values(), default=0)
    largest_right = max(right_counts.values(), default=0)
    largest_total = max(
        (left_counts[root] + right_counts[root] for root in roots), default=0
    )
    return {
        "component_count": len(roots),
        "isolated_left_node_count": sum(
            left_counts[root] == 1 and right_counts[root] == 0 for root in roots
        ),
        "largest_component_left_node_count": largest_left,
        "largest_component_left_node_fraction": (
            largest_left / len(left_values) if left_values else None
        ),
        "largest_component_source_document_count": largest_right,
        "largest_component_source_document_fraction": (
            largest_right / len(right_values) if right_values else None
        ),
        "largest_component_total_node_count": largest_total,
    }


def concentration_summary(
    *,
    label: str,
    evidence_path: Path,
    document_type: str,
    allowed_pairs: set[str] | None,
    filter_pair_to_query: dict[str, str] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    pair_to_query = dict(filter_pair_to_query or {})
    all_pairs = set(allowed_pairs or [])
    pair_documents: dict[str, set[str]] = defaultdict(set)
    evidence_row_count = 0
    selected_document_row_count = 0
    filtered_out_row_count = 0
    nonselected_reference_row_count = 0

    with open_dict_reader(evidence_path) as reader:
        fields = reader.fieldnames or []
        pair_field = choose_field(fields, ["canonical_pair_key", "pair_key"], "evidence")
        query_field = choose_field(
            fields,
            ["inchikey_full", "query_compound_inchikey_full"],
            "evidence",
        )
        reference_type_field = choose_field(fields, ["ref_id_type", "reference_type"], "evidence")
        reference_field = choose_field(fields, ["ref_id", "reference_id"], "evidence")
        for row in reader:
            evidence_row_count += 1
            pair = row[pair_field].strip()
            query = row[query_field].strip()
            if not pair or not query:
                raise ValueError("Evidence contains an empty relation or query identifier")
            if allowed_pairs is not None and pair not in allowed_pairs:
                filtered_out_row_count += 1
                continue
            all_pairs.add(pair)
            if pair in pair_to_query and pair_to_query[pair] != query:
                raise ValueError("An evidence relation maps to multiple query compounds")
            pair_to_query[pair] = query
            if row[reference_type_field].strip().upper() != document_type.upper():
                nonselected_reference_row_count += 1
                continue
            document = row[reference_field].strip()
            if not document:
                raise ValueError("Selected source-document row has an empty document identifier")
            if document_type.upper() == "PMID" and not document.isdigit():
                raise ValueError("A PMID source-document value is not numeric")
            pair_documents[pair].add(document)
            selected_document_row_count += 1

    if allowed_pairs is not None and set(pair_to_query) != allowed_pairs:
        missing = len(allowed_pairs.difference(pair_to_query))
        raise ValueError(f"Pair filter has {missing} relations without a query mapping")
    if not all_pairs:
        raise ValueError(f"Cohort {label} contains no relations")
    all_queries = set(pair_to_query.values())
    documents = set().union(*(pair_documents.values())) if pair_documents else set()
    pair_document_edges = {
        (pair, document) for pair, docs in pair_documents.items() for document in docs
    }
    query_documents: dict[str, set[str]] = defaultdict(set)
    for pair, docs in pair_documents.items():
        query_documents[pair_to_query[pair]].update(docs)
    query_document_edges = {
        (query, document) for query, docs in query_documents.items() for document in docs
    }
    pair_degree = [len(pair_documents.get(pair, set())) for pair in all_pairs]
    query_degree = [len(query_documents.get(query, set())) for query in all_queries]
    document_pair_degree = Counter(
        document for _, document in pair_document_edges
    )
    document_query_degree = Counter(
        document for _, document in query_document_edges
    )
    edge_count = len(pair_document_edges)
    sorted_document_edges = sorted(document_pair_degree.values(), reverse=True)
    hhi = (
        sum((degree / edge_count) ** 2 for degree in sorted_document_edges)
        if edge_count
        else None
    )
    top_one_percent_count = max(1, math.ceil(0.01 * len(sorted_document_edges))) if sorted_document_edges else 0

    summary = {
        "cohort": label,
        "source_document_type": document_type.upper(),
        "evidence_row_count": evidence_row_count,
        "selected_source_document_row_count": selected_document_row_count,
        "filtered_out_evidence_row_count": filtered_out_row_count,
        "nonselected_reference_row_count": nonselected_reference_row_count,
        "relation_count": len(all_pairs),
        "query_count": len(all_queries),
        "unique_source_document_count": len(documents),
        "relation_source_document_edge_count": edge_count,
        "relation_count_with_source_document": sum(value > 0 for value in pair_degree),
        "relation_count_without_source_document": sum(value == 0 for value in pair_degree),
        "query_count_with_source_document": sum(value > 0 for value in query_degree),
        "query_count_without_source_document": sum(value == 0 for value in query_degree),
        "source_documents_per_relation_distribution": distribution(pair_degree),
        "source_documents_per_query_distribution": distribution(query_degree),
        "relations_per_source_document_distribution": distribution(
            list(document_pair_degree.values())
        ),
        "queries_per_source_document_distribution": distribution(
            list(document_query_degree.values())
        ),
        "largest_source_document_relation_count": max(sorted_document_edges, default=0),
        "largest_source_document_relation_fraction": (
            max(sorted_document_edges, default=0) / len(all_pairs) if all_pairs else None
        ),
        "largest_source_document_edge_share": (
            max(sorted_document_edges, default=0) / edge_count if edge_count else None
        ),
        "top5_source_document_edge_share": (
            sum(sorted_document_edges[:5]) / edge_count if edge_count else None
        ),
        "top10_source_document_edge_share": (
            sum(sorted_document_edges[:10]) / edge_count if edge_count else None
        ),
        "top1pct_source_document_edge_share": (
            sum(sorted_document_edges[:top_one_percent_count]) / edge_count if edge_count else None
        ),
        "source_document_edge_hhi": hhi,
        "source_document_effective_number": 1.0 / hhi if hhi else None,
        "source_document_relation_degree_gini": gini_nonnegative(sorted_document_edges),
        "relation_source_component_summary": component_summary(
            all_pairs, documents, pair_document_edges, "relation"
        ),
        "query_source_component_summary": component_summary(
            all_queries, documents, query_document_edges, "query"
        ),
    }
    internal = {
        "pairs": all_pairs,
        "documents": documents,
        "pair_documents": pair_documents,
    }
    return summary, internal


def audit_row_status(label: str, path: Path) -> list[dict[str, Any]]:
    counts: Counter[tuple[str, str, str]] = Counter()
    pair_sets: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    with open_dict_reader(path) as reader:
        fields = reader.fieldnames or []
        status_field = choose_field(
            fields,
            ["asof_cutoff_status", "row_temporal_status", "temporal_row_status", "row_status"],
            "row status ledger",
        )
        pair_field = choose_field(fields, ["canonical_pair_key", "pair_key"], "row status ledger")
        date_source_field = next(
            (field for field in ["publication_date_source", "date_source"] if field in fields),
            None,
        )
        precision_field = next(
            (field for field in ["publication_date_precision", "date_precision"] if field in fields),
            None,
        )
        for row in reader:
            status = row[status_field].strip()
            if status not in LOCKED_ROW_STATUSES:
                raise ValueError(f"Unrecognized row status outside locked vocabulary: {status!r}")
            date_source = row[date_source_field].strip() if date_source_field else "not_recorded"
            precision = row[precision_field].strip() if precision_field else "not_recorded"
            key = (status, date_source or "missing", precision or "missing")
            counts[key] += 1
            pair_sets[key].add(row[pair_field].strip())
    return [
        {
            "cohort": label,
            "row_status": status,
            "date_source": date_source,
            "date_precision": precision,
            "source_row_count": count,
            "distinct_relation_count": len(pair_sets[(status, date_source, precision)]),
        }
        for (status, date_source, precision), count in sorted(counts.items())
    ]


def main() -> int:
    args = parser().parse_args()
    require_protocol_lock(args.protocol_lock)
    output_dir = require_new_output_dir(args.output_dir)
    cohorts = labeled_paths(args.cohort, "cohort")
    filters = labeled_paths(args.pair_filter, "pair filter")
    status_ledgers = labeled_paths(args.row_status_ledger, "row status ledger")
    unknown_filters = set(filters).difference(cohorts)
    if unknown_filters:
        raise ValueError(f"Pair filters have no matching cohort: {sorted(unknown_filters)}")

    input_records = []
    summaries = []
    internals: dict[str, dict[str, Any]] = {}
    for label, evidence_path in sorted(cohorts.items()):
        input_records.append(input_descriptor(f"{label}_source_evidence", evidence_path))
        allowed_pairs = None
        pair_to_query = None
        if label in filters:
            input_records.append(input_descriptor(f"{label}_relation_filter", filters[label]))
            allowed_pairs, pair_to_query = load_pair_filter(filters[label])
        summary, internal = concentration_summary(
            label=label,
            evidence_path=evidence_path,
            document_type=args.document_type,
            allowed_pairs=allowed_pairs,
            filter_pair_to_query=pair_to_query,
        )
        summaries.append(summary)
        internals[label] = internal

    pairwise = []
    labels = sorted(internals)
    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            left_data = internals[left]
            right_data = internals[right]
            shared_documents = left_data["documents"].intersection(right_data["documents"])
            left_pairs_shared = sum(
                bool(documents.intersection(shared_documents))
                for documents in left_data["pair_documents"].values()
            )
            right_pairs_shared = sum(
                bool(documents.intersection(shared_documents))
                for documents in right_data["pair_documents"].values()
            )
            pairwise.append(
                {
                    "left_cohort": left,
                    "right_cohort": right,
                    "shared_source_document_count": len(shared_documents),
                    "left_source_document_overlap_fraction": (
                        len(shared_documents) / len(left_data["documents"])
                        if left_data["documents"]
                        else None
                    ),
                    "right_source_document_overlap_fraction": (
                        len(shared_documents) / len(right_data["documents"])
                        if right_data["documents"]
                        else None
                    ),
                    "left_relation_count_with_cross_cohort_source": left_pairs_shared,
                    "left_relation_fraction_with_cross_cohort_source": (
                        left_pairs_shared / len(left_data["pairs"]) if left_data["pairs"] else None
                    ),
                    "right_relation_count_with_cross_cohort_source": right_pairs_shared,
                    "right_relation_fraction_with_cross_cohort_source": (
                        right_pairs_shared / len(right_data["pairs"]) if right_data["pairs"] else None
                    ),
                }
            )

    status_rows = []
    for label, path in sorted(status_ledgers.items()):
        input_records.append(input_descriptor(f"{label}_row_status_ledger", path))
        status_rows.extend(audit_row_status(label, path))

    summary_payload = {
        "audit_id": AUDIT_ID,
        "protocol_id": PROTOCOL_ID,
        "source_document_type": args.document_type.upper(),
        "cohorts": summaries,
        "cross_cohort_source_overlap": pairwise,
        "interpretation_boundary": (
            "shared source documents and connected components quantify provenance dependence; "
            "they do not constitute source-document-disjoint external validation"
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    summary_name = "source_concentration_aggregate_summary.json"
    cohort_table_name = "source_concentration_by_cohort.tsv"
    pairwise_name = "cross_cohort_source_overlap.tsv"
    write_json_new(output_dir / summary_name, summary_payload)

    cohort_rows = []
    for item in summaries:
        cohort_rows.append(
            {
                "cohort": item["cohort"],
                "relation_count": item["relation_count"],
                "query_count": item["query_count"],
                "unique_source_document_count": item["unique_source_document_count"],
                "relation_source_document_edge_count": item[
                    "relation_source_document_edge_count"
                ],
                "relation_count_without_source_document": item[
                    "relation_count_without_source_document"
                ],
                "largest_source_document_relation_count": item[
                    "largest_source_document_relation_count"
                ],
                "largest_source_document_relation_fraction": item[
                    "largest_source_document_relation_fraction"
                ],
                "top10_source_document_edge_share": item[
                    "top10_source_document_edge_share"
                ],
                "source_document_edge_hhi": item["source_document_edge_hhi"],
                "source_document_relation_degree_gini": item[
                    "source_document_relation_degree_gini"
                ],
                "relation_source_component_count": item[
                    "relation_source_component_summary"
                ]["component_count"],
                "largest_relation_component_fraction": item[
                    "relation_source_component_summary"
                ]["largest_component_left_node_fraction"],
                "query_source_component_count": item["query_source_component_summary"][
                    "component_count"
                ],
                "largest_query_component_fraction": item[
                    "query_source_component_summary"
                ]["largest_component_left_node_fraction"],
            }
        )
    write_tsv_new(output_dir / cohort_table_name, list(cohort_rows[0]), cohort_rows)
    pairwise_fields = [
        "left_cohort",
        "right_cohort",
        "shared_source_document_count",
        "left_source_document_overlap_fraction",
        "right_source_document_overlap_fraction",
        "left_relation_count_with_cross_cohort_source",
        "left_relation_fraction_with_cross_cohort_source",
        "right_relation_count_with_cross_cohort_source",
        "right_relation_fraction_with_cross_cohort_source",
    ]
    write_tsv_new(output_dir / pairwise_name, pairwise_fields, pairwise)
    output_names = [summary_name, cohort_table_name, pairwise_name]
    if status_rows:
        status_name = "row_date_precision_attrition.tsv"
        write_tsv_new(output_dir / status_name, list(status_rows[0]), status_rows)
        output_names.append(status_name)
    manifest = finalize_manifest(
        output_dir=output_dir,
        audit_id=AUDIT_ID,
        script_path=Path(__file__),
        inputs=input_records,
        output_names=output_names,
        extra={"cohort_count": len(summaries), "row_status_ledger_count": len(status_ledgers)},
    )
    write_json_new(output_dir / "run_manifest.json", manifest)
    print(f"{AUDIT_ID}: wrote aggregate-only summaries for {len(summaries)} cohorts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
