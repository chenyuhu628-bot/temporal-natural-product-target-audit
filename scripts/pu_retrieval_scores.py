"""Feature construction and leakage-safe score transfer for PU baselines."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn.feature_extraction.text import TfidfVectorizer


def morgan_fingerprints(smiles_by_compound: dict[str, str], radius: int, n_bits: int) -> dict[str, object]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fingerprints: dict[str, object] = {}
    for compound, smiles in smiles_by_compound.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse locked baseline structure: {compound}")
        fingerprints[compound] = generator.GetFingerprint(molecule)
    return fingerprints


def build_train_maps(
    train_pairs: list[dict[str, str]], target_index: dict[str, int], weights: dict[str, float],
) -> tuple[dict[str, list[tuple[int, float]]], np.ndarray]:
    by_compound: dict[str, list[tuple[int, float]]] = defaultdict(list)
    popularity = np.zeros(len(target_index), dtype=np.float32)
    for row in train_pairs:
        target = row["uniprot_canonical_accession"]
        if target not in target_index:
            raise ValueError(f"Historical target absent from candidate universe: {target}")
        weight = weights[row["best_strict_evidence_tier_v1_1"]]
        target_idx = target_index[target]
        by_compound[row["inchikey_full"]].append((target_idx, weight))
        popularity[target_idx] += weight
    return dict(by_compound), popularity


def tanimoto_transfer_scores(
    query_compound: str,
    train_by_compound: dict[str, list[tuple[int, float]]],
    fingerprints: dict[str, object],
    target_count: int,
) -> np.ndarray:
    scores = np.zeros(target_count, dtype=np.float32)
    train_compounds = list(train_by_compound)
    if not train_compounds:
        return scores
    similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[query_compound], [fingerprints[item] for item in train_compounds])
    for compound, similarity in zip(train_compounds, similarities):
        for target_idx, weight in train_by_compound[compound]:
            scores[target_idx] = max(scores[target_idx], float(similarity) * weight)
    return scores


def build_sequence_kmer_matrix(target_ids: list[str], sequences: dict[str, str]) -> tuple[TfidfVectorizer, object]:
    if set(target_ids).difference(sequences):
        raise ValueError("Candidate target catalogue has an accession without a C31 sequence")
    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(3, 3), lowercase=False, norm="l2", dtype=np.float32)
    return vectorizer, vectorizer.fit_transform([sequences[target] for target in target_ids])


def sequence_transfer_scores(
    query_compound: str,
    train_by_compound: dict[str, list[tuple[int, float]]],
    sequence_matrix: object,
    target_count: int,
) -> np.ndarray:
    scores = np.zeros(target_count, dtype=np.float32)
    known = train_by_compound.get(query_compound, [])
    if not known:
        return scores
    indices = [target for target, _ in known]
    weights = np.asarray([weight for _, weight in known], dtype=np.float32)
    similarities = (sequence_matrix @ sequence_matrix[indices].T).toarray().astype(np.float32, copy=False)
    return np.max(similarities * weights[np.newaxis, :], axis=1)


def software_versions() -> dict[str, str]:
    import numpy
    import sklearn
    return {"rdkit": rdBase.rdkitVersion, "numpy": numpy.__version__, "scikit_learn": sklearn.__version__}

