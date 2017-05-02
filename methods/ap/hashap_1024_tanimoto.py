#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

__license__ = "X11"

METADATA = {
    "id": "method_rdkit_hashap_1024_tanimoto",
    "representation": "hashap_1024",
    "similarity": "tanimoto"
}


def _compute_fingerprint(molecule):
    return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
        molecule, nBits=1024)


def _compute_similarity(left, right):
    return DataStructs.TanimotoSimilarity(left, right)


def create_model(train_ligands, train_decoys):
    model = []
    for molecule in train_ligands:
        model.append({
            "name": molecule.GetProp("_Name"),
            "fingerprint": _compute_fingerprint(molecule)
        })
    model_information = {}
    return model, model_information


def compute_score(model, molecule):
    fingerprint = _compute_fingerprint(molecule)
    similarities = [_compute_similarity(fingerprint, item["fingerprint"])
                    for item in model]
    max_score = max(similarities)
    index_of_max_score = similarities.index(max_score)
    closest_molecule = model[index_of_max_score]
    return {
        "value": max_score,
        "info": {
            "closest": closest_molecule["name"]
        }
    }


def compute_similarity(left, right):
    return _compute_similarity(_compute_fingerprint(left),
                               _compute_fingerprint(right))
