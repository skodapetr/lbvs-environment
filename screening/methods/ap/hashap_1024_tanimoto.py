#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_hashap_1024_tanimoto',
    'representation': 'hashap_1024',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    if len(ligands) == 0:
        raise Exception('Missing ligands!')
    return [rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
        molecule, nBits=1024) for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
        query, nBits=1024)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
