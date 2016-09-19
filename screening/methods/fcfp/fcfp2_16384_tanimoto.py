#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rdkit
from rdkit.Chem import AllChem
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_fcfp2_16384_tanimoto',
    'representation': 'fcfp2_16384',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    return [AllChem.GetMorganFingerprintAsBitVect(molecule, 1, nBits=16384,
                                                  useFeatures=True)
            for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(query, 1, nBits=16384,
                                                        useFeatures=True)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
