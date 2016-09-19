#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rdkit
from rdkit.Chem import AllChem
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_ecfc6_tanimoto',
    'representation': 'ecfc6',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    return [AllChem.GetMorganFingerprint(molecule, 3)
            for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = AllChem.GetMorganFingerprint(query, 3)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
