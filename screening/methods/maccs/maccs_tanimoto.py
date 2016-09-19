#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rdkit
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_maccs_tanimoto',
    'representation': 'maccs',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    if len(ligands) == 0:
        raise Exception('Missing ligands!')
    return [MACCSkeys.GenMACCSKeys(molecule)
            for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = MACCSkeys.GenMACCSKeys(query)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
