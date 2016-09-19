#!/usr/bin/env python
# -*- coding: utf-8 -*-

from rdkit.Chem.AtomPairs import Pairs
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_ap_tanimoto',
    'representation': 'ap',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    if len(ligands) == 0:
        raise Exception('Missing ligands!')
    return [Pairs.GetAtomPairFingerprint(molecule)
            for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = Pairs.GetAtomPairFingerprint(query)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
