#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rdkit
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Torsions

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_tt_tanimoto',
    'representation': 'tt',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    if len(ligands) == 0:
        raise Exception('Missing ligands!')
    return [Torsions.GetTopologicalTorsionFingerprintAsIntVect(molecule)
            for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = Torsions.GetTopologicalTorsionFingerprintAsIntVect(query)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
