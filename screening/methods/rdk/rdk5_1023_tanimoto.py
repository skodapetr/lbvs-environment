#!/usr/bin/env python
# -*- coding: utf-8 -*-

import rdkit
from rdkit import DataStructs
from rdkit import Chem

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_rdk5_1023_tanimoto',
    'representation': 'rdk5_1023',
    'similarity': 'tanimoto'
}


# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    if len(ligands) == 0:
        raise Exception('Missing ligands!')
    return [Chem.RDKFingerprint(molecule, maxPath=5, fpSize=1023,
                                nBitsPerHash=2) for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = Chem.RDKFingerprint(query, maxPath=5, fpSize=1023,
                                      nBitsPerHash=2)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])

# endregion
