#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Implementation of fingerprints as presented as HI-DS 2015.
Fragments:
    - paths of size 3
Descriptors:
    - field descriptor (3-path) hybridization
Fingerprint:
    - count vector, no hashing
"""

import rdkit
from rdkit import Chem

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    'id': 'method_rdkit_p3_f3_hybridization_tanimoto_like',
    'representation': 'p3_f3_hybridization',
    'similarity': 'tanimoto_like'
}


# region Representation calculation

def _pattern_fragment_generator(molecule, pattern):
    """

    :param molecule:
    :param pattern:
    :return: A generator of fragment for given molecule.
    """
    matches = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
    for atoms in matches:
        yield atoms


def _calculate_atom_code(atom_idx, field_map):
    """Create labeling for given atom from the given fields.

    :param atom_idx:
    :param field_map:
    :return:
    """
    # There may not be a value for every atom, in such case we will use '_'.
    # This might be as we use bigger fragment to build field then to
    # construct fingerprint.
    if atom_idx in field_map:
        return field_map[atom_idx]
    else:
        return None


def _calculate_fragment_code(atoms, field_map):
    """Compute fragment code.

    :param atoms:
    :param field_map:
    :return: Score for given fragment as a vector of values.
    """
    codes = [""] * len(atoms)
    for i in range(0, len(atoms)):
        codes[i] = _calculate_atom_code(atoms[i], field_map)
    # canonization of the code vector
    beg = 0
    end = len(codes) - 1
    while beg < end:
        if codes[beg] > codes[end]:
            codes.reverse()
            break
        elif codes[beg] == codes[end]:
            beg += 1
            end -= 1
        else:
            break
    return codes


def _calculate_field_hybrid_descriptor(molecule, fragments):
    """Create field hybridization descriptor for each atom.

    :param molecule:
    :param fragments:
    :return: Dictionary <atom, label value>
    """
    value_map = {}
    counter_map = {}
    for atoms in fragments:
        # Compute average field value.
        average_value = 0.0
        for atom in atoms:
            average_value += molecule.GetAtomWithIdx(atom).GetHybridization()
        average_value /= len(atoms)
        # Save value to atoms.
        for atom in atoms:
            if atom not in value_map:
                value_map[atom] = 0
                counter_map[atom] = 0
            value_map[atom] += average_value
            counter_map[atom] += 1
    # Modify values as averages.
    for key in value_map.keys():
        value_map[key] /= float(counter_map[key])
    return value_map


def compute_representation(molecule):
    """Calculate and return representation of given molecule.
    The output is dictionary of objects, where each object has a value (id)
    and counter, which represent number of fragment instances.

    :param molecule:
    :return: Dictionary <label, {value, counter}>
    """
    fragments_field = list(_pattern_fragment_generator(molecule, "*~*~*"))
    field_map = _calculate_field_hybrid_descriptor(molecule, fragments_field)
    fragments = fragments_field
    representation = {}
    total_counter = 0
    for atoms in fragments:
        fragment = _calculate_fragment_code(atoms, field_map)
        fragment_as_str = str(fragment)
        if fragment_as_str not in representation:
            representation[fragment_as_str] = {
                'value': fragment,
                'instances': [],
                'counter': 0
            }
        representation[fragment_as_str]['counter'] += 1
        # Contains information about atoms, that form the fragment.
        representation[fragment_as_str]['instances'].append({
            'atoms': atoms
        })
        total_counter += 1
    return representation


def _fragment_score_label_similarity(left_fragment, right_fragment):
    """Given two fragments return their similarity.

    :param left_fragment:
    :param right_fragment:
    :return:
    """
    # Topology must be the same.
    if not left_fragment['value'] == right_fragment['value']:
        return 0
    similarity = 0
    for index in range(1, len(left_fragment)):
        try:
            similarity += left_fragment[index] == right_fragment[index]
        except Exception as ex:
            raise ex
    return similarity


def similarity_tannimoto(left, right):
    """

    :param left: Dictionary <key, value>
    :param right: Dictionary <key, value>
    :return:
    """
    in_common = 0
    total = 0
    for left_fragment_key in left:
        for right_fragment_key in right:
            left_fragment = left[left_fragment_key]
            right_fragment = right[right_fragment_key]
            if not left_fragment['value'] == right_fragment['value']:
                continue
            # Add shared number (ie. min) to in_common - intersenction.
            in_common += min(left_fragment['counter'],
                             right_fragment['counter'])
            # Add to total counter - union.
            total += left_fragment['counter'] + right_fragment['counter']
    if not total == 0:
        # They have no fragments, this might happen if the molecules
        # are too small..
        return in_common / (float)(total)
    else:
        return 0


# endregion

# region Interface method implementation

def create_model(ligands, decoys):
    """Create and return query model.
    """
    return [compute_representation(molecule) for molecule in ligands]


def compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = compute_representation(query)
    return max([similarity_tannimoto(fingerprint, ligand)
                for ligand in model])

# endregion
