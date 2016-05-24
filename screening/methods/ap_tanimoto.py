#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import json
import rdkit
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

metadata = {
    # Must be same as the file name without an extension.
    'name': 'ap_tanimoto'
}


# region Similarity method implementation

def _create_model(ligands, decoys):
    """Create and return query model.
    """
    return [rdMolDescriptors.GetAtomPairFingerprint(molecule)
            for molecule in ligands]


def _compute_score(model, query):
    """Compute and return activity likeness for given query molecule.
    """
    fingerprint = rdMolDescriptors.GetAtomPairFingerprint(query)
    return max([DataStructs.TanimotoSimilarity(fingerprint, ligand)
                for ligand in model])


# endregion

# region Shared code

def screening(selection, molecules):
    # Select train set.
    test_ligands = [molecules[item['name']]
                    for item in selection['data']['train']['ligands']
                    if item['name'] in molecules]
    test_decoys = [molecules[item['name']]
                   for item in selection['data']['train']['decoys']
                   if item['name'] in molecules]
    model = _create_model(test_ligands, test_decoys)
    # Perform screening.
    logging.info('Screening ...')
    scores = []
    for item in selection['data']['test']:
        if not item['name'] in molecules:
            continue
        scores.append({
            'molecule': item['name'],
            'score': _compute_score(model, molecules[item['name']])
        })
    logging.info('Screening ... done')
    # Write results.
    result = {
        'selection': selection['info'],
        'method': metadata,
        'scores': scores
    }
    return result


def load_and_screen(data_directory, input_file, output_path):
    with open(input_file, 'r') as stream:
        selection = json.load(stream)
    source_files = [data_directory + selection['info']['dataset'] +
                    '/molecules/sdf/' + x + '.sdf'
                    for x in selection['files']]
    molecules = {}
    logging.info('Loading molecules ...')
    for path in source_files:
        for molecule in rdkit.Chem.SDMolSupplier(str(path)):
            if molecule is None:
                logging.error("Can't load molecule.")
                continue
            molecules[molecule.GetProp('_Name')] = molecule
    logging.info('Loading molecules ... done : %d', len(molecules))
    result = screening(selection, molecules)
    #
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    with open(output_path, 'w') as stream:
        json.dump(result, stream)
    logging.info('Results have been writen')


def _read_configuration():
    parser = argparse.ArgumentParser(
        description='Perform a simple virtual screening.')
    parser.add_argument('-i', type=str, dest='input',
                        help='Path to selection file.',
                        required=True)
    parser.add_argument('-d', type=str, dest='dataset',
                        help='Path to root fo the dataset directory.',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='Path to output file.',
                        required=True)
    return vars(parser.parse_args())


def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    config = _read_configuration()
    load_and_screen(config['dataset'], config['input'], config['output'])


if __name__ == '__main__':
    _main()

# endregion
