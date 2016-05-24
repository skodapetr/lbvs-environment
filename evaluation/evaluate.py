#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Compute performance measures for results of LBVS.

Use file system to create locks for processed data, Thanks to that
it can be executed multiple times in parallel.
"""

import os
import re
import json
import logging
import shutil
from rdkit.ML.Scoring import Scoring

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

_root_path = os.path.dirname(os.path.realpath(__file__)) + '/../'


def read_json_file(path):
    """Read and return content of given JSON file.

    :param path:
    :return:
    """
    with open(path, 'r') as stream:
        return json.load(stream)


def recursive_scan(directory_path):
    """Scan and return all *.json files in directory recursively.

    :param directory_path:
    :return:
    """
    result = []
    for name in os.listdir(directory_path):
        path = directory_path + '/' + name
        if os.path.isdir(path):
            result.extend(recursive_scan(path))
            continue
        if os.path.isfile(path) and path.endswith('.json'):
            result.append(path)
    return result


def parse_path(output_root, input_path):
    """Parse path for information about method, dataset, selection.

    :param output_root:
    :param input_path:
    :return:
    """
    name = re.split(r'/|\\', input_path[len(output_root):])
    name = [x for x in name if not x == '']
    return {
        'method': name[0],
        'dataset': name[1],
        'selection': name[2],
        'target': name[3][0:name[3].rfind('_')]
    }


def evaluate_file(file_path, dataset_path, context):
    """Compute and return evaluation for given file.

    :param file_path:
    :param dataset_path:
    :param context:
    :return:
    :return:
    """
    with open(file_path, 'r') as stream:
        try:
            score = json.load(stream)
        except Exception as ex:
            print('Error: ', file_path)
            raise ex
    # Load names of ligands and decoys.
    name_path = dataset_path + \
                score['selection']['dataset'] + \
                '/molecules/name/' + \
                score['selection']['target']

    if name_path not in context:
        context[name_path] = {
            'ligands': read_json_file(name_path + '_Ligands.json'),
            'decoys': read_json_file(name_path + '_Decoys.json')
        }
    # Set activity.
    ligands = context[name_path]['ligands']
    decoys = context[name_path]['decoys']
    for molecule_score in score['scores']:
        is_ligand = molecule_score['molecule'] in ligands
        is_decoy = molecule_score['molecule'] in decoys
        if is_ligand and not is_decoy:
            molecule_score['activity'] = True
        elif not is_ligand and is_decoy:
            molecule_score['activity'] = False
        else:
            raise Exception('Molecule is not active nor decoy: '
                            + molecule_score['molecule'] + ' path: '
                            + file_path)
    if len(score['scores']) == 0:
        logging.warning('Empty file: %s', file_path)
        raise Exception()

    # Sort based on the similarity.
    sorted_score = sorted(score['scores'],
                          key=lambda m: float(m['score']),
                          reverse=True)
    # Compute statistics.
    params = [0.005, 0.01, 0.02, 0.05]
    auc = Scoring.CalcAUC(sorted_score, 'activity')
    ef = Scoring.CalcEnrichment(sorted_score, 'activity', params)
    ef_dict = {}
    for index in range(0, len(ef)):
        ef_dict[params[index]] = ef[index]

    bedroc = {}
    for value in params:
        bedroc[value] = Scoring.CalcBEDROC(sorted_score, 'activity', value)
    rie = {}
    for value in params:
        rie[value] = Scoring.CalcRIE(sorted_score, 'activity', value)
    # Write result to a file.
    return {
        'method': score['method'],
        'selection': score['selection'],
        'metadata': {
            'tested_molecules': len(score['scores'])
        },
        'AUC': auc,
        'EF': ef_dict,
        'BEDROC': bedroc,
        'RIE': rie
    }


def evaluate_target(evaluation_root, input_root, dataset_path, files):
    """Evaluate given files for first discovered target.

    :param evaluation_root: Evaluation directory where output files are stored.
    :param input_root: Prefix of files, used to determine relative path.
    :param dataset_path:
    :param files:
    :return: List of files that are in evaluated state after this method.
    """
    output_path = None
    locked_info = None
    processed = []
    context = {}
    data = []
    for file_path in files:
        output_info = parse_path(input_root, file_path)
        # Try to lock the directory.
        if locked_info is None:
            output_path = evaluation_root + '/' + \
                          output_info['method'] + '/' + \
                          output_info['dataset'] + '/' + \
                          output_info['selection'] + '/' + \
                          output_info['target'] + '/'
            try:
                os.makedirs(output_path)
                locked_info = output_info
                logging.info('Evaluating: ' +
                             output_info['dataset'] + '/' +
                             output_info['selection'] + '/' +
                             output_info['target'])
            except OSError:
                # Already taken.
                processed.append(file_path)
                continue
        else:
            if not output_info == locked_info:
                # We are not working on this one now.
                continue
        # Load, process and save.
        data.append(evaluate_file(file_path, dataset_path, context))
        processed.append(file_path)
    # Store all data into as single file.
    with open(output_path + '/results.json', 'w') as stream:
        json.dump(data, stream)
    logging.info('Evaluation done')
    return processed


def _main():
    input_root = _root_path + '/results/output/'
    evaluation_root = _root_path + '/results/evaluation/'
    dataset_path = _root_path + '/datasets/'
    #
    logging.info('Searching files ...')
    files = recursive_scan(input_root)
    logging.info('Searching files ... done : %d', len(files))
    #
    while True:
        processed = evaluate_target(evaluation_root,
                                    input_root,
                                    dataset_path, files)
        for x in processed:
            files.remove(x)
        if len(processed) == 0 or len(files) == 0:
            break


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    #
    _main()
