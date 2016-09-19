#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Iterate over results in data/screening and perform evaluation.

The output is stored into data/evaluation.
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


def find_json_files(root):
    """Find all JSON files recursively.

    :param root:
    :return:
    """
    result = []
    for file_name in os.listdir(root):
        file_path = root + '/' + file_name
        if os.path.isdir(file_path):
            result.extend(find_json_files(file_path))
            continue
        if file_name.endswith('.json'):
            result.append(file_path)
    return result


def find_selection_files(root):
    """Find all JSON selection files recursively.

    :param root:
    :return:
    """
    logging.info('Searching for selections ...')
    result = []
    for file_name in os.listdir(root):
        file_path = root + '/' + file_name
        if os.path.isdir(file_path):
            result.extend(find_selection_files(file_path))
            continue
        if file_name.endswith('.json') and file_name.startswith('selection_'):
            result.append(file_path)
    logging.info('Searching for selections ... done')
    return result


def find_results():
    """Search and return paths to all screening results.

    :return:
    """
    logging.info('Searching for results ...')
    data_path = os.path.dirname(os.path.realpath(__file__)) + \
                '/../data/screening/'
    results = find_json_files(data_path)
    logging.info('Searching for results ... done')
    return results


def load_groups():
    """Find load and return all groups.

    :return:
    """
    logging.info('Loading groups ...')
    datasets_path = os.path.dirname(os.path.realpath(__file__)) + \
                    '/../data/datasets'
    paths = [file for file in find_json_files(datasets_path)
             if file.endswith('group.json')]

    results = {}
    for path in paths:
        with open(path) as stream:
            data = json.load(stream)
            id = data['metadata']['dataset'] + \
                 data['metadata']['selection'] + \
                 data['metadata']['id']
            results[id] = data

    logging.info('Loading groups ... done')
    return results


def evaluate(sorted_scores):
    params = [0.005, 0.01, 0.02, 0.05]
    auc = Scoring.CalcAUC(sorted_scores, 'activity')
    ef = {}
    bedroc = {}
    rie = {}
    for value in params:
        ef[value] = Scoring.CalcEnrichment(sorted_scores, 'activity', [value])[
            0]
        bedroc[value] = Scoring.CalcBEDROC(sorted_scores, 'activity', value)
        rie[value] = Scoring.CalcRIE(sorted_scores, 'activity', value)
    return {
        'auc': auc,
        'ef': ef,
        'bedroc': bedroc,
        'rie': rie
    }


def toCsvLine(obj):
    values = [obj['dataset'], obj['selection'], obj['group'],
              obj['instance'], obj['method'], obj['auc'],
              obj['ef'][0.01], obj['ef'][0.02], obj['ef'][0.05]]
    return [str(item) for item in values]


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    # Search for group files and results.
    results = find_results()
    # selections = find_selections()
    groups = load_groups()
    #
    logging.info('Evaluating ...')
    eval_list = []
    for result_path in results:
        with open(result_path) as stream:
            result = json.load(stream)
        group_id = result['dataset'] + \
                   result['selection'] + \
                   result['group']
        group = groups[group_id]
        # group_path = find_group_definition(result['selection'])
        # with open(group_path) as stream:
        #     group = json.load(stream)
        #
        ligand_names = set()
        for item in group['data']['ligands']:
            ligand_names.add(item['id'])
        decoys_names = set()
        for item in group['data']['decoys']:
            decoys_names.add(item['id'])
        # Add activity to molecules.
        for score in result['scores']:
            is_ligand = score['molecule'] in ligand_names
            is_decoy = score['molecule'] in decoys_names
            if is_ligand and not is_decoy:
                score['activity'] = 1
            elif not is_ligand and is_decoy:
                score['activity'] = 0
            else:
                raise Exception('Missing record for molecule: ' + score['name'])
        # Sort based on the similarity.
        sorted_scores = sorted(result['scores'],
                               key=lambda m: float(m['score']),
                               reverse=True)
        # Evaluate and store the result.
        eval = evaluate(sorted_scores)
        eval.update({
            'group': result['group'],
            'dataset': result['dataset'],
            'selection': result['selection'],
            'instance': result['instance'],
            'method': result['method']
        });
        eval_list.append(eval)
    logging.info('Evaluating ... done')
    # Write the results CSV file.
    logging.info('Writing results ...')
    header = ['dataset', 'selection', 'group', 'instance', 'method', 'auc',
              'ef-0.01', 'ef-0.02', 'ef-0.05']

    eval_dir = os.path.dirname(os.path.realpath(__file__)) + \
               '/../data/evaluation/'
    os.makedirs(eval_dir, exist_ok=True)
    with open(eval_dir + 'results.csv', 'w') as stream:
        stream.write('"')
        stream.write('","'.join(header))
        stream.write('"\n')
        for item in eval_list:
            stream.write('"')
            stream.write('","'.join(toCsvLine(item)))
            stream.write('"\n')
    logging.info('Writing results ... done')


if __name__ == '__main__':
    main()
