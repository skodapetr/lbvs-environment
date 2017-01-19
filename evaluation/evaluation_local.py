#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Evaluate results of the screening.

Iterate over results in data/screening and perform evaluation. The output is
stored into data/evaluation.
"""

import os
import argparse
import json
import logging
import csv
import time
import functools
from rdkit.ML.Scoring import Scoring

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'


# region Loading files

def find_json_files(root):
    """Recursively find all JSON files.

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


def find_results():
    """Search and return paths to all screening results.

    :return:
    """
    logging.info('Searching for results ...')
    data_path = os.path.dirname(os.path.realpath(__file__)) + \
                '/../data/screening/'
    results = find_json_files(data_path)
    logging.info('Searching for results ... done (%d)', len(results))
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
            group = ''
            if 'id' in data['metadata']:
                group = data['metadata']['id']
            if 'group' in data['metadata']:
                group = data['metadata']['group']
            id = data['metadata']['dataset'] + data['metadata']['selection'] + \
                 group
            results[id] = data
    logging.info('Loading groups ... done (%d)', len(results))
    return results


# endregion

# region CSV read and write

def line_to_entry(columns, row):
    """Convert and return line CSV object into dictionary.

    :param columns:
    :param row:
    :return:
    """
    output = {}
    for i in range(0, len(columns)):
        output[columns[i]] = row[i]
    return output


def header():
    """Return a header line.

    :return:
    """
    return [
        'path', 'dataset', 'selection', 'group', 'instance', 'method',
        'actives', 'decoys',
        'w-auc', 'w-ef-0.005', 'w-ef-0.01', 'w-ef-0.02', 'w-ef-0.05',
        'w-rie-0.005', 'w-rie-0.01', 'w-rie-0.02', 'w-rie-0.05',
        'w-bedroc-0.005', 'w-bedroc-0.01', 'w-bedroc-0.02', 'w-bedroc-0.05',
        'b-auc', 'b-ef-0.005', 'b-ef-0.01', 'b-ef-0.02', 'b-ef-0.05',
        'b-rie-0.005', 'b-rie-0.01', 'b-rie-0.02', 'b-rie-0.05',
        'b-bedroc-0.005', 'b-bedroc-0.01', 'b-bedroc-0.02', 'b-bedroc-0.05'
    ]


def entry_to_line(entry):
    """Convert entry to the line.

    :param entry:
    :return:
    """
    values = []
    for key in header():
        values.append(entry[key])
    return [str(item) for item in values]


# endregion

# region Score comparison

def score_compare_worst(left, right):
    if left['score'] == right['score']:
        # Decide based on activity.
        if left['activity'] and right['activity']:
            return 0
        elif left['activity']:
            # Only left is active.
            return -1
        else:
            return 1
    else:
        if left['score'] > right['score']:
            return 1
        else:
            return -1


def score_compare_best(left, right):
    if left['score'] == right['score']:
        # Decide based on activity.
        if left['activity'] and right['activity']:
            return 0
        elif left['activity']:
            # Only left is active.
            return 1
        else:
            return -1
    else:
        if left['score'] > right['score']:
            return 1
        else:
            return -1


# endregion

def compute_measures(scores, params=[0.005, 0.01, 0.02, 0.05]):
    """Evaluate scores and return object with performance metrics,

    :param scores:
    :param params:
    :return:
    """
    output = {
        'auc': Scoring.CalcAUC(scores, 'activity')
    }
    for value in params:
        output['ef-' + str(value)] = \
            Scoring.CalcEnrichment(scores, 'activity', [value])[0]
        output['bedroc-' + str(value)] = \
            Scoring.CalcBEDROC(scores, 'activity', value)
        output['rie-' + str(value)] = \
            Scoring.CalcRIE(scores, 'activity', value)
    return output


def evaluate_result_file(result_path, groups):
    """Evaluate a single result file and return the result object.

    result_path:
    :param groups:
    :return: None in case of an error.
    """
    with open(result_path) as stream:
        try:
            result = json.load(stream)
        except:
            logging.exception('Invalid result file: %s', result_path)
            return None
    group_id = result['dataset'] + result['selection'] + result['group']
    if group_id not in groups:
        logging.warning('Missing group.')
        return None
    group = groups[group_id]
    ligand_names = set()
    for item in group['data']['ligands']:
        ligand_names.add(item['id'])
    decoys_names = set()
    for item in group['data']['decoys']:
        decoys_names.add(item['id'])
    # Add activity to molecules.
    missing_activity = False
    actives_count = 0
    decoys_count = 0
    for score in result['scores']:
        is_ligand = score['molecule'] in ligand_names
        is_decoy = score['molecule'] in decoys_names
        if is_ligand and not is_decoy:
            score['activity'] = 1
            actives_count += 1
        elif not is_ligand and is_decoy:
            score['activity'] = 0
            decoys_count += 1
        else:
            logging.error('Missing: "%s" in "%s"',
                          score['molecule'], result_path)
            missing_activity = True
    if missing_activity:
        return None
    # Evaluate.
    eval_best = compute_measures(sorted(result['scores'],
                                        key=functools.cmp_to_key(
                                            score_compare_best),
                                        reverse=True))
    eval_worst = compute_measures(sorted(result['scores'],
                                         key=functools.cmp_to_key(
                                             score_compare_worst),
                                         reverse=True))
    # Store results to object.
    return {
        'group': result['group'],
        'dataset': result['dataset'],
        'selection': result['selection'],
        'instance': result['instance'],
        'method': result['method'],
        'actives': actives_count,
        'decoys': decoys_count,
        'w-auc': eval_worst['auc'],
        'w-ef-0.005': eval_worst['ef-0.005'],
        'w-ef-0.01': eval_worst['ef-0.01'],
        'w-ef-0.02': eval_worst['ef-0.02'],
        'w-ef-0.05': eval_worst['ef-0.05'],
        'w-rie-0.005': eval_worst['rie-0.005'],
        'w-rie-0.01': eval_worst['rie-0.01'],
        'w-rie-0.02': eval_worst['rie-0.02'],
        'w-rie-0.05': eval_worst['rie-0.05'],
        'w-bedroc-0.005': eval_worst['bedroc-0.005'],
        'w-bedroc-0.01': eval_worst['bedroc-0.01'],
        'w-bedroc-0.02': eval_worst['bedroc-0.02'],
        'w-bedroc-0.05': eval_worst['bedroc-0.05'],
        'b-auc': eval_best['auc'],
        'b-ef-0.005': eval_best['ef-0.005'],
        'b-ef-0.01': eval_best['ef-0.01'],
        'b-ef-0.02': eval_best['ef-0.02'],
        'b-ef-0.05': eval_best['ef-0.05'],
        'b-rie-0.005': eval_worst['rie-0.005'],
        'b-rie-0.01': eval_worst['rie-0.01'],
        'b-rie-0.02': eval_worst['rie-0.02'],
        'b-rie-0.05': eval_worst['rie-0.05'],
        'b-bedroc-0.005': eval_worst['bedroc-0.005'],
        'b-bedroc-0.01': eval_worst['bedroc-0.01'],
        'b-bedroc-0.02': eval_worst['bedroc-0.02'],
        'b-bedroc-0.05': eval_worst['bedroc-0.05']
    }


def read_results(path):
    """Read and return existing results.

    :param path:
    :return:
    """
    results = []
    paths = set()
    with open(path) as stream:
        csv_reader = csv.reader(stream, delimiter=',', quotechar='"')
        first_line = True
        for row in csv_reader:
            if first_line:
                columns = row
                first_line = False
                continue
            entry = line_to_entry(columns, row)
            results.append(entry)
            paths.add(entry['path'])
    return results, paths


def load_configuration():
    """Load and return the configuration.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Evaluate runs on local machine.')
    parser.add_argument('-t', type=int, dest='time', required=False, default=-1,
                        help='Time limit in minutes.')
    return vars(parser.parse_args())


def main():
    configuration = load_configuration()
    time_end = -1
    if not configuration['time'] == -1:
        time_end = time.time() + (configuration['time'] * 60)
    # Search for group files and results.
    results = find_results()
    groups = load_groups()
    # Load existing results if they exists.
    eval_dir = os.path.dirname(os.path.realpath(__file__)) + \
               '/../data/evaluation/'
    if not os.path.exists(eval_dir):
        os.makedirs(eval_dir)
    results_path = eval_dir + 'results.csv'
    eval_list = []
    eval_paths = set()
    if os.path.exists(results_path):
        logging.info('Reading results ...')
        eval_list, eval_paths = read_results(results_path)
        logging.info('Reading results ... done')
    # Evaluation.
    logging.info('Evaluating ...')
    with open(results_path, 'w') as output_stream:
        output_stream.write('"')
        output_stream.write('","'.join(header()))
        output_stream.write('"\n')
        # Write old results.
        for item in eval_list:
            output_stream.write('"')
            output_stream.write('","'.join(entry_to_line(item)))
            output_stream.write('"\n')
        eval_list = []
        # Iterate over new results
        for result_path in results:
            relative_path = result_path[
                            result_path.index('/../data/screening/') + 19:]
            # Use path to determine if the file was already processed or not.
            if relative_path in eval_paths:
                continue
            result_entry = evaluate_result_file(result_path, groups)
            if result_entry is None:
                continue
            result_entry['path'] = relative_path
            # Write to file.
            output_stream.write('"')
            output_stream.write('","'.join(entry_to_line(result_entry)))
            output_stream.write('"\n')
            if -1 < time_end < time.time():
                logging.info('timeout')
                break
    logging.info('Evaluating ... done')


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%H:%M:%S')
    main()
