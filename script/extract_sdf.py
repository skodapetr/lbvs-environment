#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract molecule SDF files for selected tests.
This enable easy integration of the defined benchmarking data with other tools.

Usage example:
 python extract_sdf.py -i 10.1021%2Fci200412p/random_00_05_100_20_3900/5HT1A_Agonist -o ../output/ --flat
 python extract_sdf.py -i 10.1021%2Fci200412p/random_00_05_100_20_3900 -o ../output/ --flat

"""

import os
import json
import argparse
import logging

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

# Test in single groups often share molecule files, so caching can
# save us reading the same files over and over again.
molecules_cache = {}


def read_sdf(path):
    global molecules_cache
    if path in molecules_cache:
        return molecules_cache[path]
    if len(molecules_cache) > 2:
        molecules_cache = {}

    molecules = {}
    mol = ''
    mol_id = ''
    with open(path) as stream:
        for line in stream:
            if mol == '':
                mol_id = line.strip('\n')
            mol += line
            if line.startswith('$$$$'):
                molecules[mol_id] = mol
                mol = ''
    #
    molecules_cache[path] = molecules
    return molecules


def output_path_directory(output, dataset, selection, group, instance, is_flat):
    """

    :param output:
    :param dataset:
    :param selection:
    :param group:
    :param instance:
    :param is_flat:
    :return: Output path.
    """
    if is_flat:
        return output + '/' + dataset + '-' + selection + '-' + group + \
               '/' + instance + '/'
    else:
        return output + '/' + dataset + '/' + selection + '/' + group + \
               '/' + instance + '/'


def input_path_directories(input_split):
    """

    :param input_split:
    :return: Entities to load.
    """
    base_dir = os.path.dirname(os.path.realpath(__file__)) + '/../'
    base_dir += '/data/datasets/' + input_split[0] + \
                '/selections/' + input_split[1]
    if len(input_split) == 3:
        return [{
            'dataset': input_split[0],
            'selection': input_split[1],
            'group': input_split[2],
            'path': base_dir + '/' + input_split[2]
        }]
    # Scan for groups.
    inputs = []
    for group_name in os.listdir(base_dir):
        inputs.append({
            'dataset': input_split[0],
            'selection': input_split[1],
            'group': group_name,
            'path': base_dir + '/' + group_name
        })
    return inputs


def extract(dataset, input_path, output_path, group_actives,
            group_inactives):
    with open(input_path) as f:
        split_data = json.load(f)
        [fn_inactives, fn_actives] = split_data['data']['files']

    sdf_base_path = os.path.dirname(os.path.realpath(__file__)) + '/../' + \
                    '/data/datasets/' + dataset + \
                    '/molecules/sdf/'

    os.makedirs(output_path, exist_ok=True)

    molecules = read_sdf(sdf_base_path + fn_actives + '.sdf')
    fns = [output_path + '/' + fn_actives + '_train.sdf',
           output_path + '/' + fn_actives + '_test.sdf']
    with open(fns[0], 'w') as fta:
        for molecule in split_data['data']['train']['ligands']:
            fta.write(molecules[molecule['name']])
    with open(fns[1], 'w') as fta:
        for molecule in split_data['data']['test']:
            if molecule['name'] in group_actives:
                fta.write(molecules[molecule['name']])

    molecules = read_sdf(sdf_base_path + '/' + fn_inactives + '.sdf')
    fns = [output_path + '/' + fn_inactives + '_train.sdf',
           output_path + '/' + fn_inactives + '_test.sdf']
    with open(fns[0], 'w') as fta:
        for molecule in split_data['data']['train']['decoys']:
            fta.write(molecules[molecule['name']])
    with open(fns[1], 'w') as fta:
        for molecule in split_data['data']['test']:
            if molecule['name'] in group_inactives:
                fta.write(molecules[molecule['name']])


def extract_group(input, output, is_flat):
    # The group file contains definition of activity for the test set.
    with open(input['path'] + '/group.json') as fg:
        group_data = json.load(fg)
        group_actives = set(
            [item['name'] for item in group_data['data']['ligands']])
        group_inactives = set(
            [item['name'] for item in group_data['data']['decoys']])

    for file in os.listdir(input['path']):
        if not file.startswith('s_'):
            continue
        file_name = file[0:file.rfind('.')]
        output_dir = output_path_directory(output, input['dataset'],
                                           input['selection'], input['group'],
                                           file_name, is_flat)
        # Process instance.
        logging.info('Extracting (%s,%s,%s,%s) ...', input['dataset'],
                     input['selection'], input['group'], file_name)
        extract(input['dataset'], input['path'] + '/' + file,
                output_dir, group_actives, group_inactives)

    pass


def main(input, output, flat):
    for item in input_path_directories(input.split('/')):
        extract_group(item, output, flat)
    return


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%H:%M:%S')

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', required=True,
                        help='Input in format: dataset/selection/group'
                             ' or dataset/selection')
    parser.add_argument('-o', dest='output', required=True,
                        help='Root output directory.')
    parser.add_argument('--flat', dest='flat', action='store_true',
                        help='Use flat output directory structure.')
    args = vars(parser.parse_args())

    main(args['input'], args['output'], args['flat'])
