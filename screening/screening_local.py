#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Perform screening on the local machine.

"""

import argparse
import logging
import os
import json

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'


def import_methods(methods_names):
    """Import given methods specified by method names an return them.

    :param methods_names:
    :return:
    """
    import imp
    method_directory = os.path.dirname(os.path.realpath(__file__)) + '/methods/'
    methods_modules = {}
    for method in methods_names:
        module = imp.load_source(method, method_directory + method)
        methods_modules[method] = module
    return methods_modules


def load_molecules(selection_path, selection, cache):
    """Load molecules for given selection into cache.

    :param selection_path:
    :param selection:
    :param cache:
    :return:
    """
    import rdkit
    molecule_path = selection_path[0:selection_path.rfind('selections')] \
                    + '/molecules/sdf/'
    molecules_file = selection['data']['files']
    source_files = [molecule_path + item + '.sdf' for item in molecules_file]
    # We use caching on the level of set of files.
    source_files.sort()
    if cache['files'] == source_files:
        return
    logging.info('Loading molecules ... ')
    cache['files'] = source_files
    molecules = {}
    for path in source_files:
        for molecule in rdkit.Chem.SDMolSupplier(str(path)):
            if molecule is None:
                logging.error("Can't load molecule.")
                continue
            molecules[molecule.GetProp('_Name')] = molecule
    # Store results.
    cache['molecules'] = molecules
    logging.info('Loading molecules ... done : %d', len(molecules))


def screen_selection_method(module, output_path, selection, cache):
    """Use methods from given module to perform LBVS.

    :param module:
    :param output_path:
    :param selection:
    :param cache:
    :return:
    """
    # Perform screening.
    logging.info('Screening ...')
    # def screening(selection, molecules):
    molecules = cache['molecules']
    # Select train set.
    test_ligands = [molecules[item['id']]
                    for item in selection['data']['train']['ligands']
                    if item['id'] in molecules]
    test_decoys = [molecules[item['id']]
                   for item in selection['data']['train']['decoys']
                   if item['id'] in molecules]
    model = module.create_model(test_ligands, test_decoys)
    # Perform screening.
    logging.info('Screening ...')
    scores = []
    for item in selection['data']['test']:
        if not item['id'] in molecules:
            continue
        scores.append({
            'molecule': item['id'],
            'score': module.compute_score(model, molecules[item['id']])
        })
    logging.info('Screening ... done')
    # Write results.
    results = {
        'dataset': selection['metadata']['dataset'],
        'instance': selection['metadata']['id'],
        'selection': selection['metadata']['selection'],
        'group': selection['metadata']['group'],
        'method': module.metadata['id'],
        'scores': scores
    }
    # return result
    # results = method.screening(selection, cache['molecules'])
    logging.info('Screening ... done')
    # Save results
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    with open(output_path, 'w') as stream:
        json.dump(results, stream)


def screening(methods, selections):
    """Perform screening be importing the Python script.

    :param methods:
    :param selections:
    :return:
    """
    methods_modules = import_methods(methods)
    # Use caching for working with molecules.
    cache = {
        'molecules': {},
        # Cache for loaded files.
        'files': {}
    }
    data_directory = os.path.dirname(os.path.realpath(__file__)) + '/../data/'
    for selection_path in selections:
        logging.info('Selection: %s', selection_path)
        # TODO Add error check here
        with open(selection_path, 'r') as stream:
            selection = json.load(stream)
        # TODO Check for output, so we do not screen the same data twice !
        load_molecules(selection_path, selection, cache)
        for method_name in methods:
            # Remove the .py extension from name before use.
            id = method_name[0:-3] + '/' + selection['metadata']['id']
            output_path = data_directory + '/screening/' + id + '.json'
            lock_path = data_directory + '/lock/' + id
            # Use directory based lock.
            try:
                os.makedirs(lock_path)
            except Exception as ex:
                # Other thread is computing this representation.
                return
            #
            try:
                screen_selection_method(methods_modules[method_name],
                                        output_path, selection, cache)
            finally:
                # Release lock.
                os.rmdir(lock_path)
    pass


def find_selection_files(root):
    """Find all JSON selection files recursively.

    :param root: Can be directory or file.
    :return:
    """
    if os.path.isfile(root):
        file_name = os.path.basename(root)
        if file_name.endswith('.json') and file_name.startswith('s_'):
            return [root]
        else:
            raise Exception('No selection file found!')
    #
    result = []
    for file_name in os.listdir(root):
        file_path = root + '/' + file_name
        if os.path.isdir(file_path):
            result.extend(find_selection_files(file_path))
            continue
        if file_name.endswith('.json') and file_name.startswith('s_'):
            result.append(file_path)
    if len(result) == 0:
        raise Exception('No selection file found!')
    return result


def load_configuration():
    """Load and return the configuration.

    :return:
    """
    default_path_to_output = os.path.dirname(os.path.realpath(__file__)) + \
                             '/../data/output/screening/'
    parser = argparse.ArgumentParser(
        description='Perform screening on local machine.')
    parser.add_argument('-i', type=str, dest='input', required=True,
                        help='Relative path from script file to selections.')
    parser.add_argument('-m', type=str, dest='methods', required=True,
                        help='Comma separated names of methods files.')
    return vars(parser.parse_args());


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    #
    config = load_configuration()
    # Prepare list of methods.
    methods = []
    method_directory = os.path.dirname(os.path.realpath(__file__)) + '/methods/'
    for x in config['methods'].split(','):
        if x.endswith('*'):
            # Scan for files.
            for file_name in os.listdir(method_directory + x[0:-1]):
                if file_name.lower().endswith('.py'):
                    methods.append(x[0:-1] + file_name)
        else:
            methods.append(x + '.py')
    # Get selections and methods to screen.
    screening(methods, find_selection_files(config['input']))


if __name__ == '__main__':
    main()
