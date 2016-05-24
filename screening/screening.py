#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Perform simple screening of selected dataset with selected method.

Expect datasets in the datasets directory, method in the methods
directory and store results to results/output.
"""

import argparse
import logging
import os
import re
import json

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

_root_path = os.path.dirname(os.path.realpath(__file__)) + '/../'


def create_output_path(root, selection, method):
    """Create and return output path for given selection and configuration.

    :param rppt:
    :param selection:
    :return:
    """
    return root + '/results/output/' + \
           method + '/' + \
           selection['info']['dataset'] + '/' + \
           selection['info']['selection'] + '/' + \
           selection['info']['target'] + '_' + \
           selection['info']['index'] + '.json'


def list_selection(root, dataset_filter, selection_filter, target_filter):
    """Return paths to molecules directories based on given parameters.

    :param root:
    :param dataset_filter:
    :param selection_filter:
    :param target_filter:
    :return:
    """
    dataset_filter = re.compile(dataset_filter)
    datasets = [dataset for dataset in os.listdir(root + '/datasets/')
                if dataset_filter.match(dataset)]
    #
    selections = []
    selection_filter = re.compile(selection_filter)
    target_filter = re.compile(target_filter)
    for dataset in datasets:
        selections_root = root + '/datasets/' + dataset + '/selections'
        for selection in os.listdir(selections_root):
            if not selection_filter.match(selection):
                continue
            for target in os.listdir(selections_root + '/' + selection):
                if target_filter.match(target):
                    selections.append(root + '/datasets/' +
                                      dataset + '/selections/' +
                                      selection + '/' +
                                      target)
    return selections


def cmd_screening(root, methods, selections):
    """Call screening script on every execution.

    :param root:
    :param methods:
    :param selections:
    :return:
    """
    import subprocess
    for item in selections:
        for file_name in os.listdir(item):
            selection_path = item + '/' + file_name
            with open(selection_path, 'r') as stream:
                selection = json.load(stream)
            # Execute screening method.
            for method in methods:
                output_path = create_output_path(root, selection, method[:-3])
                method_path = root + '/screening/methods/' + method
                thread = subprocess.Popen(
                    ['python',
                     method_path,
                     '-i', selection_path,
                     '-d', root + '/datasets/',
                     '-o', output_path],
                    shell=True)
                thread.wait()


def import_screening(root, methods, selections):
    """Perform screening be importing the Python script.

    :param root:
    :param methods:
    :param selections:
    :return:
    """
    import imp
    import rdkit
    # Import methods.
    methods_modules = {}
    for method in methods:
        print('method:', root + '/screening/methods/' + method)
        module = imp.load_source(method, root +
                                 '/screening/methods/' + method)
        methods_modules[method] = module
    # Use caching for working with molecules.
    cache = {
        'molecules': {},
        'files': {}
    }
    for item in selections:
        for file_name in os.listdir(item):
            with open(item + '/' + file_name, 'r') as stream:
                try:
                    selection = json.load(stream)
                except Exception as ex:
                    print('Error on file: ', item + '/' + file_name)
                    raise ex
            # Check if we need to screen this.
            missing_execution = False
            for method in methods:
                output_path = create_output_path(root, selection,
                                                 method[:-3])
                if not os.path.exists(output_path):
                    missing_execution = True
            if not missing_execution:
                continue
            # Load molecules.
            source_files = [root + '/datasets/' +
                            selection['info']['dataset'] + '/molecules/sdf/' +
                            item + '.sdf'
                            for item in selection['files']]
            if not cache['files'] == source_files:
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
            # Perform screening.
            for method in methods:
                module = methods_modules[method]
                output_path = create_output_path(root, selection,
                                                 method[:-3])
                lock_path = output_path.replace('.json', '-lock')

                if os.path.exists(output_path):
                    continue

                # Use directory based lock.
                try:
                    os.makedirs(lock_path)
                except Exception as ex:
                    # Other thread is computing this representation.
                    continue

                logging.info('Screening: %s %s \n -> %s', method,
                             selection['info']['dataset'] + '/' +
                             selection['info']['selection'] + '/' +
                             selection['info']['target'],
                             output_path)

                results = module.screening(selection, cache['molecules'])
                # Save results
                if not os.path.exists(os.path.dirname(output_path)):
                    os.makedirs(os.path.dirname(output_path))

                with open(output_path, 'w') as stream:
                    json.dump(results, stream)
                    # TODO Report screening start/end on the target.

                # Release lock.
                os.rmdir(lock_path)


def list_methods(root, method_filter):
    """Return name of methods.

    :param root:
    :param method_filter:
    :return:
    """
    methods = []
    method_filter = re.compile(method_filter)
    for method in os.listdir(root + '/screening/methods'):
        if method.endswith('.py') and method_filter.match(method):
            methods.append(method)
    return methods


def read_configuration():
    """Read and return command line arguments.

    """
    parser = argparse.ArgumentParser('Execute given method on given data')
    parser.add_argument('-r', type=str, dest='root',
                        help='Root of the LBVS-Environment directory.',
                        required=False)
    parser.add_argument('-f', type=str, dest='input',
                        help='File with targets to screen',
                        required=False)
    parser.add_argument('-m', type=str, dest='method_filter',
                        help='Regular expression filter for methods directory,'
                             'without the extension.',
                        required=False, default='.*')
    parser.add_argument('-M', type=str, dest='methods',
                        help='Comma separated names of methods to test. '
                             'If used the -m parameter is ignored.',
                        required=False)
    # Multiple levels for data input.
    parser.add_argument('-d', type=str, dest='dataset_filter',
                        help='Regular expression filter for dataset to use.',
                        required=False, default='.*')
    parser.add_argument('-t', type=str, dest='target_filter',
                        help='Regular expression filter for selections to use.',
                        required=False, default='.*')
    parser.add_argument('-s', type=str, dest='selection_filter',
                        help='Regular expression filter for selections to use.',
                        required=False, default='.*')
    parser.add_argument('--type', type=str, dest='type',
                        help='Execution type "cmd" or "import"',
                        required=False, default='import')
    parser.add_argument('--interactive', dest='interactive',
                        action='store_true',
                        help='It true load files to screen from input. '
                             'Filters are still applied.')
    return vars(parser.parse_args())


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')

    config = read_configuration()

    root = config.get('root', None)
    if root is None:
        root = _root_path

    # Read files to screen.
    if config['input'] is not None:
        with open(config['input'], 'r') as stream:
            selections = [root + '/datasets/' + line.strip() for line in stream]
    else:
        selections = list_selection(root,
                                    config.get('dataset_filter', '.*'),
                                    config.get('selection_filter', '.*'),
                                    config.get('target_filter', '.*'))

    if config['methods'] is not None:
        methods = [x + '.py' for x in config['methods'].split(',')]
    else:
        methods = list_methods(root, config.get('method_filter', '.*'))

    if config['type'] == 'cmd':
        cmd_screening(root, methods, selections)
    else:
        import_screening(root, methods, selections)


if __name__ == '__main__':
    main()
