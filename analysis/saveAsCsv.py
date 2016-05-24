#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import logging
import argparse

import numpy

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

_root_path = os.path.dirname(os.path.realpath(__file__)) + '/../'


def scan(directory_path):
    """Scan and return all *.json files in directory recursively.

    :param directory_path:
    :return:
    """
    result = []
    for name in os.listdir(directory_path):
        path = directory_path + '/' + name
        if os.path.isdir(path):
            result.extend(scan(path))
            continue
        if os.path.isfile(path) and path.endswith('.json'):
            result.append(path)
    return result


def reorganize_runs_data(input_data):
    """Put runs on the same data together.

    For each method, and selection we get a list of properties from
    the executed runs.
    :param input_data:
    :return:
    """
    data = {}
    counter = 0
    for item in input_data:
        counter += 1
        if counter % 1000 == 0:
            logging.info('%d/%d', counter, len(input_data))
        #
        try:
            selection = item['method']['name'] + \
                        item['selection']['selection'] + \
                        item['selection']['dataset'] + \
                        item['selection']['target']
        except KeyError as ex:
            print('Missing key in:', item)
            raise ex

        if selection not in data:
            data[selection] = {
                'method': item['method'],
                'selection': item['selection'],
                # We use array to enable access from other properties.
                'count': [1],
                'BEDROC': [item['BEDROC']],
                'RIE': [item['RIE']],
                'AUC': [item['AUC']],
                'EF': [item['EF']],
                'metadata': [item['metadata']]
            }
            del data[selection]['selection']['index'];
        else:
            data[selection]['count'][0] += 1
            data[selection]['AUC'].append(item['AUC'])
            data[selection]['BEDROC'].append(item['BEDROC'])
            data[selection]['RIE'].append(item['RIE'])
            data[selection]['EF'].append(item['EF'])
            data[selection]['metadata'].append(item['metadata'])
    return data


# region Aggregation function definition

def select_maximum(data):
    return max(data)


def select_minimum(data):
    return min(data)


def select_mean(data):
    return numpy.mean(data)


aggregate_data = {
    'max': select_maximum,
    'min': select_minimum,
    'avg': select_mean
}


# endregion Aggregation function definition

def merge_runs(data, property):
    """Merge data from multiple runs.

    The output is data summary for each combination of selection and method.
    :param data:
    :param property: Information about property to pick.
    :return:
    """

    def select_property(object, property_name):
        """Select value of given name from given object.

        :param object:
        :param property_name:
        :return:
        """
        if property_name is None:
            return object
        #
        keys = [item.replace('\\', '') for item in property_name.split('[^\].')]
        value = object
        for key in keys:
            value = value[key]
        return value

    merged_data = {}
    for key in data:
        item = data[key]
        # Read values - first we use name to select object on read data
        # from. Next we use selector to fine tune the selection (EF, RIE, ..)
        values = [select_property(item, property['selector'])
                  for item in select_property(data[key], property['name'])
                  ]
        # Perform aggregation.
        value = aggregate_data[property['type']](values)
        #
        merged_data[key] = {
            'method': item['method'],
            'count': item['count'],
            'selection': item['selection'],
            'value': value
        }
    return merged_data

def read_configuration():
    """Read and return command line arguments.

    """
    parser = argparse.ArgumentParser('Execute given method on given data')
    parser.add_argument('-f', type=str, dest='input',
                        help='File with targets to screen.',
                        required=False)
    parser.add_argument('-o', type=str, dest='output',
                        help='Path to output CSV file.',
                        required=True)

    parser.add_argument('-n', type=str, dest='name',
                        help='Name of property.',
                        required=True)
    parser.add_argument('-s', type=str, dest='selector',
                        help='Selector in the property.',
                        required=False)
    parser.add_argument('-t', type=str, dest='type',
                        help='Aggregation type.',
                        required=True)

    return vars(parser.parse_args())


def _main():
    evaluation_root = _root_path + '/results/evaluation/'

    config = read_configuration()

    # config = {
    #     # 'file': './../collections/auc_8.0-8.5.dat', # 1
    #     'input': './../collections/auc_8.5-9.0.dat',  # OK
    #     # 'file': './../collections/auc_9.0-9.5.dat', # 4
    #     # 'file': './../collections/auc_9.8-1.0.dat', # 2
    #     # 'file': './../collections/all.dat', # 7
    #     'output': '../performance.csv'
    # }

    # Read files to process.
    if config['input'] is None:
        # Use all files.
        files = scan(evaluation_root)
    else:
        # Read filter.
        with open(config['input'], 'r') as stream:
            files_filter = [line.strip() for line in stream.readlines()]
        files = []
        for file_name in os.listdir(evaluation_root):
            for item in files_filter:
                file_path = evaluation_root + '/' + file_name + '/' + \
                            item + '/results.json'
                # Remove '/selections/' directory from the path.
                file_path = file_path.replace('/selections/', '/')
                # Test if the path exits.
                if os.path.exists(file_path):
                    files.append(file_path)

    # Load files with results.
    data = []
    logging.info('Loading files ...')
    for file in files:
        with open(file, 'r') as stream:
            data.extend(json.load(stream))
    logging.info('Loading files ... done')

    # Work with the data.
    logging.info('Merging data ...')
    data = reorganize_runs_data(data)
    data = merge_runs(data, config)
    logging.info('Merging data ... done')

    # Prepare data.
    def get(data_object, key):
        """Return object under given key, create new if non exists.

        :param data_object:
        :param key:
        :return:
        """
        if key not in data_object:
            data_object[key] = {}
        return data_object[key]

    methods = set()
    prepared_data = {}
    for item in data.values():
        method = item['method']['name']
        selection_dataset = item['selection']['dataset']
        selection_target = item['selection']['target']
        selection_selection = item['selection']['selection']
        # Store data.
        selection = get(get(get(prepared_data, selection_dataset),
                            selection_target), selection_selection)
        selection[method] = item['value']
        methods.add(method)

    # Prepare CSV
    header = ['""', '""', '""']
    lines = {}
    for dataset in sorted(prepared_data.keys()):
        for molecules in sorted(prepared_data[dataset].keys()):
            for selection in sorted(prepared_data[dataset][molecules].keys()):
                selections = prepared_data[dataset][molecules][selection]
                header[0] += ',"' + dataset + '"'
                header[1] += ',"' + molecules + '"'
                header[2] += ',"' + selection + '"'
                for method in methods:
                    if method not in lines:
                        lines[method] = '"' + method + '"'
                    lines[method] += ',"'
                    if method in selections:
                        lines[method] += str(selections[method])
                    else:
                        lines[method] += ''
                    lines[method] += '"'

    # Write to file.
    with open(config['output'], 'w') as stream:
        for line in header:
            stream.write(line)
            stream.write('\n')
        for key in sorted(lines.keys()):
            line = lines[key]
            stream.write(line)
            stream.write('\n')


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    #
    _main()
