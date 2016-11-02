#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Create a summary from the evaluation.

Read data from the results in ../data/evaluation/results.csv. Performs
aggregation and print the output CSV file to standard output.
"""

import os
import logging
import csv
import argparse
import statistics
import sys

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

AGGREGATION_METHODS = {
    'mean': statistics.mean
}


def initialize_values(values, row):
    """Initialize values.

    :param values:
    :param row:
    :return:
    """
    for item in values:
        item['index'] = row.index(item['name'])
        if item['index'] == -1:
            raise Exception('Missing column ' + item['value'])


def load_configuration():
    """Load and return the configuration.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Perform screening on local machine.')
    parser.add_argument('-m', type=str, dest='methods', required=False,
                        help='Method filter, comma separated list.')
    parser.add_argument('-s', type=str, dest='selection', required=False,
                        help='Selection filter, comma separated list.')
    parser.add_argument('-c', type=str, dest='collection', required=False,
                        help='Path to collection file that is used as a filter')
    parser.add_argument('-v', type=str, dest='values', required=True,
                        help='Value name, comma separated list.')
    parser.add_argument('-a', type=str, dest='aggregation', required=True,
                        help='Aggregation functions for values.')
    return vars(parser.parse_args())


def main():
    config = load_configuration()
    # Prepare filters.
    method_filter = None
    if config['methods'] is not None:
        method_filter = config['methods'].split(',')
    selection_filter = None
    if config['selection'] is not None:
        selection_filter = config['selection'].split(',')
    collection_filter = None
    if config['collection'] is not None:
        collection_filter = []
        with open(config['collection']) as input_stream:
            for row in input_stream:
                filter = row.rstrip().replace('selections/', '')
                collection_filter.append(filter)
                print(collection_filter)
                break
    config['values'] = config['values'].split(',')
    config['aggregation'] = config['aggregation'].split(',')
    columns = []
    for index in range(0, len(config['values'])):
        columns.append({
            'name': config['values'][index],
            'method': AGGREGATION_METHODS[config['aggregation'][index]],
            'method_name': config['aggregation'][index]
        })
    # Load data.
    logging.info('Loading data ...')
    base_path = os.path.dirname(os.path.realpath(__file__)) + '/../'
    eval_path = base_path + 'data/evaluation/results.csv'
    records = {}
    with open(eval_path) as input_stream:
        csv_reader = csv.reader(input_stream, delimiter=',', quotechar='"')
        first_line = True
        for row in csv_reader:
            if first_line:
                initialize_values(columns, row)
                first_line = False
                continue
            path = row[1]
            dataset = row[1]
            selection = row[2]
            group = row[3]
            method = row[5]
            # Test filters.
            if method_filter is not None and method not in method_filter:
                continue
            if selection_filter is not None and selection not in selection_filter:
                continue
            if collection_filter is not None:
                collection_match = False
                for item in collection_filter:
                    if path in item:
                        collection_match = True
                        break
                if not collection_match:
                    continue
            # Store values.
            id = dataset + selection + method + group
            if id not in records:
                records[id] = {
                    'dataset': dataset,
                    'selection': selection,
                    'group': group,
                    'method': method,
                    'values': []
                }
            records[id]['values'].append(
                [float(row[v['index']]) for v in columns])
    logging.info('Loading data ... done')
    # Perform aggregation over data.
    import json
    rows = []
    # Header
    header = ['dataset', 'selection', 'group', 'method']
    for index in range(0, len(columns)):
        column = columns[index]
        header.append(column['name'] + '_' + column['method_name'])
    rows.append(header)
    # Data
    for record in records.values():
        new_row = [record['dataset'], record['selection'],
                   record['group'], record['method']]
        for index in range(0, len(columns)):
            column = columns[index]
            new_row.append(column['method'](
                [item[index] for item in record['values']]))
        rows.append(new_row)
    # Print
    csv_writer = csv.writer(sys.stdout, delimiter=',', quotechar='"',
                            lineterminator='\n')
    for row in rows:
        csv_writer.writerow(row)
    exit()


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%H:%M:%S')
    main()
