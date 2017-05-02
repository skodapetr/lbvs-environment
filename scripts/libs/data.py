#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Provide ways to loading molecules for the test instances.
"""

import os
import re

__license__ = 'X11'

__DATA_DIRECTORY = os.path.dirname(os.path.realpath(__file__)) + '/../../data/'
__DATASET_DIRECTORY = __DATA_DIRECTORY + 'datasets/'
__MOLECULES_FILE_CACHE = {}


class DatasetReference(object):
    """Definition of a dataset reference.
    """

    def __init__(self, dataset, selection, group):
        self.dataset = dataset
        self.selection = selection
        self.group = group


class Molecules(object):
    """Object with molecules for a single test instance.
    """

    def __init__(self):
        self.test = []
        self.train = {
            "actives": [],
            "inactives": []
        }
        self.validation = {
            "actives": [],
            "inactives": []
        }


def list_datasets(as_path=False):
    """

    :param as_path:
    :return: Datasets in the platform.
    """
    datasets = [name for name in os.listdir(__DATASET_DIRECTORY)
                if os.path.isdir(__DATASET_DIRECTORY + name)]
    if as_path:
        return [__DATASET_DIRECTORY + name for name in datasets]
    else:
        return datasets


def list_selections(dataset, as_path=False):
    """

    :param dataset: Name of the dataset.
    :param as_path:
    :return: Selections in the dataset.
    """
    directory = __DATASET_DIRECTORY + dataset + '/selections/'
    selections = [name for name in os.listdir(directory)
                  if os.path.isdir(directory + name)]
    if as_path:
        return [directory + name for name in selections]
    else:
        return selections


def list_groups(dataset, selection, as_path=False):
    """

    :param dataset: Name of the dataset.
    :param selection: Name of the selection.
    :param as_path:
    :return: Groups in given selection and datasets.
    """
    directory = __DATASET_DIRECTORY + dataset + '/selections/' + selection + '/'
    if as_path:
        return [directory + name for name in os.listdir(directory)]
    else:
        return os.listdir(directory)


def list_instances_from_reference(
        dataset_reference, as_path=False):
    return list_instances(dataset_reference.dataset,
                          dataset_reference.selection,
                          dataset_reference.group, as_path)


def list_instances(dataset, selection, group, as_path=False):
    """

    :param dataset: Name of the dataset.
    :param selection: Name of the selection.
    :param group: Name of the group.
    :param as_path:
    :return: Instances for given dataset, selection and group.
    """
    directory = __DATASET_DIRECTORY + dataset + '/selections/' + \
                selection + '/' + group + '/'
    instances_names = [name for name in os.listdir(directory)
                       if name.startswith("s_")]
    if as_path:
        return [directory + name for name in instances_names]
    else:
        return instances_names


def __load_molecules(path):
    """

    :param path:
    :return: Valid molecules from given file.
    """
    global __MOLECULES_FILE_CACHE
    if path in __MOLECULES_FILE_CACHE:
        return __MOLECULES_FILE_CACHE[path]
    if len(__MOLECULES_FILE_CACHE) > 2:
        __MOLECULES_FILE_CACHE = {}
    import rdkit
    from rdkit import Chem
    molecules = [molecule for molecule in rdkit.Chem.SDMolSupplier(str(path))
                 if molecule is not None]
    __MOLECULES_FILE_CACHE[path] = molecules
    return molecules


def load_molecules(dataset_reference, instance_data):
    """

    :param dataset_reference:
    :param instance_data: Data of the instance.
    :return:
    """
    sdf_directory = __DATASET_DIRECTORY + dataset_reference.dataset + \
                    '/molecules/sdf/'
    molecules = {}
    for file in instance_data['data']['files']:
        sdf_path = sdf_directory + file + '.sdf'
        for molecule in __load_molecules(sdf_path):
            molecules[molecule.GetProp('_Name')] = molecule

    result = Molecules()
    for item in instance_data['data']['test']:
        result.test.append(molecules[item['name']])
    for item in instance_data['data']['train']['decoys']:
        result.train["inactives"].append(molecules[item['name']])
    for item in instance_data['data']['train']['ligands']:
        result.train["actives"].append(molecules[item['name']])
    if "validation" in instance_data['data']:
        for item in instance_data['data']['validation']['decoys']:
            result.validation["inactives"].append(molecules[item['name']])
        for item in instance_data['data']['validation']['ligands']:
            result.validation["actives"].append(molecules[item['name']])

    return result


def resolve(dataset_filter='.*', selection_filter='.*', group_filter='.*'):
    """

    :param dataset_filter:
    :param selection_filter:
    :param group_filter:
    :return: Array of matches to given filters.
    """
    result = []
    re_dataset = re.compile(dataset_filter)
    re_selection = re.compile(selection_filter)
    re_group = re.compile(group_filter)
    for dataset in list_datasets():
        if not re_dataset.match(dataset):
            continue
        for selection in list_selections(dataset):
            if not re_selection.match(selection):
                continue
            for group in list_groups(dataset, selection):
                if not re_group.match(group):
                    continue
                result.append(DatasetReference(dataset, selection, group))
    return result


def dataset_to_path(dataset_reference):
    """

    :param dataset_reference:
    :return: Path to group directory.
    """
    return __DATASET_DIRECTORY + dataset_reference.dataset + '/selections/' + \
           dataset_reference.selection + '/' + dataset_reference.group


def list_collections(as_path=False):
    """

    :param as_path:
    :return: List of collections.
    """
    if as_path:
        return [__DATA_DIRECTORY + name
                for name in os.listdir(__DATA_DIRECTORY + '/collections/')]
    else:
        return os.listdir(__DATA_DIRECTORY + '/collections/')


def list_datasets_for_collection(collection, default_selection=None):
    """

    :param collection:
    :param default_selection:
    :return: Groups of datasets.
    """
    collection_dir = __DATA_DIRECTORY + '/collections/' + collection + '/'
    result = {}
    for name in os.listdir(collection_dir):
        datasets_in_collection = []
        result[name] = datasets_in_collection
        with open(collection_dir + name) as stream:
            for line in stream:
                line = line.rstrip().split(',')
                datasets_in_collection.append([line[0], line[1], line[2]])
    return result


if __name__ == '__main__':
    raise Exception('This module should be used only as a library!')
