#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Let user download selected datasets.

"""

import os
import json
import zipfile
from urllib import request

__author__ = 'Petr Škoda'
__license__ = 'X11'
__email__ = 'skoda@ksi.mff.cuni.cz'

_base_uri = 'http://siret.ms.mff.cuni.cz/skoda/vs-datasets/'
_root_path = os.path.dirname(os.path.realpath(__file__)) + '/../'


def create_if_not_exists(path):
    """Create directory if does not exists.

    :param path:
    :return:
    """
    if not os.path.exists(path):
        os.makedirs(path)


def download_basic_info():
    """Download and return basic information about datasets.

    :return:
    """
    info = {}
    print('Downloading basic info ... ', end='')
    url = _base_uri + '/info.json'
    root_info = json.loads(request.urlopen(url).read().decode('utf-8'))
    for item in root_info['datasets']:
        url = _base_uri + item['ref'] + '/info.json'
        item_content = request.urlopen(url).read().decode('utf-8')
        info_path = _root_path + 'datasets/' + item['dir'] + '/info.json'
        create_if_not_exists(os.path.dirname(info_path))
        info[item['dir']] = {
            'ref': item['ref'],
            'dir': item['dir'],
            'value': json.loads(item_content)
        }
        with open(info_path, 'w') as stream:
            json.dump(info[item['dir']]['value'], stream, indent=1)
    print('done')
    return info


def select_items(items, single=True):
    """Let user select items and retutn their indexes.

    :param items:
    :param single:
    :return:
    """
    if single:
        print('Select single value:')
    else:
        print('Select a multiples value (comma separated list):')
    counter = 0
    for name in items:
        counter += 1
        print(' ', counter, ') ', name, sep='')
    while True:
        print(' >', end='')
        try:
            value = input()
            if single:
                value = [value]
            else:
                value = value.split(',')
            value = [int(item) - 1 for item in value]
            # Check that all indexes are in bounds.
            for index in value:
                if index < 0 or index >= len(items):
                    raise IndexError()
        except:
            print('Invalid input!')
            continue
        return value


def select_datasets(info, single=True):
    """Let user select one or more datasets.

    :param info:
    :param single:
    :return:
    """
    data = []
    data_names = []
    for key in sorted(info.keys()):
        item = info[key]
        name = item['value']['dataset']['name']
        if name == '':
            name = item['value']['dataset']['id']
        data.append(key)
        data_names.append(name)
    selected = select_items(data_names, single)
    return [data[index] for index in selected]


def download_and_unpack(url, target_path):
    """Download file from given URI and unpack it to given directory.

    :param url:
    :param target_path:
    :return:
    """
    # Prepare output and temp directory.
    create_if_not_exists(target_path)
    temp_file_path = _root_path + '/temp/download.tmp'
    create_if_not_exists(os.path.dirname(temp_file_path))
    # Print info.
    print('URL:', url)
    print(' ->', target_path)
    # Download.
    print(' Downloading ', end="", flush=True)
    with request.urlopen(url) as response, \
            open(temp_file_path, 'wb') as output_stream:
        total_size = int(response.headers['Content-Length'])
        # Download with progress bar.
        step_size = total_size / 20.0
        next_step = step_size
        position = 0
        for content in response:
            position += len(content)
            output_stream.write(content)
            if position > next_step:
                next_step += step_size
                print('.', end="", flush=True)
    print(' done\n Unzipping ... ', end='')
    with zipfile.ZipFile(temp_file_path) as zip_file:
        zip_file.extractall(target_path)
    os.remove(temp_file_path)
    print('done')


def download_molecules(info):
    """Download molecules for selected datasets.

    :param info:
    :return:
    """
    datasets = select_datasets(info, False)
    for key in datasets:
        item = info[key]
        url = _base_uri + item['ref'] + '/molecules/sdf.zip'
        path = _root_path + 'datasets/' + item['dir'] + '/molecules/'
        download_and_unpack(url, path)


def download_selections(info):
    """Download selections for selected datasets.

    :param info:
    :return:
    """
    datasets = select_datasets(info, False)
    for key in datasets:
        item = info[key]
        path = _root_path + 'datasets/' + item['dir'] + '/selections/'
        for selection in item['value']['selections']:
            url = _base_uri + item['ref'] + '/selections/' + \
                  selection['name'] + '.zip'
            download_and_unpack(url, path)


def download_all(info):
    """Download all datasets and selections.

    :param info:
    :return:
    """
    for item in info.values():
        # Download dataset.
        url = _base_uri + item['ref'] + '/molecules/sdf.zip'
        path = _root_path + 'datasets/' + item['dir'] + '/molecules/'
        download_and_unpack(url, path)
        # Download definitions.
        path = _root_path + 'datasets/' + item['dir'] + '/selections/'
        for selection in item['value']['selections']:
            url = _base_uri + item['ref'] + '/selections/' + \
                  selection['name'] + '.zip'
            download_and_unpack(url, path)


def main():
    info = download_basic_info()
    # Main menu.
    options = ['Download molecules for datasets',
               'Download selections for datasets',
               'Download all',
               'Exit']
    while True:
        option = select_items(options)[0]
        if option == 0:
            download_molecules(info)
        elif option == 1:
            download_selections(info)
        elif option == 2:
            download_all(info)
        elif option == 3:
            break


if __name__ == '__main__':
    main()
