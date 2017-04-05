#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
"""
Download data used in the platform.

Usage example:
    python download_data.py
"""

import os
import json
import zipfile
from urllib import request

__license__ = "X11"

BASE_URL = "http://siret.ms.mff.cuni.cz/skoda/vs-datasets/"


def main():
    info = download_basic_info()
    options = ["Download concrete dataset molecules",
               "Download concrete dataset selections",
               "Download collections",
               "Download all",
               "Exit"]
    while True:
        option = menu_select_items(options)[0]
        if option == 0:
            download_molecules(info)
        elif option == 1:
            download_selections(info)
        elif option == 2:
            download_collections()
        elif option == 3:
            download_all(info)
        elif option == 4:
            break


def download_basic_info():
    print("Downloading basic info ... ", end="")
    url = BASE_URL + "/info.json"
    info = json.loads(request.urlopen(url).read().decode("utf-8"))
    print("done")
    return info


def menu_select_items(items, single=True):
    if single:
        print("Select single value:")
    else:
        print("Select a multiples value (comma separated list):")
    counter = 0
    for name in items:
        counter += 1
        print(" ", counter, ") ", name, sep="")
    while True:
        print(" > ", end="")
        try:
            value = input()
            if single:
                value = [value]
            else:
                value = value.split(",")
            value = [int(item) - 1 for item in value]
            # Check that all indexes are in bounds.
            for index in value:
                if index < 0 or index >= len(items):
                    raise IndexError()
        except:
            print("Invalid input!")
            continue
        return value


def download_molecules(info):
    datasets = menu_select_datasets(info, False)
    for item in datasets:
        url = BASE_URL + item["ref"] + "/molecules/sdf.zip"
        path = get_dataset_root() + item["dir"] + "/molecules/"
        download_and_unpack(url, path)


def menu_select_datasets(info, single=True):
    data = []
    data_names = []
    for dataset in info["datasets"]:
        data.append(dataset)
        data_names.append(dataset["name"])
    selected = menu_select_items(data_names, single)
    return [data[index] for index in selected]


def get_dataset_root():
    return os.path.dirname(os.path.realpath(__file__)) + "/../data/datasets/"


def download_selections(info):
    datasets = menu_select_datasets(info, False)
    for item in datasets:
        path = get_dataset_root() + item["dir"] + "/selections/"
        for selection in item["selections"]:
            url = BASE_URL + item["ref"] + "/selections/" + \
                  selection + ".zip"
            download_and_unpack(url, path)


def download_collections():
    url = BASE_URL + "/collections.zip"
    path = os.path.dirname(os.path.realpath(__file__)) + "/../data/"
    print(url, path)
    download_and_unpack(url, path)


def download_all(info):
    for item in info["datasets"]:
        # Download molecules.
        url = BASE_URL + item["ref"] + "/molecules/sdf.zip"
        path = get_dataset_root() + item["dir"] + "/molecules/"
        download_and_unpack(url, path)
        # Download selections.
        path = get_dataset_root() + item["dir"] + "/selections/"
        for selection in item["selections"]:
            url = BASE_URL + item["ref"] + "/selections/" + \
                  selection + ".zip"
            download_and_unpack(url, path)
    # Download collections
    download_collections()


def download_and_unpack(url, target_path):
    """Download file from given URI and unpack it to given directory.

    :param url:
    :param target_path:
    :return:
    """
    # Prepare output and temp directory.
    temp_file_path = os.path.dirname(os.path.realpath(__file__)) + \
                     "/../data/temp/download.tmp"
    create_directory(target_path)
    create_directory(os.path.dirname(temp_file_path))
    # Print info.
    print("URL:", url)
    print(" ->", target_path)
    # Download.
    print(" Downloading ", end="", flush=True)
    with request.urlopen(url) as response, \
            open(temp_file_path, "wb") as output_stream:
        total_size = int(response.headers["Content-Length"])
        # Download with progress bar.
        step_size = total_size / 20.0
        next_step = step_size
        position = 0
        for content in response:
            position += len(content)
            output_stream.write(content)
            if position > next_step:
                next_step += step_size
                print(".", end="", flush=True)
    print(" done\n Unzipping ... ", end="")
    with zipfile.ZipFile(temp_file_path) as zip_file:
        zip_file.extractall(target_path)
    os.remove(temp_file_path)
    print("done")


def create_directory(path):
    if os.path.exists(path):
        return
    os.makedirs(path)


if __name__ == "__main__":
    main()
