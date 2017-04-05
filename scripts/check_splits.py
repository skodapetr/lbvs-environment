#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Check splits.
"""

import argparse
import logging
import json

from libs import data

__license__ = "X11"


def main():
    init_logging()
    args = parse_args()
    for reference in list_datasets(args["input"], args["collection"]):
        check(reference)


def init_logging():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] - %(message)s",
        datefmt="%H:%M:%S")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", required=False,
                        help="Dataset name in format: dataset/selection/group")
    parser.add_argument("-c", dest="collection", required=False,
                        help="Name of a collection.")
    return vars(parser.parse_args())


def list_datasets(input_filters, collection_name):
    if input_filters is not None:
        filters = input_filters["input"].split(",")
        return data.resolve(filters[0], filters[1], filters[2])
    if collection_name is not None:
        datasets = []
        collections = data.list_datasets_for_collection(collection_name)
        for collection_group in collections.values():
            for item in collection_group:
                datasets.extend(data.resolve(item[0], item[1], item[2]))
        return datasets
    logging.error("Missing -i or -c option.")
    exit()


def check(reference):
    # dataset, selection, group
    group = load_group(reference)

    print(group)

    exit()


def load_group(reference):
    path = data.dataset_to_path(reference) + "/group.json"
    with open(path, "r") as stream:
        return json.load(stream)


if __name__ == "__main__":
    main()
