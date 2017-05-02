#!/usr/bin/env python
# -*- coding: utf-8 -*-


import json
import csv
import os
import logging
import gzip

__license__ = "X11"


def init_logging():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%H:%M:%S')


def create_directory(path):
    if not os.path.exists(path) and not path == "":
        os.makedirs(path)


def create_parent_directory(path):
    parent_directory = os.path.dirname(path)
    if not os.path.exists(parent_directory) and not parent_directory == "":
        os.makedirs(parent_directory)


def read_json(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as stream:
            return json.load(stream)
    else:
        with open(path, "r") as stream:
            return json.load(stream)


def write_json(path, object_to_write):
    create_parent_directory(path)
    if path.endswith(".gz"):
        with gzip.open(path, "wt") as stream:
            json.dump(object_to_write, stream, indent=2)
    else:
        with open(path, "w") as stream:
            json.dump(object_to_write, stream, indent=2)


def read_csv_as_object(path):
    """
    Read CSV lines as objects.
    """
    results = []
    with open(path) as stream:
        reader = csv.reader(stream, delimiter=",", quotechar='"')
        header = next(reader)
        for row in reader:
            new_object = {}
            for index in range(0, len(row)):
                new_object[header[index]] = row[index]
            results.append(new_object)
    return results


if __name__ == "__main__":
    raise Exception("This module can be used only as a library!")
