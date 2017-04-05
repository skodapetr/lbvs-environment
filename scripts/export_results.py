#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Print results as a CSV file.

Usage example:
    python export_results.py
    python export_results.py --columns
"""

import os
import json
import gzip
import csv
import sys
import argparse

__license__ = "X11"


def main():
    args = parse_args()
    writer = csv.writer(sys.stdout, delimiter=",", quotechar="\"",
                        quoting=csv.QUOTE_ALL, lineterminator="\n")
    files = recursive_find_json_gz_files(get_evaluation_path())
    if args["columns"]:
        row_writer = MethodColumnWriter()
        row_writer.write(files, writer)
    else:
        row_writer = MethodRowWriter()
        row_writer.write(files, writer)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--columns", dest="columns", action="store_true",
                        help="Print methods as columns")
    return vars(parser.parse_args())


def recursive_find_json_gz_files(directory):
    files = []
    for file_name in os.listdir(directory):
        path = directory + "/" + file_name
        if os.path.isdir(path):
            files.extend(recursive_find_json_gz_files(path))
        if file_name.endswith(".json.gz"):
            files.append(path)
    return files


def get_evaluation_path():
    return os.path.dirname(os.path.realpath(__file__)) + \
           "/../data/evaluation"


class MethodRowWriter:
    def write_header(self):
        header = ["method", "dataset", "selection", "group", "instance",
                  "auc", "bedroc_20"]
        self.writer.writerow(header)

    def evaluation_to_line(self, evaluation):
        data = evaluation["data"]
        metadata = evaluation["metadata"]
        return [
            metadata["method"],
            metadata["dataset"],
            metadata["selection"],
            metadata["group"],
            metadata["instance"],
            data["w-auc"],
            data["w-bedroc_20"]
        ]

    def write_line(self, path):
        with gzip.open(path, "rt") as stream:
            data = json.load(stream)
        for evaluation in data:
            line = self.evaluation_to_line(evaluation)
            self.writer.writerow(line)

    def write(self, files, writer):
        self.writer = writer
        self.write_header()
        for path in files:
            self.write_line(path)


class MethodColumnWriter:
    def extract_method_names(self):
        names = set()
        for item in self.data:
            names.add(item[1])
        self.methods = list(names)

    def write_header(self):
        header = ["dataset", "selection", "group", "instance"]
        for method in self.methods:
            header.append(method + "-auc")
            header.append(method + "-bedroc_20")
        self.writer.writerow(header)

    def evaluation_to_record(self, evaluation):
        data = evaluation["data"]
        metadata = evaluation["metadata"]
        id = metadata["dataset"] + metadata["selection"] + \
             metadata["group"] + metadata["instance"]
        return [
            id,
            metadata["method"],
            metadata["dataset"],
            metadata["selection"],
            metadata["group"],
            metadata["instance"],
            data["w-auc"],
            data["w-bedroc_20"]
        ]

    def load_file(self, path):
        with gzip.open(path, "rt") as stream:
            data = json.load(stream)
        for evaluation in data:
            self.data.append(self.evaluation_to_record(evaluation))

    def load_files(self, files):
        self.data = []
        for file in files:
            self.load_file(file)

    def create_line(self, data_item):
        return [data_item[2], data_item[3], data_item[4], data_item[5]] + \
               ([0] * (2 * len(self.methods)))

    def prepare_lines(self):
        lines = {}
        for item in self.data:
            id = item[0]
            if id not in lines:
                lines[id] = self.create_line(item)
            line = lines[id]
            index = self.methods.index(item[1])
            # 4 is offset (dataset, selection, group, instance)
            line[4 + 2 * index] = item[6]
            line[4 + (2 * index) + 1] = item[7]
        self.lines = sorted(
            lines.values(), key=lambda m: m[0] + m[1] + m[2] + m[3])

    def print_lines(self):
        for line in self.lines:
            self.writer.writerow(line)

    def write(self, files, writer):
        self.writer = writer
        self.load_files(files)
        self.extract_method_names()
        self.write_header()
        self.prepare_lines()
        self.print_lines()


if __name__ == "__main__":
    main()
