#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract SDF files with train/test actives/decoys from the platform.

Usage example:
    python extract_sdf.py -i {GROUP_TO_EXTRACT} -o {OUTPUT_DIR}
    python extract_sdf.py -c {COLLECTION_TO_EXTRACT} -o {OUTPUT_DIR}
"""

import os
import argparse
import logging
import gzip
import json

from libs import data

__license__ = "X11"


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%H:%M:%S')

    args = parse_args()
    output_directory = args["output"]
    export_options = create_export_options(args)

    for reference in list_datasets(args["input"], args["collection"]):
        extract_group(reference, output_directory, export_options)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", required=False,
                        help="Dataset name in format: dataset/selection/group")
    parser.add_argument("-c", dest="collection", required=False,
                        help="Name of a collection.")
    parser.add_argument("-o", dest="output", required=True,
                        help="Output directory.")
    parser.add_argument("--flat", dest="flat", action="store_true",
                        help="Use flat output directory structure.")
    parser.add_argument("--gzip", dest="gzip", action="store_true",
                        help="Use gzip.")
    parser.add_argument("--info", dest="include_info", action="store_true",
                        help="Save JSON with dataset reference.")
    return vars(parser.parse_args())


def create_export_options(args):
    return {
        "flat": args["flat"],
        "gzip": args["gzip"],
        "info": args["include_info"]
    }


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


def extract_group(reference, output_directory, options):
    path = data.dataset_to_path(reference)
    actives, inactives = load_group_file(path)
    for file in os.listdir(path):
        if not file.startswith("s_"):
            continue
        instance_name = file[0:file.rfind(".")]
        extract_instance(reference, instance_name, actives, inactives,
                         output_directory, options)


def load_group_file(path_to_group):
    group_data = read_json(path_to_group + "/group.json")
    group_actives = set(
        [item["name"] for item in group_data["data"]["ligands"]])
    group_inactives = set(
        [item["name"] for item in group_data["data"]["decoys"]])
    return group_actives, group_inactives


def read_json(path):
    with open(path) as stream:
        return json.load(stream)


def extract_instance(reference, instance, actives, inactives,
                     output_directory, options):
    output_dir = output_path_directory(
        output_directory, reference, instance, options["flat"])
    logging.info("Extracting (%s,%s,%s,%s) ...", reference.dataset,
                 reference.selection, reference.group, instance)

    instance_file = data.dataset_to_path(reference) + \
                    "/" + instance + ".json"

    extract_instance_sdf(reference, instance_file,
                         output_dir, actives, inactives,
                         options["gzip"])

    if options["info"]:
        save_info(output_dir, reference, instance_file)


def output_path_directory(output, reference, instance, is_flat):
    if is_flat:
        return output + "/" + reference.dataset + "-" + reference.selection + \
               "-" + reference.group + "/" + instance + "/"
    else:
        return output + "/" + reference.dataset + "/" + reference.selection + \
               "/" + reference.group + "/" + instance + "/"


def extract_instance_sdf(reference, instance_file, output_path,
                         group_actives, group_inactives, gzip_output):
    split_data = read_json(instance_file)

    global MOLECULE_SOURCE
    molecules = MOLECULE_SOURCE.read_sdf_for_instance(reference, split_data)

    os.makedirs(output_path, exist_ok=True)
    path_prefix = output_path + "/" + reference.group

    train_actives_path = path_prefix + "_actives_train.sdf"
    with open_write_stream(train_actives_path, gzip_output) as stream:
        for molecule in split_data["data"]["train"]["ligands"]:
            stream.write(molecules[molecule["name"]])

    test_actives_path = path_prefix + "_actives_test.sdf"
    with open_write_stream(test_actives_path, gzip_output) as stream:
        for molecule in split_data["data"]["test"]:
            if molecule["name"] in group_actives:
                stream.write(molecules[molecule["name"]])

    if "validation" in split_data["data"]:
        validation_actives_path = path_prefix + "_actives_validation.sdf"
        with open_write_stream(validation_actives_path, gzip_output) as stream:
            for molecule in split_data["data"]["validation"]["ligands"]:
                stream.write(molecules[molecule["name"]])

    train_inactives_path = path_prefix + "_decoys_train.sdf"
    with open_write_stream(train_inactives_path, gzip_output) as stream:
        for molecule in split_data["data"]["train"]["decoys"]:
            stream.write(molecules[molecule["name"]])

    test_inactives_path = path_prefix + "_decoys_test.sdf"
    with open_write_stream(test_inactives_path, gzip_output) as stream:
        for molecule in split_data["data"]["test"]:
            if molecule["name"] in group_inactives:
                stream.write(molecules[molecule["name"]])

    if "validation" in split_data["data"]:
        validation_inactives_path = path_prefix + "_decoys_validation.sdf"
        with open_write_stream(validation_inactives_path,
                               gzip_output) as stream:
            for molecule in split_data["data"]["validation"]["decoys"]:
                stream.write(molecules[molecule["name"]])


def open_write_stream(path, gzip_output):
    if gzip_output:
        return gzip.open(path + ".gz")
    else:
        return open(path)


class MoleculeSource(object):
    def __init__(self):
        self.cache = {}

    def read_sdf_for_instance(self, reference, split_data):
        sdf_base_path = os.path.dirname(os.path.realpath(__file__)) + "/../" + \
                        "/data/datasets/" + reference.dataset + \
                        "/molecules/sdf/"

        files = [sdf_base_path + relative_path + ".sdf"
                 for relative_path in split_data["data"]["files"]]

        molecules = {}
        for file in files:
            molecules.update(self._read_molecule_file(file))

        return molecules

    def _get_molecule_for_file(self, path):
        if path in self.cache:
            return self.cache[path]
        if len(self.cache) > 2:
            self.cache = {}
        molecules = self._read_molecule_file(path)
        self.cache[path] = molecules
        return molecules

    def _read_molecule_file(self, path):
        molecules = {}
        mol = ""
        mol_id = ""
        with open(path) as stream:
            for line in stream:
                if mol == "":
                    mol_id = line.strip("\n")
                mol += line
                if line.startswith("$$$$"):
                    molecules[mol_id] = mol
                    mol = ""
        return molecules


MOLECULE_SOURCE = MoleculeSource()


def save_info(output_dir, reference, instance):
    write_json(output_dir + "/metadata.json", {
        "dataset": reference.dataset,
        "group": reference.group,
        "instance": instance,
        "selection": reference.selection
    })


def write_json(path, object):
    create_directory(os.path.dirname(path))
    with open(path, "w") as stream:
        json.dump(object, stream, indent=2)


def create_directory(path):
    if not os.path.exists(path) and not path == "":
        os.makedirs(path)


if __name__ == "__main__":
    main()
