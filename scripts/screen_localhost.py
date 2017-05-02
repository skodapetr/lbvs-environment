#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Perform virtual screening with given methods on given data.

Usage example:
    screen_localhost.py -c 02 -m tt/tt_tanimoto
"""

import imp
import os
import logging
import argparse

from libs import data
from libs import core as core_libs

__license__ = "X11"


def main():
    core_libs.init_logging()
    args = parse_args()
    dataset_references = get_datasets_to_screen(args)
    methods_list = args["methods"].split(",")
    screen(methods_list, dataset_references, args["overwrite"])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="collection",
                        help="Collection name (directories in data/collection)")
    parser.add_argument("-d", dest="dataset",
                        help="Dataset name in format: dataset/selection/group")
    parser.add_argument("-m", dest="methods", required=True,
                        help="Comma separated list of methods: ap/ap_tanimoto")
    parser.add_argument("--overwrite", dest="overwrite", action="store_true",
                        help="Overwrite existing output.")
    return vars(parser.parse_args())


def get_datasets_to_screen(args):
    if args["collection"] is not None:
        return load_collection(args["collection"])
    elif args["dataset"] is not None:
        dataset = args["dataset"].split("/")
        return data.resolve(dataset[0], dataset[1], dataset[2])
    else:
        print("Either dataset [-d] of collection [-c] must be provided.")
        exit()


def load_collection(collection_name):
    dataset_references = []
    collections_map = data.list_datasets_for_collection(collection_name)
    for name in collections_map.keys():
        for dataset in collections_map[name]:
            dataset_references.extend(data.resolve(
                dataset[0], dataset[1], dataset[2]))
    return dataset_references


def screen(method_names, dataset_references, overwrite_existing):
    logging.info("Loading methods ...")
    method_modules = [load_method_module(name) for name in method_names]
    logging.info("Screening ...")
    for method in method_modules:
        for reference in dataset_references:
            screen_dataset(method, reference, overwrite_existing, True)
    logging.info("Screening ... done")


def load_method_module(method_name):
    method_directory = os.path.dirname(os.path.realpath(__file__)) + \
                       "/../methods/"
    method_module_path = method_directory + method_name + ".py"
    return imp.load_source(method_name, method_module_path)


def screen_dataset(method, dataset_reference, update_existing,
                   include_additional_info):
    """
    :param method:
    :param dataset_reference:
    :param update_existing: If true overwrite existing results.
    :param include_additional_info: True to include model and score information.
    :return:
    """
    input_directory = data.dataset_to_path(dataset_reference)
    for instance_file in os.listdir(input_directory):
        if not instance_file.startswith("s_"):
            continue
        instance_path = input_directory + "/" + instance_file
        screen_split(method, dataset_reference, update_existing,
                     include_additional_info, instance_path)


def screen_split(method, dataset_reference, update_existing,
                 include_additional_info, instance_path):
    output_root_directory = os.path.dirname(os.path.realpath(__file__)) + \
                            "/../data/screening/"
    instance_data = core_libs.read_json(instance_path)

    method_id = method.METADATA["id"]
    instance_id = create_selection_id(instance_data["metadata"])
    output_path = output_root_directory + method_id + "/" + \
                  instance_id + ".json.gz"

    if os.path.exists(output_path) and not update_existing:
        return

    molecules = data.load_molecules(dataset_reference, instance_data)
    execute_screening(method, molecules, output_path, instance_data,
                      include_additional_info)


def create_selection_id(metadata):
    return metadata["dataset"] + "/" + metadata["selection"] + "/" + \
           metadata["group"] + "/" + metadata["instance"]


def execute_screening(method, molecules, output_path,
                      instance_data, include_additional_info):
    logging.info("Screening %s/%s/%s/%s",
                 instance_data["metadata"]["dataset"],
                 instance_data["metadata"]["selection"],
                 instance_data["metadata"]["group"],
                 instance_data["metadata"]["instance"])
    model, model_info = method.create_model(
        molecules.train["actives"], molecules.train["inactives"])
    scores = []
    for molecule in molecules.test:
        computed_score = method.compute_score(model, molecule)
        score_object = {
            "molecule": molecule.GetProp("_Name"),
            "score": computed_score["value"]
        }
        if include_additional_info and "info" in computed_score:
            score_object["info"] = computed_score["info"]
        scores.append(score_object)
    results = {
        "metadata": create_screening_metadata(instance_data, method),
        "data": {
            "scores": scores
        }
    }
    if include_additional_info:
        results["data"]["model"] = model_info

    core_libs.write_json(output_path, results)


def create_screening_metadata(instance_data, method):
    return {
        "dataset": instance_data["metadata"]["dataset"],
        "instance": instance_data["metadata"]["instance"],
        "selection": instance_data["metadata"]["selection"],
        "group": instance_data["metadata"]["group"],
        "method": method.METADATA["id"]
    }


if __name__ == "__main__":
    main()
