#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create splits for given datasets.
"""

import os
import math
import json

import eu_skodape_libs.rdkit.core as rdkit_utils


def main():
    options, inputs = read_configuration()
    for item in inputs:
        create_splits(item["dataset"], item["group"], options,
                      generate_buckets_fixed_size_folding)


def read_configuration():
    options = {
        "train": {
            "actives": 1,
            "decoys": 1
        },
        "validation": {
            "actives": 2,
            "decoys": 10
        },
        "test": {
            "actives": 3,
            "decoys": 50
        },
        "selection": "cv_00_10_100_20_1000_30_5000",
        # Buckets are used to split molecules.
        "bucket": {
            "actives": 10,
            "decoys": 100
        }
    }
    inputs = [
        {
            "dataset": "10.1021%2Fci200412p",
            "group": "5HT2B_Antagonist"
        },
        {
            "dataset": "10.1021%2Fci200412p",
            "group": "5HT2C_Antagonist"
        },
        {
            "dataset": "10.1021%2Fci200412p",
            "group": "ADA2A_Antagonist"
        },
        {
            "dataset": "10.1021%2Fci2005274",
            "group": "CDK2"
        },
        {
            "dataset": "10.1021%2Fci5005515",
            "group": "HDAC01"
        },
        {
            "dataset": "10.1021%2Fjm500132p",
            "group": "PXR_Agonist"
        }
    ]
    return options, inputs


def create_splits(dataset, group, options, bucket_generator):
    actives = load_actives(dataset, group)
    decoys = load_decoys(dataset, group)
    index = 0
    used_actives = set()
    used_decoys = set()
    for split_data in generate_by_folding(
            actives, decoys, options, bucket_generator):
        index += 1
        split = {
            "metadata": {
                "dataset": dataset,
                "group": group,
                "instance": str(index).zfill(3),
                "selection": options["selection"]
            },
            "data": {
                "train": {
                    "ligands": names_to_refs(
                        split_data["train"]['ligands']),
                    "decoys": names_to_refs(
                        split_data["train"]['decoys'])
                },
                "validation": {
                    "ligands": names_to_refs(
                        split_data["validation"]['ligands']),
                    "decoys": names_to_refs(
                        split_data["validation"]['decoys'])
                },
                "test": names_to_refs(split_data["test"]["decoys"] +
                                      split_data["test"]["ligands"]),
                "files": [
                    group + "_Ligands",
                    group + "_Decoys"
                ]
            }
        }
        used_actives.update(split_data["test"]["ligands"])
        used_actives.update(split_data["validation"]["ligands"])
        used_actives.update(split_data["train"]["ligands"])
        used_decoys.update(split_data["test"]["decoys"])
        used_decoys.update(split_data["validation"]["decoys"])
        used_decoys.update(split_data["train"]["decoys"])
        #
        save_split_file(split)
    print("Split count: ", index)
    create_group_file(dataset, group, options,
                      list(used_actives),
                      list(used_decoys))
    if not len(used_decoys.intersection(used_actives)) == 0:
        print("Molecules: ", used_decoys.intersection(used_actives))
        raise Exception("Invalid data!")


def load_actives(dataset, group):
    path = get_actives_path(dataset, group)
    return rdkit_utils.load_molecules(path)


def names_to_refs(names):
    return [{"name": name} for name in names]


def get_actives_path(dataset, group):
    return get_dataset_root() + "datasets/" + dataset + "/molecules/sdf/" + \
           group + "_Ligands.sdf"


def load_decoys(dataset, group):
    path = get_decoys_path(dataset, group)
    return rdkit_utils.load_molecules(path)


def get_decoys_path(dataset, group):
    return get_dataset_root() + "datasets/" + dataset + "/molecules/sdf/" + \
           group + "_Decoys.sdf"


def get_dataset_root():
    return os.path.dirname(os.path.realpath(__file__)) + "/../data/"


def generate_by_folding(actives, decoys, options, bucket_generator):
    actives_bucket = split_to_buckets(actives, options["bucket"]["actives"])
    decoys_bucket = split_to_buckets(decoys, options["bucket"]["decoys"])

    for bucket_info in bucket_generator(
            len(actives_bucket), len(decoys_bucket), options):
        train_actives = molecules_to_name(merge_selected(
            actives_bucket, bucket_info["train"]["actives"]))
        train_decoys = molecules_to_name(merge_selected(
            decoys_bucket, bucket_info["train"]["decoys"]))

        validation_actives = molecules_to_name(merge_selected(
            actives_bucket, bucket_info["validation"]["actives"]))
        validation_decoys = molecules_to_name(merge_selected(
            decoys_bucket, bucket_info["validation"]["decoys"]))

        test_actives = molecules_to_name(merge_selected(
            actives_bucket, bucket_info["test"]["actives"]))
        test_decoys = molecules_to_name(merge_selected(
            decoys_bucket, bucket_info["test"]["decoys"]))

        yield {
            "test": {
                "decoys": train_decoys,
                "ligands": train_actives
            },
            "validation": {
                "decoys": validation_decoys,
                "ligands": validation_actives
            },
            "train": {
                "decoys": test_decoys,
                "ligands": test_actives
            }
        }


def split_to_buckets(molecules, bucket_size):
    bucket_count = int(math.floor(len(molecules) / bucket_size))
    buckets = []
    for bucket_index in range(0, bucket_count):
        begin = bucket_size * bucket_index
        end = bucket_size * (bucket_index + 1)
        bucket = []
        for index in range(begin, end):
            bucket.append(molecules[index])
        buckets.append(bucket)
    return buckets


def generate_buckets_fixed_size_folding(actives_len, decoy_len, options):
    train_a = options["train"]["actives"]
    train_d = options["train"]["decoys"]

    valid_a = options["validation"]["actives"]
    valid_d = options["validation"]["decoys"]

    test_a = options["test"]["actives"]
    test_d = options["test"]["decoys"]

    for selected_index in range(0, actives_len, train_a + valid_a):
        train_a_range = create_range(selected_index, train_a)
        valid_a_range = create_range(selected_index + train_a, valid_a)
        for decoys_index in range(0, decoy_len, train_d + valid_d):
            train_d_range = create_range(decoys_index, train_d)
            valid_d_range = create_range(decoys_index + train_d, valid_d)

            unused_actives = diff(range(0, actives_len),
                                  merge(train_a_range, valid_a_range))

            unused_decoys = diff(range(0, decoy_len),
                                 merge(train_d_range, valid_d_range))

            yield {
                "train": {
                    "actives": train_a_range,
                    "decoys": train_d_range
                },
                "validation": {
                    "actives": valid_a_range,
                    "decoys": valid_d_range
                },
                "test": {
                    "actives": unused_actives[0:test_a],
                    "decoys": unused_decoys[0:test_d]
                }
            }


def create_range(start, size):
    return list(range(start, start + size))


def diff(left, right):
    second = set(right)
    return [item for item in left if item not in second]


def merge(left, right):
    output = []
    output.extend(left)
    output.extend(right)
    return output


def merge_selected(buckets, to_merge):
    output = []
    for index in range(0, len(buckets)):
        if not index in to_merge:
            continue
        output.extend(buckets[index])
    return output


def molecules_to_name(molecules):
    return [mol.GetProp("_Name") for mol in molecules]


def save_split_file(split_data):
    path = get_split_path(split_data)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as stream:
        json.dump(split_data, stream, indent=2)


def get_split_path(split_data):
    dataset = split_data["metadata"]["dataset"]
    selection = split_data["metadata"]["selection"]
    group = split_data["metadata"]["group"]
    instance = split_data["metadata"]["instance"]
    return get_dataset_root() + "datasets/" + dataset + "/selections/" + \
           selection + "/" + group + "/s_" + instance + ".json"


def create_group_file(dataset, group, options, actives, decoys):
    path = get_group_path(dataset, group, options["selection"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as stream:
        json.dump({
            "metadata": {
                "dataset": dataset,
                "group": group,
                "selection": options["selection"]
            },
            "data": {
                "decoys": names_to_refs(decoys),
                "ligands": names_to_refs(actives)
            }
        }, stream, indent=2)


def get_group_path(dataset, group, selection):
    return get_dataset_root() + "datasets/" + dataset + "/selections/" + \
           selection + "/" + group + "/group.json"


if __name__ == "__main__":
    main()
