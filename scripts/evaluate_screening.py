#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Evaluate performance of screening.

Usage example:
    python evaluate_screening.py
"""

import os
import functools
import logging

from rdkit.ML.Scoring import Scoring

import eu_skodape_libs.core as core_libs

__license__ = "X11"


def main():
    core_libs.init_logging()
    for screen_entry in list_groups():
        eval_path = get_evaluation_path(screen_entry)
        if os.path.exists(eval_path):
            continue
        evaluate(screen_entry, eval_path)


def list_groups():
    screening_dir = get_data_directory() + "screening/"
    for method in os.listdir(screening_dir):
        for dataset in os.listdir(screening_dir + "/" + method):
            dataset_directory = screening_dir + "/" + method + "/" + dataset
            for selection in os.listdir(dataset_directory):
                for group in os.listdir(dataset_directory + "/" + selection):
                    group_directory = dataset_directory + "/" + \
                                      selection + "/" + group
                    yield {
                        "path": group_directory,
                        "dataset": dataset,
                        "selection": selection,
                        "group": group,
                        "method": method
                    }


def get_evaluation_path(screen_entry):
    return get_data_directory() + "/evaluation/" + \
           screen_entry["method"] + "/" + \
           screen_entry["dataset"] + "/" + \
           screen_entry["selection"] + "/" + \
           screen_entry["group"] + ".json.gz"


def get_data_directory():
    return os.path.dirname(os.path.realpath(__file__)) + \
           "/../data/"


def evaluate(screen_entry, output_path):
    logging.info("Screening %s/%s/%s",
                 screen_entry["dataset"],
                 screen_entry["selection"],
                 screen_entry["group"])
    core_libs.create_directory(output_path)
    ligand_names, decoys_names = load_definition(screen_entry)
    results = []
    for screen_output in list_screen_outputs(screen_entry):
        evaluation = evaluate_screening(
            ligand_names, decoys_names, screen_output)
        results.append(evaluation)
    core_libs.write_json(output_path, results)


def load_definition(screen_entry):
    definition_path = get_data_directory() + "/datasets/" + \
                      screen_entry["dataset"] + "/selections/" + \
                      screen_entry["selection"] + "/" + \
                      screen_entry["group"] + "/group.json"

    group = core_libs.read_json(definition_path)

    ligand_names = set()
    for item in group["data"]["ligands"]:
        ligand_names.add(item["name"])

    decoys_names = set()
    for item in group["data"]["decoys"]:
        decoys_names.add(item["name"])

    return ligand_names, decoys_names


def list_screen_outputs(screen_entry):
    screen_directory = get_data_directory() + "/screening/" + \
                       screen_entry["method"] + "/" + \
                       screen_entry["dataset"] + "/" + \
                       screen_entry["selection"] + "/" + \
                       screen_entry["group"] + "/"
    for file_name in os.listdir(screen_directory):
        file_path = screen_directory + file_name
        yield core_libs.read_json(file_path)


def evaluate_screening(ligand_names, decoys_names, screen):
    scores = screen["data"]["scores"]
    set_activity(ligand_names, decoys_names, scores)
    result = {}
    result.update(evaluate_best_case(scores))
    result.update(evaluate_worst_case(scores))
    return {
        "data": result,
        "metadata": screen["metadata"]
    }


def set_activity(ligand_names, decoys_names, scores):
    for score in scores:
        is_ligand = score["molecule"] in ligand_names
        is_decoy = score["molecule"] in decoys_names
        if is_ligand and not is_decoy:
            score["activity"] = 1
        elif not is_ligand and is_decoy:
            score["activity"] = 0
        else:
            raise Exception("Missing: '" + score["molecule"] + "'")


def evaluate_best_case(scores):
    scores_best = sorted(scores,
                         key=functools.cmp_to_key(score_compare_best),
                         reverse=True)
    result_object = {}
    for label in METRICS.keys():
        result_object["b-" + label] = METRICS[label](scores_best)
    return result_object


def evaluate_worst_case(scores):
    scores_worst = sorted(scores,
                          key=functools.cmp_to_key(score_compare_worst),
                          reverse=True)
    result_object = {}
    for label in METRICS.keys():
        result_object["w-" + label] = METRICS[label](scores_worst)
    return result_object


# region Score comparison

def score_compare_worst(left, right):
    if left["score"] == right["score"]:
        # Decide based on activity.
        if left["activity"] and right["activity"]:
            return 0
        elif left["activity"]:
            # Only left is active.
            return -1
        else:
            return 1
    else:
        if left["score"] > right["score"]:
            return 1
        else:
            return -1


def score_compare_best(left, right):
    if left["score"] == right["score"]:
        # Decide based on activity.
        if left["activity"] and right["activity"]:
            return 0
        elif left["activity"]:
            # Only left is active.
            return 1
        else:
            return -1
    else:
        if left["score"] > right["score"]:
            return 1
        else:
            return -1


# endregion

# region Metrics

METRICS = {
    "auc": lambda s: Scoring.CalcAUC(s, "activity"),
    "bedroc_20": lambda s: Scoring.CalcBEDROC(s, 'activity', 20)
}

# endregion


if __name__ == "__main__":
    main()
