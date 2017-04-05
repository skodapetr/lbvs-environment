#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Use given similarity method to produce a similarity matrix.

By default compute similarities between molecules of given files. If given
files A B, computes distances AxB. If more then two files are given program
fail.

If the --pair option is given compute distances between all molecules, ie.
(A,B,C)x(A,B,C). This can take long for bigger number of molecules.

Usage example:
    python compute_similarities.py -i A.sdf B.sdf -m tt/tt_tanimoto
"""

import imp
import sys
import numpy
import argparse
import os
import csv

import eu_skodape_libs.rdkit.core as rdkit_utils

__license__ = "X11"


def main():
    args = parse_args()
    module = load_method_module(args["method"])
    if args["pair"]:
        molecules = load_molecules(args["input"])
        similarity_matrix = compute_pair_similarity(molecules, module)
        print_square_matrix(molecules, similarity_matrix, args["name"])
    else:
        if not len(args["input"]) == 2:
            print("Expected number of files 2 given ", len(args["input"]))
            print("Two files must be given if --pair option is not used.")
            exit()
        left_molecules = load_molecules([args["input"][0]])
        right_molecules = load_molecules([args["input"][1]])
        similarity_matrix = compute_between_group_similarity(
            left_molecules, right_molecules, module)
        print_matrix(left_molecules, right_molecules,
                     similarity_matrix, args["name"])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="input", nargs='+',
                        help="File or a directory with molecule files.")
    parser.add_argument("-m", dest="method", required=True,
                        help="Method name.")
    parser.add_argument("--name", dest="name", default="_Name",
                        help="Used molecule name")
    parser.add_argument("--pair", dest="pair", action="store_true",
                        help="Compute pairwise similarities.")
    return vars(parser.parse_args())


def load_method_module(method_name):
    method_directory = os.path.dirname(os.path.realpath(__file__)) + \
                       "/../../methods/"
    method_module_path = method_directory + method_name + ".py"
    return imp.load_source(method_name, method_module_path)


def load_molecules(paths):
    molecules = []
    for path in paths:
        if not os.path.exists(path):
            raise Exception("Path '" + path + "' does not exists.")
        if os.path.isfile(path):
            molecules.extend(rdkit_utils.load_molecules(path))
        else:
            paths_to_load = [path + "/" + file_name
                             for file_name in os.listdir(path)]
            molecules.extend(rdkit_utils.load_molecules(paths_to_load))
    return molecules


def get_header(molecules, molecule_name_property):
    return [""] + [mol.GetProp(molecule_name_property) for mol in molecules]


def compute_pair_similarity(molecules, module):
    return compute_between_group_similarity(molecules, molecules, module)


def compute_between_group_similarity(molecules_x, molecules_y, module):
    matrix = numpy.zeros((len(molecules_x), len(molecules_y)))
    for y in range(0, len(molecules_y)):
        left = molecules_y[y]
        for x in range(0, len(molecules_x)):
            right = molecules_x[x]
            similarity = module.compute_similarity(left, right)
            matrix[x, y] = similarity
    return matrix


def print_square_matrix(molecules, similarity_matrix, molecule_name_property):
    writer = csv.writer(sys.stdout, delimiter=",", quotechar="\"",
                        quoting=csv.QUOTE_ALL, lineterminator="\n")
    writer.writerow(get_header(molecules, molecule_name_property))
    for y in range(0, len(molecules)):
        left = molecules[y]
        line = [left.GetProp(molecule_name_property)]
        for x in range(0, len(molecules)):
            line.append(similarity_matrix[y, x])
        writer.writerow(line)


def print_matrix(molecules_x, molecules_y, similarity_matrix,
                 molecule_name_property):
    writer = csv.writer(sys.stdout, delimiter=",", quotechar="\"",
                        quoting=csv.QUOTE_ALL, lineterminator="\n")
    writer.writerow(get_header(molecules_x, molecule_name_property))
    for y in range(0, len(molecules_y)):
        left = molecules_y[y]
        line = [left.GetProp(molecule_name_property)]
        for x in range(0, len(molecules_x)):
            line.append(similarity_matrix[x, y])
        writer.writerow(line)


if __name__ == "__main__":
    main()
