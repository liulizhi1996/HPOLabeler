#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Implementation of Neighbor method.

It propagates HPO terms from neighbors in PPI network.
"""
import json
from collections import defaultdict
from functools import reduce
from src.utils.ontology import HumanPhenotypeOntology
from src.utils.file_reader import load_protein, load_annotation


def neighbor_scoring(network, test_proteins, train_annotation):
    """Scoring function of Neighbor method.
    :param network: protein-protein interaction network, like
        { protein1: { protein_a: score1a, ... }, ... }
    :param test_proteins: list of proteins in test set
        [ protein1, protein2, ... ]
    :param train_annotation: HPO annotations of training set
        { protein1: [ hpo_term1, hpo_term2, ... ], ... }
    :return: predictive score, like
        { protein1: { hpo_term1: score1, ... }, ... }
    """
    scores = defaultdict(dict)
    for protein in test_proteins:
        if protein in network:
            hpo_terms = reduce(lambda a, b: a | b,
                               [set(train_annotation.get(neighbour, set()))
                                for neighbour in network[protein]])
            normalizer = sum(network[protein].values())
            for hpo_term in hpo_terms:
                scores[protein][hpo_term] = sum(
                    [(hpo_term in train_annotation.get(neighbour, set())) *
                     network[protein][neighbour]
                     for neighbour in network[protein]]) / normalizer
    return scores


def add_weight(network):
    """Add weight to the edges in the network.
    :param network: PPI network (but scores are 0/1)
        { protein1: { protein1a: score1a, protein1b: score1b, ... },
          protein2: { protein2a: score2a, protein2b: score2b, ... },
          ... }
    :return: weighted PPI network, dict
        { protein1: { protein1a: score1a, protein1b: score1b, ... },
          protein2: { protein2a: score2a, protein2b: score2b, ... },
          ... }
    """
    ppi = defaultdict(dict)
    for protein_a in network:
        neighbour_a = set(network[protein_a].keys()) | {protein_a}
        for protein_b in network[protein_a]:
            neighbour_b = set(network[protein_b].keys()) | {protein_b}
            score = ((2 * len(neighbour_a & neighbour_b)) /
                     (len(neighbour_a - neighbour_b) +
                      2 * len(neighbour_a & neighbour_b) + 1)) * \
                    ((2 * len(neighbour_a & neighbour_b)) /
                     (len(neighbour_b - neighbour_a) +
                      2 * len(neighbour_a & neighbour_b) + 1))
            ppi[protein_a][protein_b] = score
            ppi[protein_b][protein_a] = score
    return ppi


if __name__ == "__main__":
    with open("../../../config/basic/neighbor/neighbor_COXPRESdb.json") as fp:
        config = json.load(fp)

    # load PPI network
    with open(config["network"]["path"]) as fp:
        network = json.load(fp)
    # actually customized for BioGRID
    if config["network"]["type"] == "unweighted":
        network = add_weight(network)

    # load proteins in training set and test set
    ltr_proteins = load_protein(config["protein_list"]["ltr"])
    test_proteins = load_protein(config["protein_list"]["test"])

    # load HPO
    ontology = HumanPhenotypeOntology(config["ontology"]["path"],
                                      version=config["ontology"]["version"])
    # load HPO annotations of training set
    train_annotation = load_annotation(config["annotation"], ontology, ns="all")

    # scoring by Neighbor method
    ltr_result = neighbor_scoring(network, ltr_proteins, train_annotation)
    test_result = neighbor_scoring(network, test_proteins, train_annotation)

    # write into file
    with open(config["result"]["ltr"], 'w') as fp:
        json.dump(ltr_result, fp, indent=2)
    with open(config["result"]["test"], 'w') as fp:
        json.dump(test_result, fp, indent=2)
