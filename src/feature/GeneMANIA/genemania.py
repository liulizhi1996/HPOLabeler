#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct GeneMANIA PPI network.

You can download identifier mapping file in
    http://genemania.org/data/current/Homo_sapiens/identifier_mappings.txt
and combined network data in
    http://genemania.org/data/current/Homo_sapiens.COMBINED/
    COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt
"""
import json
from collections import defaultdict


def ensembl2uniprot(file_path):
    """Map Ensembl gene id to UniProt accession.
    :param file_path: path to mapping file provided by GeneMANIA website:
        http://genemania.org/data/current/Homo_sapiens/identifier_mappings.txt
    :return: dict, key: Ensembl gene id, value: UniProt accession
        { ENSG1: uniprot_ac1, ENSG2: uniprot_ac2, ... }
    """
    mapping = dict()
    with open(file_path) as fp:
        for line in fp:
            if line.startswith("Preferred_Name"):
                continue
            ensembl_id, name, source = line.strip().split('\t')
            if "Uniprot ID" in source:
                mapping[ensembl_id] = name
    return mapping


def get_genemania_network(path_to_network, path_to_mapping):
    """Construct GeneMANIA PPI network.
    :param path_to_network: path to combined GeneMANIA network data
    :param path_to_mapping: path to identifier mapping file
    :return: dict, PPI network
        { protein1: { protein1a: score1a, protein1b: score1b, ... },
          protein2: { protein2a: score2a, protein2b: score2b, ... },
          ... }
    """
    network = defaultdict(dict)
    mapping = ensembl2uniprot(path_to_mapping)
    with open(path_to_network) as fp:
        for line in fp:
            if line.startswith("Gene_A"):
                continue
            gene1, gene2, score = line.strip().split()
            try:    # if no matched accession found, pass it
                protein1 = mapping[gene1]
                protein2 = mapping[gene2]
            except KeyError:
                continue
            score = float(score)
            network[protein1][protein2] = network[protein2][protein1] = score
    return network


if __name__ == "__main__":
    with open("../../../config/feature/GeneMANIA/genemania.json") as fp:
        config = json.load(fp)

    # get PPI network
    network = get_genemania_network(config["network"], config["mapping"])
    # write into file
    with open(config["feature"], 'w') as fp:
        json.dump(network, fp, indent=2)
