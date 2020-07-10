#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct Protein-Protein Interaction network from HIPPIE database.

Download data set from
    http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
with "HIPPIE tab format".
Then download UniProt Entry name to ID mapping file from
    https://www.uniprot.org/mapping/
"""
import json
from collections import defaultdict


def get_mapping(path_to_file):
    """Generate mapping from UniProt Entry name to ID.
    :param path_to_file: path to the mapping file
    :return: dict, like
        { name1: ac1, name2: ac2, ... }
    """
    name2id = dict()
    with open(path_to_file) as fp:
        for line in fp:
            if line.startswith("Entry"):
                continue
            accession, name = line.strip().split('\t')
            name2id[name] = accession
    return name2id


def get_network(path_to_file, name2id):
    """Construct PPI network.
    :param path_to_file: path to the PPI data set
    :param name2id: dict, mapping from UniProt Entry name to ID
    :return: nested dict containing PPI scores, like
        { protein1: { protein1a: score1a, protein1b: score1b, ... }, ... }
    """
    interaction = defaultdict(dict)
    with open(path_to_file) as fp:
        for line in fp:
            entries = line.strip().split('\t')
            name1, name2, score = entries[0], entries[2], entries[4]
            if name1 in name2id and name2 in name2id:
                protein1, protein2 = name2id[name1], name2id[name2]
                interaction[protein1][protein2] = float(score)
                interaction[protein2][protein1] = float(score)
    return interaction


if __name__ == "__main__":
    with open("../../../config/feature/HIPPIE/hippie.json") as fp:
        config = json.load(fp)

    # get UniProt name to ID mapping
    mapping = get_mapping(config["mapping"])
    # get PPI network
    network = get_network(config["network"], mapping)
    # write into file
    with open(config["output"], 'w') as fp:
        json.dump(network, fp, indent=2)
