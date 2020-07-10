#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct BioGRID PPI network.

You can download PPI network data from one directory you're interested in
    https://downloads.thebiogrid.org/BioGRID/Release-Archive
and download file with extension ".tab2" with "ALL" sources.
Meanwhile, the identifier mapping file can be downloaded from
    https://downloads.thebiogrid.org/Download/BioGRID/
    External-Database-Builds/UNIPROT.tab.txt
But there are no guarantees of version here.
"""
import json
from collections import defaultdict


def biogrid2uniprot(file_path):
    """Map BioGRID id to UniProt id.
    :param file_path: path to mapping file provided by BioGRID website
        https://downloads.thebiogrid.org/Download/BioGRID/
        External-Database-Builds/UNIPROT.tab.txt
        NOTE: no guarantees of version
    :return: dict, key: BioGRID id, value: UniProt accession
        { biogrid_ac1: uniprot_ac1, biogrid_ac2: uniprot_ac2, ... }
    """
    mapping = dict()
    with open(file_path) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            uniprot_id, biogrid_id, *_ = line.split()
            mapping[biogrid_id] = uniprot_id
    return mapping


def get_biogrid_network(path_to_network, path_to_mapping):
    """Construct BioGRID PPI network.
    :param path_to_network: path to BioGRID network data
        You can open https://downloads.thebiogrid.org/BioGRID/Release-Archive
        and choose a directory you want and download file with the extension
        ".tab2" with "ALL" sources
    :param path_to_mapping: path to mapping file provided by BioGRID website
    :return: dict, PPI network
        { protein1: { protein1a: score1a, protein1b: score1b, ... },
          protein2: { protein2a: score2a, protein2b: score2b, ... },
          ... }
    """
    network = defaultdict(dict)
    mapping = biogrid2uniprot(path_to_mapping)
    with open(path_to_network) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            entries = line.strip().split('\t')
            biogrid_a, biogrid_b = entries[3], entries[4]
            organism_a, organism_b = entries[15], entries[16]
            if organism_a != '9606' or organism_b != '9606':    # not human's
                continue
            try:  # if no matched accession found, pass it
                protein_a = mapping[biogrid_a]
                protein_b = mapping[biogrid_b]
            except KeyError:
                continue
            # BioGRID interaction doesn't provide confidence score (mostly),
            # so we construct unweighted graph here
            network[protein_a][protein_b] = network[protein_b][protein_a] = 1
    return network


if __name__ == "__main__":
    with open("../../../config/feature/BioGRID/biogrid.json") as fp:
        config = json.load(fp)

    # get PPI network
    network = get_biogrid_network(config["network"], config["mapping"])
    # write into file
    with open(config["feature"], 'w') as fp:
        json.dump(network, fp, indent=2)
