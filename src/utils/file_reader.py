#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Readers of files with different formats.
"""
import json
from collections import defaultdict
from src.utils.ontology import get_root, get_subontology


def gene2uniprot(file_path, gene_column, uniprot_column):
    """Mapping entrez gene id to uniprot protein id.

    :param file_path: path to mapping file
    :param gene_column: the column index of gene id
    :param uniprot_column: the column index of uniprot id
    :return: a dict with key being gene id and value being list of uniprot ids
    { gene_id: [uniprot_id1, uniprot_id2, ...] }
    """
    gene_to_protein = defaultdict(list)
    with open(file_path) as fp:
        for line in fp:
            if line.startswith("y"):    # omit the header line
                continue
            entries = line.strip().split('\t')
            # multi-genes mapped to the same protein
            if ',' in entries[gene_column]:
                genes = entries[gene_column].split(',')
                protein = entries[uniprot_column]
                for gene in genes:
                    gene_to_protein[gene].append(protein)
            # one gene mapped to one protein
            else:
                gene, protein = entries[gene_column], entries[uniprot_column]
                gene_to_protein[gene].append(protein)
    return gene_to_protein


def load_annotation(file_path, ontology, ns="all"):
    """Get propagated HPO annotations.
    :param file_path: path to raw annotation
    :param ontology: instance of HumanPhenotypeOntology
    :param ns: namespace, valid values includes:
        - "all": return all annotations
        - specific namespace (e.g. "pa"): only return annotations of
            specified sub-ontology
    :return: annotations after propagation
    { protein1: [ hpo_term1, hpo_term2, ... ] ... }
    """
    # load raw annotations without propagation
    with open(file_path) as fp:
        leaf_annotation = json.load(fp)
    # filter annotations of specified namespace
    if ns != "all":
        filtered_leaf_annotation = defaultdict(list)
        for protein in leaf_annotation:
            filtered = list(filter(lambda t: ontology[t].ns == ns,
                                   leaf_annotation[protein]))
            if len(filtered) > 0:
                filtered_leaf_annotation[protein] = filtered
        leaf_annotation = filtered_leaf_annotation
    # propagate annotations and discard roots (All & sub-ontology)
    propagated_annotation = dict()
    for protein in leaf_annotation:
        propagated_annotation[protein] = list(
            ontology.transfer(leaf_annotation[protein]) - {get_root()} -
            set(get_subontology(ontology.version)))
    return propagated_annotation


def load_protein(file_path):
    """Return protein list.
    :param file_path: path to protein list file
    :return: list of proteins
    [ protein1, protein2, ... ]
    """
    with open(file_path) as fp:
        protein_list = json.load(fp)
    return protein_list


def load_feature(file_path):
    """Load features into a dict.
    :param file_path: path to feature file
    :return: dict,
    { protein1: { feature1: score1, feature2: score2, ... } ... }
    """
    with open(file_path) as fp:
        feature = json.load(fp)
    return feature


def load_result(file_path):
    """Load prediction results.
    :param file_path: path to prediction result file
    :return: dict, like
    { protein1: { hpo_term1: score1, hpo_term2: score2, ... } ... }
    """
    with open(file_path) as fp:
        result = json.load(fp)
    return result


def load_label_list(file_path):
    """Load list of HPO terms.
    :param file_path: path to label list file, one HPO term a line
    :return: list, like [ hpo_term1, hpo_term2, ...]
    """
    with open(file_path) as fp:
        label_list = json.load(fp)
    return label_list
