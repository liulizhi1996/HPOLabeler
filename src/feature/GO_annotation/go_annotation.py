#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construct GO annotation features.

Here, we first extract raw annotation from GOA
(see https://www.ebi.ac.uk/GOA/downloads) and then propagate via Gene Ontology
downloaded from http://geneontology.org/page/download-ontology.
"""
import json
from collections import defaultdict
from src.utils.gene_ontology import GeneOntology, get_subontology, get_ns_id


# Reliable annotation evidence codes
EvidenceCode = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC']


def goa(file_path):
    """Extract raw GO annotations from GOA.
    :param file_path: path to annotation (.gaf) file
        Please visit https://www.ebi.ac.uk/GOA/downloads and open GOA ftp site,
        enter "HUMAN" directory, download goa_human.gaf.gz.
    :return: dict, key: protein, value: list of GO terms
    { protein1: [ go_term1, go_term2, ... ], ... }
    """
    annotation = defaultdict(list)
    with open(file_path) as fp:
        for line in fp:
            if line.startswith('!'):
                continue
            entries = line.strip().split('\t')
            qualifier = entries[3]
            evidence_code = entries[6]
            db_object_type = entries[11]
            taxon = entries[12]
            if "NOT" not in qualifier and evidence_code in EvidenceCode and \
               db_object_type == "protein" and "taxon:9606" in taxon:
                protein = entries[1]
                go_term = entries[4]
                annotation[protein].append(go_term)
    return annotation


def propagate_go_annotation(leaf_annotation, ontology, ns="all"):
    """Return propagated GO annotations (or in specified sub-ontology).
    :param leaf_annotation: raw annotations, it's a dict:
        { protein1: [ go_term1, go_term2, ... ], ... }
    :param ontology: instance of GeneOntology
    :param ns: valid values:
        - "all": return all GO terms in Gene Ontology
        - specified namespace: only return terms in this sub-ontology
    :return: propagated GO annotations, dict:
    { protein1: [ go_term1, go_term2, ... ], ... }
    """
    # discard obsolete GO terms
    existed_leaf_annotation = defaultdict(list)
    for protein in leaf_annotation:
        existed_leaf_annotation[protein] = list(
            filter(lambda t: t in ontology, leaf_annotation[protein]))
    leaf_annotation = existed_leaf_annotation
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
            ontology.transfer(leaf_annotation[protein]) -
            set(get_subontology()))
    return propagated_annotation


def get_feature(annotation):
    """Transform annotation to feature structure like
    { protein1: { go_term1: 1, go_term2: 1, ... }, ... }
    :param annotation: dict, like
        { protein1: [ go_term1, go_term2, ... ], ... }
    :return: dict, but in different structure like
    { protein1: { go_term1: 1, go_term2: 1, ... }, ... }
    """
    feature = defaultdict(dict)
    for protein in annotation:
        for term in annotation[protein]:
            feature[protein][term] = 1
    return feature


if __name__ == "__main__":
    with open("../../../config/feature/GO_annotation/go_annotation.json") as fp:
        config = json.load(fp)

    # load GOA annotation
    goa_annotation = goa(config["goa"])
    # load GO
    ontology = GeneOntology(config["ontology"])
    for ns in get_ns_id():
        # propagate annotations
        full_annotation = propagate_go_annotation(goa_annotation, ontology,
                                                  ns=ns)
        # transform to feature structure
        annotation_feature = get_feature(full_annotation)
        # write into file
        with open(config["feature"][ns], 'w') as fp:
            json.dump(annotation_feature, fp, indent=2)
