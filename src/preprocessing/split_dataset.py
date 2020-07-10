#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Split protein list into three parts: train set, ltr set, and test set.
We will get three annotation datasets, three protein lists, and term list.
Besides, we will split HPO terms into several groups according to frequency.
"""
import json
from collections import defaultdict
from functools import reduce
from src.utils.ontology import HumanPhenotypeOntology
from src.utils.ontology import get_root, get_subontology


with open("../../config/preprocessing/split_dataset.json") as fp:
    config = json.load(fp)

# load various versions of HPO annotations
with open(config["raw_annotation"]["time0"]) as fp:
    annotation_t0 = json.load(fp)
with open(config["raw_annotation"]["time1"]) as fp:
    annotation_t1 = json.load(fp)
with open(config["raw_annotation"]["time2"]) as fp:
    annotation_t2 = json.load(fp)

# load various versions of HPO
ontology_t0 = HumanPhenotypeOntology(config["ontology"]["time0"]["path"],
                                     version=config["ontology"]["time0"]["version"])
ontology_t1 = HumanPhenotypeOntology(config["ontology"]["time1"]["path"],
                                     version=config["ontology"]["time1"]["version"])
ontology_t2 = HumanPhenotypeOntology(config["ontology"]["time2"]["path"],
                                     version=config["ontology"]["time2"]["version"])

# split proteins in train set, ltr set and test set
train_protein = list(annotation_t0.keys())
ltr_protein = list(set(annotation_t1.keys()) - set(annotation_t0.keys()))
test_protein = list(set(annotation_t2.keys()) - set(annotation_t1.keys()))

# HPO annotations for training
train_annotation = annotation_t0

# HPO annotations for learning to rank
ltr_annotation = defaultdict(list)
for protein in ltr_protein:
    for hpo_term in annotation_t1[protein]:
        # if "veteran" HPO term, just copy down
        if hpo_term in ontology_t0:
            ltr_annotation[protein].append(hpo_term)
        # else if has old, alternative HPO terms, replace it
        elif hpo_term in ontology_t1.alt_ids:
            for alternative in ontology_t1.alt_ids[hpo_term]:
                print("Replace (%s, %s) to (%s, %s)" % (protein, hpo_term, protein, alternative))
                ltr_annotation[protein].append(alternative)
        # if not found, then discard
        else:
            print("Discard (%s, %s)" % (protein, hpo_term))

# HPO annotations for evaluation
test_annotation = defaultdict(list)
for protein in test_protein:
    for hpo_term in annotation_t2[protein]:
        # if "veteran" HPO term, just copy down
        if hpo_term in ontology_t0:
            test_annotation[protein].append(hpo_term)
        # else if has old, alternative HPO terms, replace it
        elif hpo_term in ontology_t2.alt_ids:
            for alternative in ontology_t2.alt_ids[hpo_term]:
                print("Replace (%s, %s) to (%s, %s)" % (protein, hpo_term, protein, alternative))
                test_annotation[protein].append(alternative)
        # if not found, then discard
        else:
            print("Discard (%s, %s)" % (protein, hpo_term))

# output annotation
with open(config["processed_annotation"]["train"], 'w') as fp:
    json.dump(train_annotation, fp, indent=2)
with open(config["processed_annotation"]["ltr"], 'w') as fp:
    json.dump(ltr_annotation, fp, indent=2)
with open(config["processed_annotation"]["test"], 'w') as fp:
    json.dump(test_annotation, fp, indent=2)

# output protein lists
with open(config["protein_list"]["train"], 'w') as fp:
    json.dump(list(train_annotation.keys()), fp, indent=2)
with open(config["protein_list"]["ltr"], 'w') as fp:
    json.dump(list(ltr_annotation.keys()), fp, indent=2)
with open(config["protein_list"]["test"], 'w') as fp:
    json.dump(list(test_annotation.keys()), fp, indent=2)

# merge three annotations into one dict
combined_annotation = {**train_annotation,
                       **ltr_annotation,
                       **test_annotation}

# propagate annotations and discard roots (All & sub-ontology)
propagated_combined_annotation = dict()
for protein in combined_annotation:
    propagated_combined_annotation[protein] = list(
        ontology_t0.transfer(combined_annotation[protein]) - {get_root()} -
        set(get_subontology(ontology_t0.version)))

# write down used HPO terms into file
with open(config["term_list"], 'w') as fp:
    json.dump(list(reduce(lambda a, b: set(a) | set(b),
                          propagated_combined_annotation.values())),
              fp, indent=2)

# split HPO terms according to its frequency
# 1. count the frequency of HPO terms
#    term_counts = { hpo_term1: cnt1, hpo_term2: cnt2, ... }
term_counts = defaultdict(int)
for protein in propagated_combined_annotation:
    for term in propagated_combined_annotation[protein]:
        term_counts[term] += 1
# 2. write term list each interval
for interval in config["frequency"]:
    with open(interval["path"], 'w') as fp:
        json.dump([term for term in term_counts
                   if interval["low"] <= term_counts[term] <= interval["high"]],
                  fp, indent=2)
