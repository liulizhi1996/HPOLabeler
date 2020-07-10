#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Naive model: count frequency of HPO terms in given HPO annotations.
                   | {pro in training proteins | annotation[pro][term] = 1} |
S(protein, term) = ----------------------------------------------------------
                                 size of training proteins
"""
import json
from collections import defaultdict
from src.utils.ontology import HumanPhenotypeOntology, get_ns_id
from src.utils.file_reader import load_annotation, load_protein


class Naive:
    """Definition of naive model.

    Attributes:
        - _frequency (private): the frequency of HPO terms in train annotations
    """
    def __init__(self):
        self._frequency = defaultdict(float)

    def fit(self, annotation):
        """Count the frequency of HPO terms according to the given annotations.
        :param annotation: selected HPO annotations with structure like
            { protein1: [ hpo_term1, hpo_term2, ... ] ... }
        :return: None
        """
        for protein in annotation:
            for hpo_term in annotation[protein]:
                self._frequency[hpo_term] += 1
        n_protein = len(annotation)
        for hpo_term in self._frequency:
            self._frequency[hpo_term] /= n_protein

    def predict(self, protein_list):
        """Assign frequency to each protein in protein_list.
        :param protein_list: list of test proteins
        :return: predicted frequency of HPO terms:
        { protein1: { hpo_term1: score1, hpo_term2: score2, ... }, ... }
        """
        prediction = dict()
        for protein in protein_list:
            prediction[protein] = self._frequency
        return prediction


if __name__ == "__main__":
    with open("../../../config/basic/naive/naive.json") as fp:
        config = json.load(fp)

    # load HPO
    ontology = HumanPhenotypeOntology(config["ontology"]["path"],
                                      version=config["ontology"]["version"])
    # get namespace id list
    ns_id = get_ns_id(version=config["ontology"]["version"])

    # load training set and test set
    ltr_protein_list = load_protein(config["protein_list"]["ltr"])
    test_protein_list = load_protein(config["protein_list"]["test"])

    ltr_result = defaultdict(dict)
    test_result = defaultdict(dict)
    for ns in ns_id:
        # load propagated HPO annotations of specified sub-ontology
        hpo_annotation = load_annotation(config["annotation"], ontology, ns)

        predictor = Naive()
        predictor.fit(hpo_annotation)
        ltr_predict_result = predictor.predict(ltr_protein_list)
        test_predict_result = predictor.predict(test_protein_list)

        # combine the prediction
        for protein in ltr_predict_result:
            for hpo_term in ltr_predict_result[protein]:
                ltr_result[protein][hpo_term] = ltr_predict_result[protein][hpo_term]
        for protein in test_predict_result:
            for hpo_term in test_predict_result[protein]:
                test_result[protein][hpo_term] = test_predict_result[protein][hpo_term]

    # write into file
    with open(config["result"]["ltr"], 'w') as fp:
        json.dump(ltr_result, fp, indent=2)
    with open(config["result"]["test"], 'w') as fp:
        json.dump(test_result, fp, indent=2)
