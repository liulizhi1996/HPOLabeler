#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Flat classification trained and tested on each HPO term.
"""
import json
import time
from collections import defaultdict
import pandas as pd
from scipy import sparse
import numpy as np
from sklearn.linear_model import LogisticRegression
from src.utils.ontology import HumanPhenotypeOntology, get_ns_id
from src.utils.file_reader import load_protein, load_annotation, load_feature


def df_to_csr(df):
    """Convert Pandas DataFrame to SciPy sparse matrix.
    :param df: a Pandas DataFrame
    :return: the contents of the frame as a sparse SciPy CSR matrix
    """
    return sparse.csr_matrix(df.values)


class SameModel:
    """Works when only one class in ground truth.
    """
    def __init__(self):
        self.value = 0

    def fit(self, X, y):
        """Fill predictive score with the label (0/1) in label vector y.
        :param X: no use here
        :param y: label vector, can be list, numpy array or pd.DataFrame
        :return: None
        """
        if isinstance(y, pd.DataFrame):
            y = np.asarray(y)[:, 0]
        self.value = y[0]

    def predict(self, X):
        """Predict score using the label in label vector (0/1).
        :param X: feature matrix, its size of rows is useful here.
        :return: a numpy matrix of shape (n_samples, 2), the first column is
            false label's (i.e. 0), the second column is true label's (i.e. 1)
        """
        if isinstance(X, pd.DataFrame):
            X = np.asarray(X)
        return np.ones((X.shape[0], 2)) * [1 - self.value, self.value]

    def predict_proba(self, X):
        """Reture probability score of each protein on each HPO term.
        :param X: feature matrix (actually no use)
        :return: predictive scores, see predict()
        """
        return self.predict(X)


class FlatModel:
    def __init__(self, model):
        self._model = model
        self._classifiers = dict()

    def _get_model(self):
        """Return model prototype you need.
        :return: model prototype
        """
        if self._model == "lr":
            return LogisticRegression()
        else:
            raise ValueError("Can't recognize the model %s" % self._model)

    def fit(self, feature, annotation):
        """Fit the model according to the given feature and HPO annotations.

        N.B. The number of proteins in feature and annotation are MUST be the
            SAME!!!
        :param feature: features, DataFrame instance with rows being proteins
            and columns being HPO terms, the values are real number
        :param annotation: HPO annotations, DataFrame instance with rows being
            proteins and columns being HPO terms, the values are 0/1
        :return: None
        """
        assert isinstance(feature, pd.DataFrame), \
            "Argument feature must be Pandas DataFrame instance."
        assert isinstance(annotation, pd.DataFrame), \
            "Argument annotation must be Pandas DataFrame instance."
        assert feature.shape[0] == annotation.shape[0], \
            "The number of proteins in feature and annotation are must be " \
            "the same."

        X = df_to_csr(feature)
        for hpo_term in annotation.columns:
            y = np.asarray(annotation[[hpo_term]])[:, 0]
            if len(np.unique(y)) == 2:
                clf = self._get_model()
            else:
                clf = SameModel()
            clf.fit(X, y)
            self._classifiers[hpo_term] = clf
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                  "Fit", hpo_term)

    def predict(self, feature):
        """Predict scores on each HPO terms according to given features.
        :param feature: features, DataFrame instance with rows being proteins
            and columns being HPO terms, the values are real number
        :return: predictive score, dict like
        { protein1: { term1: score1, term2: score2
        """
        assert isinstance(feature, pd.DataFrame), \
            "Argument feature must be Pandas DataFrame instance."

        score = defaultdict(dict)
        protein_list = feature.axes[0].tolist()
        for hpo_term in self._classifiers:
            clf = self._classifiers[hpo_term]
            prediction = clf.predict_proba(feature)[:, 1]
            for idx, protein in enumerate(protein_list):
                score[protein][hpo_term] = prediction[idx]
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                  "Predict", hpo_term)

        return score


if __name__ == "__main__":
    with open("../../../config/basic/flat/flat_HIPPIE.json") as fp:
        config = json.load(fp)

    # load HPO
    ontology = HumanPhenotypeOntology(config["ontology"]["path"],
                                      version=config["ontology"]["version"])
    # get namespace id list
    ns_id = get_ns_id(version=config["ontology"]["version"])

    # load training set, ltr training set and test set
    train_protein_list = load_protein(config["protein_list"]["train"])
    ltr_protein_list = load_protein(config["protein_list"]["ltr"])
    test_protein_list = load_protein(config["protein_list"]["test"])

    # load features and convert them to DataFrame
    feature = load_feature(config["feature"])
    df_feature = pd.DataFrame.from_dict(feature, orient="index")
    df_feature = df_feature.fillna(0)

    combined_ltr_result = defaultdict(dict)
    combined_test_result = defaultdict(dict)
    for ns in ns_id:
        # load propagated HPO annotations of specified sub-ontology
        hpo_annotation = load_annotation(config["annotation"], ontology, ns)
        # transform it to double-layer dict, i.e.
        # { protein1: { hpo_term1: score1, hpo_term2: score2, ... }, ... }
        dict_annotation = defaultdict(dict)
        for protein in hpo_annotation:
            for hpo_term in hpo_annotation[protein]:
                dict_annotation[protein][hpo_term] = 1
        # convert annotation to DataFrame
        df_annotation = pd.DataFrame.from_dict(dict_annotation, orient="index")
        df_annotation = df_annotation.fillna(0)

        # extract training features and annotations
        train_protein_of_ns = list(set(train_protein_list) &
                                   set(df_annotation.axes[0]) &
                                   set(df_feature.axes[0]))
        train_feature = df_feature.loc[train_protein_of_ns]
        train_annotation = df_annotation.loc[train_protein_of_ns]
        # extract ltr training features and annotations
        ltr_protein_of_ns = list(set(ltr_protein_list) &
                                 set(df_feature.axes[0]))
        ltr_feature = df_feature.loc[ltr_protein_of_ns]
        # extract test features and annotations
        test_protein_of_ns = list(set(test_protein_list) &
                                  set(df_feature.axes[0]))
        test_feature = df_feature.loc[test_protein_of_ns]

        # train model and predict
        classifier = FlatModel(model=config["model"])
        classifier.fit(train_feature, train_annotation)
        ltr_result = classifier.predict(ltr_feature)
        test_result = classifier.predict(test_feature)

        # combine result into final result
        for protein in ltr_result:
            for term in ltr_result[protein]:
                combined_ltr_result[protein][term] = ltr_result[protein][term]
        for protein in test_result:
            for term in test_result[protein]:
                combined_test_result[protein][term] = test_result[protein][term]

    # write result
    with open(config["result"]["ltr"], "w") as fp:
        json.dump(combined_ltr_result, fp, indent=2)
    with open(config["result"]["test"], "w") as fp:
        json.dump(combined_test_result, fp, indent=2)
