#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaluate performance of prediction methods, including F-max, AUROC and AUPR.
"""
import json
from collections import defaultdict
import pandas as pd
import numpy as np
from sklearn import metrics
from src.utils.ontology import HumanPhenotypeOntology, get_ns_id
from src.utils.file_reader import load_annotation, load_label_list, load_result

# HPO terms' group id according to frequency
frequency_group = ["very_rare",         # 1-3
                   "rare",              # 4-10
                   "uncommon",          # 11-30
                   "common",            # 31-100
                   "very_common",       # 101-300
                   "extremely_common"]  # >300


def f_max(result, annotation, n_threshold=101):
    """Calculate F-max and the corresponding threshold
    :param result: predictive scores, DataFrame like
        { protein1: { hpo_term1: score1, ... }, ... }
    :param annotation: true HPO annotations, DataFrame like
        { protein1: { hpo_term1: 0/1, ... }, ... } (0: no, 1: yes)
    :param n_threshold: number of thresholds, default: 101 (i.e. step=0.01)
    :return: F-max and the corresponding threshold
    """
    f_max_list = list()

    # vary the threshold over the range between 0 to 1 with n_threshold steps
    for threshold in np.linspace(0., 1., n_threshold):
        n_proteins, n_proteins_covered = 0, 0
        precision_ths, recall_ths = 0, 0

        # calculate precision and recall of each protein
        for protein in annotation.axes[0]:
            # HPO terms retrieved above the given threshold
            retrieved = set(result.loc[protein][result.loc[protein] >= threshold].index)
            # HPO terms associated with the given protein
            relevant = set(annotation.loc[protein][annotation.loc[protein] == 1].index)

            if len(relevant) > 0:
                # total of proteins having HPO annotations
                n_proteins += 1
                if len(retrieved) > 0:
                    # total of proteins annotated at least on predicted HPO term
                    n_proteins_covered += 1
                    precision_ths += len(retrieved & relevant) / len(retrieved)
                recall_ths += len(retrieved & relevant) / len(relevant)

        try:
            precision_ths = precision_ths / n_proteins_covered
        except ZeroDivisionError:
            precision_ths = 0
        try:
            recall_ths /= n_proteins
        except ZeroDivisionError:
            recall_ths = 0
        try:
            f_max_ths = (2*precision_ths*recall_ths)/(precision_ths+recall_ths)
        except ZeroDivisionError:
            f_max_ths = 0
        f_max_list.append(f_max_ths)

    # search the highest F-max value
    f_max_overall = max(f_max_list)
    # get the corresponding threshold
    threshold = np.argmax(f_max_list) / (n_threshold - 1)

    return f_max_overall, threshold


def auroc(result, annotation):
    """Calculate term-centric AUROC
    :param result: predictive scores, DataFrame like
        { protein1: { hpo_term1: score1, ... }, ... }
    :param annotation: true HPO annotations, DataFrame like
        { protein1: { hpo_term1: 0/1, ... }, ... } (0: no, 1: yes)
    :return: term-centric AUROC
    """
    auc = 0
    n_terms = 0

    # scan each HPO term
    for hpo_term in annotation.axes[1]:
        y_true = annotation[hpo_term].fillna(0).tolist()
        y_pred = result[hpo_term].fillna(0).tolist()

        # if only one class
        if len(np.unique(y_true)) != 2:
            continue

        n_terms += 1
        auc += metrics.roc_auc_score(y_true, y_pred)
    # average over terms
    auc = auc / n_terms if n_terms > 0 else 0

    return auc


def aupr(result, annotation):
    """Calculate pairwise AUPR
    :param result: predictive scores, DataFrame like
        { protein1: { hpo_term1: score1, ... }, ... }
    :param annotation: true HPO annotations, DataFrame like
        { protein1: { hpo_term1: 0/1, ... }, ... } (0: no, 1: yes)
    :return: pairwise AUPR
    """
    y_true = annotation.values.reshape(-1)
    y_pred = result.values.reshape(-1)
    aupr_value = metrics.average_precision_score(y_true, y_pred)
    return aupr_value


if __name__ == "__main__":
    with open("../../config/utils/evaluation/evaluation.json") as fp:
        config = json.load(fp)

    # load HPO
    ontology = HumanPhenotypeOntology(config["ontology"]["path"],
                                      version=config["ontology"]["version"])
    # get namespace list
    ns_id = get_ns_id(version=config["ontology"]["version"])

    # HPO annotations separated into sub-ontology
    # key: namespace of sub-ontology
    # value: DataFrame, row: protein, column: HPO term, value: 0/1
    df_annotation = dict()
    for ns in ns_id + ["all"]:
        # load propagated HPO annotations of specified sub-ontology
        annotation_ns = load_annotation(config["annotation"], ontology, ns=ns)
        # transform it to double-layer dict, i.e.
        # { protein1: { hpo_term1: score1, hpo_term2: score2, ... }, ... }
        dict_annotation_ns = defaultdict(dict)
        for protein in annotation_ns:
            for hpo_term in annotation_ns[protein]:
                dict_annotation_ns[protein][hpo_term] = 1
        # convert annotation to DataFrame
        df_annotation_ns = pd.DataFrame.from_dict(dict_annotation_ns,
                                                  orient="index")
        df_annotation[ns] = df_annotation_ns.fillna(0)

    # load HPO terms list according to frequency
    frequency = dict()
    for group_id in frequency_group:
        frequency[group_id] = load_label_list(config["frequency"][group_id])

    performance = dict()
    # evaluate each result in the list
    for res in config["result"]:
        performance[res] = dict()
        # load prediction result
        result = load_result(res)
        df_result = pd.DataFrame.from_dict(result, orient="index")
        # separate into sub-ontology
        for ns in ns_id:
            if ns == "freq":
                continue

            performance[res][ns] = dict()
            # get result only in one sub-ontology
            df_result_ns = df_result.reindex_like(df_annotation[ns])
            df_result_ns = df_result_ns.fillna(0)

            # calculate F-max
            f_max_ns, threshold_ns = f_max(df_result_ns, df_annotation[ns])
            performance[res][ns]['f_max'] = f_max_ns
            performance[res][ns]['threshold'] = threshold_ns

            # calculate term-centric AUROC
            term_auroc_ns = auroc(df_result_ns, df_annotation[ns])
            performance[res][ns]['auroc'] = term_auroc_ns

            # calculate pairwise AUPR
            pair_aupr_ns = aupr(df_result_ns, df_annotation[ns])
            performance[res][ns]['aupr'] = pair_aupr_ns

        # calculate AUROC by groups
        performance[res]["frequency"] = dict()
        for group_id in frequency_group:
            # select annotations of sub-group of HPO terms
            df_annotation_freq = df_annotation["all"].reindex(
                columns=frequency[group_id], fill_value=0)
            # select prediction result of sub-group of HPO terms
            df_result_freq = df_result.reindex_like(df_annotation_freq)
            df_result_freq = df_result_freq.fillna(0)

            # calculate term-centric AUROC according to frequency
            term_auroc_freq = auroc(df_result_freq, df_annotation_freq)
            performance[res]["frequency"][group_id] = term_auroc_freq

    # pretty output
    for res in performance:
        print(res)
        for ns in ns_id:
            if ns == "freq":
                continue
            print(ns, round(performance[res][ns]['f_max'], 4),
                  round(performance[res][ns]['threshold'], 4),
                  round(performance[res][ns]['auroc'], 4),
                  round(performance[res][ns]['aupr'], 4), sep='\t')
        for group_id in frequency_group:
            print(round(performance[res]['frequency'][group_id], 4), end='\t')
        print()
