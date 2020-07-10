#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract trigram frequency of amino acid sequences.

You can get sequences from https://www.uniprot.org/mapping/.
"""
import json
from collections import defaultdict
from Bio import SeqIO

# 20 amino acids (some unknown amino acids are not shown in the following)
acids = "ACDEFGHIKLMNPQRSTVWY"
# map amino acids to integers
acids_idx = {acid: i for i, acid in enumerate(acids)}


def get_id(mers):
    """Return id of the mers.
    :param mers: n-gram of amino acid sequence
    :return: id of the given mers
    """
    index = 0
    for mer in mers:
        if mer not in acids_idx:
            return None
        index = index * len(acids) + acids_idx[mer]
    return index


def get_mer_features(seqs, k):
    """Return k-gram frequency of amino acid sequence in seqs.
    :param seqs: list of amino acid sequence
    :param k: size of gram
    :return: a dict like
    { protein1: { k-gram1: frequency1, k-gram2: frequency2, ... }, ... }
    """
    features = defaultdict(lambda: defaultdict(int))
    for seq in seqs:
        for i in range(len(seq.seq) - k):
            # slice k-gram
            mers = seq.seq[i:i+k]
            # map k-gram to id
            col_id = get_id(mers)
            if col_id is not None:
                features[seq.id][col_id] += 1
            else:       # unknown amino acid in sequence
                print(seq.id, mers)
    return features


if __name__ == "__main__":
    with open("../../../config/feature/Trigram/trigram.json") as fp:
        config = json.load(fp)

    # load amino acid sequence
    seqs = list(SeqIO.parse(config["fasta"], "fasta"))
    # extract UniPort id
    for seq in seqs:
        seq.id = seq.id.split('|')[1]
    # get features
    trigrams = get_mer_features(seqs, 3)
    # write into file
    with open(config["feature"], 'w') as fp:
        json.dump(trigrams, fp, indent=2)
