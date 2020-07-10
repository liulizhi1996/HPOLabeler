#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate the similarity of proteins by co-expression.

You can download co-expression data from
    https://coxpresdb.jp/download/
Choose "Species" to "Human" and "Method" to "R".
Besides, please download Entrez GeneID to UniProt Protein ID mapping file from
    https://www.uniprot.org/mapping/
"""
import json
import zipfile
from collections import defaultdict
import pandas as pd


def get_entrez_mapping(path_to_file):
    """Extract mapping from Entrez GeneID to UniProt Protein ID.
    :param path_to_file: path to the mapping file
        The Entrez GeneID to UniProt Protein ID mapping file is downloaded
        from https://www.uniprot.org/mapping/.
    :return: dict, key: Entrez GeneID, value: set of UniProt Protein IDs
        { gene1: { protein1a, protein1b, ... }, ... }
    """
    entrez2uniprot = defaultdict(set)
    with open(path_to_file) as fp:
        for line in fp:
            if line.startswith("Entry"):
                continue
            protein_id, gene_ids = line.split('\t')
            if len(gene_ids) == 0:
                continue
            else:
                id_list = gene_ids.strip().split(';')
                for gene_id in id_list:
                    if len(gene_id) > 0:
                        entrez2uniprot[gene_id].add(protein_id)
    return entrez2uniprot


def get_co_expression(path_to_file, entrez2uniprot):
    """Calculate proteins' similarity with respect to Co-expression.
    :param path_to_file: path to the COXPRESdb .zip file. You can download
        .zip file from https://coxpresdb.jp/download/ with "Method" being "R".
    :param entrez2uniprot: dict, map from Entrez GeneID to UniProt Protein ID,
        like { gene1: { protein1a, protein1b, ... }, ... }
    :return: dict, similarity of proteins, like
        { protein1: { protein1a: score1a, protein1b, score1b, ... }, ... }
    """
    ranks = defaultdict(dict)
    # open the .zip file
    with zipfile.ZipFile(path_to_file) as fp:
        # traverse each compressed file in the .zip file
        for filename in fp.namelist():
            # read the file
            text = fp.read(filename).decode('ascii').strip()
            # pass the directory (not file)
            if len(text) == 0:
                continue
            # extract Entrez GeneID from the filename
            gene_a = filename.split('/')[-1]
            # skip the gene without mapped protein
            if gene_a not in entrez2uniprot:
                continue
            # parse each line of the file
            lines = text.split('\n')
            for line in lines:
                gene_b, score = line.split('\t')
                for protein_b in entrez2uniprot[gene_b]:
                    for protein_a in entrez2uniprot[gene_a]:
                        ranks[protein_a][protein_b] = float(score)
    # transform the dict to DataFrame
    df_ranks = pd.DataFrame.from_dict(ranks).fillna(0)
    # normalize the score and convert it to the similarity
    # the closer the score is to 1, the more similar the two proteins are
    df_scores = pd.DataFrame(1 - df_ranks.values / df_ranks.values.max(),
                             columns=df_ranks.columns, index=df_ranks.index)
    # back to dict
    scores = df_scores.to_dict()
    return scores


if __name__ == "__main__":
    with open("../../../config/feature/COXPRESdb/coxpresdb.json") as fp:
        config = json.load(fp)

    # get Entrez GeneID to UniProt Protein ID mapping
    gene_mapping = get_entrez_mapping(config["entrez-uniprot"])
    # get Co-expression correlation scores
    coexp = get_co_expression(config["co-expression"], gene_mapping)
    # write into file
    with open(config["output"], 'w') as fp:
        json.dump(coexp, fp, indent=2)
