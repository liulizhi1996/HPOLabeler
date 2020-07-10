#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Definition of Gene Ontology.
"""
from collections import defaultdict
from src.utils.obo_parser import GODag


def get_short_ns(long_ns):
    """Return short namespace id of the given long namespace.
    :param long_ns: long name of namespace
    :return: short namespace id
    """
    if long_ns == "biological_process":
        return "bp"
    elif long_ns == "cellular_component":
        return "cc"
    elif long_ns == "molecular_function":
        return "mf"
    else:
        raise ValueError("Can't recognize %s" % long_ns)


def get_ns_id():
    """Return list of namespace id of each sub-ontology of GO.
    :return: namespace id list
    """
    return ["bp", "cc", "mf"]


def get_subontology():
    """Return term list of roots of each sub-ontology.
    :return: root terms list
    """
    return ['GO:0008150',   # biological process
            'GO:0005575',   # cellular component
            'GO:0003674']   # molecular function


def get_ns_id2hpo():
    """Return mapping of namespace to the root of sub-ontology.
    :return: mapping of namespace to GO term
    """
    return {"bp": 'GO:0008150',
            "cc": 'GO:0005575',
            "mf": 'GO:0003674'}


def get_hpo2ns_id():
    """Return mapping of namespaces of roots of sub-ontology.
    :return: mapping, key: root term of sub-ontology, value: namespace
    """
    return {'GO:0008150': "bp",
            'GO:0005575': "cc",
            'GO:0003674': "mf"}


class GOTerm(object):
    def __init__(self, go_term):
        """Definition of an GO term.
        Attributes:
            - id: accession of GO term, e.g. GO:0000005
            - parents: set, parents of this term
            - name: name of GO term
            - ns: namespace the term belongs to
            - children: set, children of this term
            - depth: the depth of term in the whole GO
        :param go_term: instance of the HPO term
        :return: None
        """
        self.id = go_term.id
        self.parents = set([p.id for p in go_term.parents])
        if hasattr(go_term, 'relationship'):
            for parent in go_term.relationship.get('part_of', set()):
                if parent.namespace == go_term.namespace:
                    self.parents.add(parent.id)
        self.name = go_term.name
        self.ns = get_short_ns(go_term.namespace)
        self.children = set()
        self.depth = 0


class GeneOntology(dict):
    """Definition of GO, it is inheritance of dict.
    Attributes:
        - alt_ids: dict, alternative ids of GO terms, like
            { go_term1: [ alt_id1, alt_id2, ... ], ... }
        - dict of GO terms, you can visit it like dict,
            e.g. ontology['GO:0000005']
    """
    def __init__(self, obo_file_path):
        """
        :param obo_file_path: path to obo file
        :return: None
        """
        super(GeneOntology, self).__init__()
        go_dag = GODag(obo_file_path, 'relationship')
        self.alt_ids = go_dag.alt_ids
        for go_id, go_term in go_dag.items():
            self[go_id] = GOTerm(go_term)
        self._get_children()
        self._get_depth()

    def _get_children(self):
        """Fill in children of each term.
        :return: None
        """
        for go_id in self:
            for parent in self[go_id].parents:
                self[parent].children.add(go_id)

    def _get_depth(self):
        """Fill in depth of each term.
        :return: None
        """
        for root in get_subontology():
            self[root].depth = 1
        now = set(get_subontology())
        while len(now) > 0:
            next = set()
            for go_term in now:
                for child in self[go_term].children:
                    if self[child].depth == 0:
                        next.add(child)
                        self[child].depth = self[go_term].depth + 1
            now = next

    def transfer(self, go_list):
        """Propagate GO terms by true-path-rule.
        :param go_list: the GO terms should be transferred
        :return: Propagated GO terms of go_list
        """
        go_list = list(filter(lambda x: x in self, go_list))
        ancestors, now = set(go_list), set(go_list)
        while len(now) > 0:
            next = set()
            for go_term in now:
                if go_term in self:
                    next |= self[go_term].parents - ancestors
            now = next
            ancestors |= now
        return ancestors

    def transfer_scores(self, term_scores):
        """Keep the consistency of predictive scores - the score of the term
        must not be greater than the score of all its child nodes.
        :param term_scores: dict, predictive scores of GO terms.
            { go_term1: score1, go_term2: score2, ... }
        :return: Scores that adhere to consistency constraints.
        """
        scores = defaultdict(float)
        for go_term in sorted(self.transfer(term_scores.keys()),
                              key=lambda x: self[x].depth, reverse=True):
            scores[go_term] = max(scores[go_term], term_scores.get(go_term, 0))
            for parent_id in self[go_term].parents:
                scores[parent_id] = max(scores[parent_id], scores[go_term])
        return scores
