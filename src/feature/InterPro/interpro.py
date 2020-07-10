#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parse the output xml InterPro matches file.

Please read https://github.com/ebi-pf-team/interproscan/wiki/HowToRun
to understand the way to run InterProScan:
    ./interproscan.sh -i /path/to/sequences.fasta -b /path/to/output_file -f XML
To figure out the format of output file, please refer to
https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats.
"""
import json
import xml.etree.ElementTree as ET


def get_ipr_annotation(path_to_ipr):
    """Parse the InterPro annotation file.
    :param path_to_ipr: path to the xml output file
    :return: protein2ipr, dict, containing the InterPro annotations
        of each protein, like
        {protein1: {ipr_1: 1, ipr_2: 1, ...}, ...}
    """
    protein2ipr = dict()
    # namespace of xml file
    namespace = {"xmlns": "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"}
    # read in xml file
    tree = ET.parse(path_to_ipr)
    root = tree.getroot()
    # parse each section of protein
    for protein_sec in root.findall("xmlns:protein", namespace):
        features = set()
        # matches section of this protein
        matches_sec = protein_sec.find("xmlns:matches", namespace)
        for match_sec in matches_sec:
            # each signature matched of this protein
            signature = match_sec.find("xmlns:signature", namespace)
            accession = signature.get("ac")
            features.add(accession)
        # get the ID of this protein
        for xref in protein_sec.findall("xmlns:xref", namespace):
            xref_id = xref.get("id").split('|')[1]
            protein2ipr[xref_id] = {ipr: 1 for ipr in features}
    return protein2ipr


if __name__ == "__main__":
    with open("../../../config/feature/InterPro/interpro.json") as fp:
        config = json.load(fp)

    # parse the xml file and obtain the InterPro annotation
    ipr_annotation = get_ipr_annotation(config["ipr_path"])

    # write into file
    with open(config["output"], 'w') as fp:
        json.dump(ipr_annotation, fp, indent=2)
