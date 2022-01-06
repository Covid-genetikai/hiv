#!/usr/bin/env python
# coding: utf-8

from Bio import GenBank
from pathlib import Path
import pandas as pd
import re, csv

data_path = "/hiv/data"
data_path = "/home/mantydze/hiv/hiv"
data_path = "data"

GENES = ['pol', 'env', 'gag', 'vpr', 'vif', 'tat', 'rev', 'vpu', 'nef']


def read_qualifiers(feature, *args):
    resp = [ None for n in args ]
    shit_names = [ f"/{n}=" for n in args ]
    for qualifier in feature.qualifiers:
        if qualifier.key in shit_names:
            i = shit_names.index(qualifier.key)
            resp[i] = qualifier.value.replace('"', '').replace("'", "")
    return tuple(resp)


def parse_record(record):

    obj = None

    if 'HIV' in record.source:

        obj = {
            "accession": record.accession[0],
            'length': len(record.sequence),
            "sequence": record.sequence
        }

        for gene in GENES:
            obj[gene + "_pro"] = None
            obj[gene + "_loc"] = None

        for feature in record.features:

            if feature.key == "source":
                obj["country"] = read_qualifiers(feature, 'country')[0]
                continue

            if feature.key == "gene":
                location = feature.location
                gene = read_qualifiers(feature, 'gene')[0]
                if gene is not None and location is not None:
                    gene = gene.lower()
                    if gene in GENES:
                        obj[gene + "_loc"] = location
                continue

            if feature.key == "CDS":

                resp = read_qualifiers(feature, 'gene', 'translation')
                gene, protein = resp[0], resp[1]

                if protein is not None and gene is not None and gene in GENES:
                    obj[gene + "_pro"] = protein
                continue

    return obj


def parse_file(filename):
    data = []
    with open(filename, "r") as handle:
        records = GenBank.parse(handle)
        while True:
            try:
                obj = parse_record(record=records.__next__())
                if obj is not None:
                    data.append(obj)

            except StopIteration as e:
                break
            except Exception as e:
                print(e, "Ignoring..")
                pass

    return data


for filename in Path(data_path).glob("*.seq"):

    print(f"Parsing {filename}")
    data = parse_file(filename)
    df = pd.DataFrame(data)

    filename = re.sub("\\.seq$",".csv", str(filename))
    df.to_csv(filename, quoting=csv.QUOTE_NONNUMERIC)
    print(f"Wrote {filename} len {len(df)}")

