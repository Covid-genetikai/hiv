#!/usr/bin/env python
# coding: utf-8

from Bio import GenBank
from pathlib import Path
import pandas as pd


data_path = "/data/hiv/data"
GENES = ['pol']#, 'env', 'gag', 'vpr', 'vif', 'tat', 'rev', 'vpu', 'nef']


def parse_record(record):
    obj = {}
    
    if 'HIV' in record.source:
        
#         if len(record.sequence) <= 200:
#             return None

        obj = {
            "accession": record.accession[0], 
            'length': len(record.sequence),
            'sequence': record.sequence       
        }

        #[Feature(key='source', location='1..1686'), Feature(key='gene', location='<1..>1686'), Feature(key='CDS', location='<1..>1686')]
        for feature in record.features:
            if feature.key == "source":
                # [Qualifier(key='/organism=', value='"Human immunodeficiency virus 1"'), Qualifier(key='/proviral', value=''), Qualifier(key='/mol_type=', value='"genomic DNA"'), Qualifier(key='/db_xref=', value='"taxon:11676"'), Qualifier(key='/country=', value='"Spain"')]
                for qualifier in feature.qualifiers:
                    if qualifier.key == "/country=":
                        obj["country"] = qualifier.value.replace('"','').replace("'", "")

            if feature.key == "gene":
                pass

            if feature.key == "CDS":
                gene = None
                protein = None
                
                for qualifier in feature.qualifiers:
                    if qualifier.key == "/gene=":
                        gene = qualifier.value.replace('"','').replace("'", "").lower()

                    if qualifier.key == "/translation=":
                        protein = qualifier.value.replace('"','').replace("'", "")
                    
                if protein and gene and gene in GENES:
                    obj[f"{gene}"] = protein    
#                     obj[f"{gene}_loc"] = feature.location.replace("<", "").replace(">", "")
    
    if len(obj.keys()) <= 4:
        # accession, length, sequence, country
        return None
    
    return obj 
    
    
def parse_file(filename):
    data = []
    with open(filename, "r") as handle:
        records = GenBank.parse(handle)
        while True:
            try:
                obj = parse_record(records.__next__())
                if obj:
                    data.append(obj)
            
            except StopIteration as e:
                break
                    
            except Exception as exc:
                print(exc)
                
    return data                


data = []
for filename in sorted(list(Path(data_path).glob("*.seq"))):
    print(f"Parsing {filename}")
    sub_data = parse_file(filename)
    if len(sub_data) > 0:
        data.extend(sub_data)
        sub_df = pd.DataFrame(sub_data)
        sub_df.to_csv(filename.with_suffix(".csv"))
        print("Done")
    else:
        print("No useful data")
        

df = pd.DataFrame(data)
df.to_csv(Path(data_path).joinpath("1-full.csv"), index = False)    


df = df[~df.duplicated(['sequence', 'pol'], keep='first')]


for gene in GENES:
    sub_df = df[["accession", gene]]
    sub_df.to_csv(Path(data_path).joinpath(f"1-{gene}.csv"), index=False)   
    
    with open(Path(data_path).joinpath(f"1-{gene}.fasta"), "w") as fasta_file:
        for _, row in df.iterrows():
            fasta_file.write(f">{row['accession']}\n")
            fasta_file.write(f"{row[gene]}\n")
    
