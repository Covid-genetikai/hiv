#!/usr/bin/env python
# coding: utf-8

# In[60]:


from Bio import GenBank
from pathlib import Path
import pandas as pd


# In[22]:


data_path = "/hiv/data"
data_path = "/home/mantydze/hiv/hiv"


# In[58]:


def parse_file(filename):
    data = []
    with open(filename, "r") as handle:
        for record in GenBank.parse(handle):

            # Mums idomus tik tiek kur yra 'HIV' SOURCE eilutej
            ## SOURCE      Human immunodeficiency virus 1 (HIV-1)
            if 'HIV' in record.source:

                obj = {
                    "accession": record.accession[0], 
                    'length': len(record.sequence),
                    "sequence": record.sequence       
                    
                }
                
                #[Feature(key='source', location='1..1686'), Feature(key='gene', location='<1..>1686'), Feature(key='CDS', location='<1..>1686')]
                for feature in record.features:
#                     print(feature.key)
                    if feature.key == "source":
                        # [Qualifier(key='/organism=', value='"Human immunodeficiency virus 1"'), Qualifier(key='/proviral', value=''), Qualifier(key='/mol_type=', value='"genomic DNA"'), Qualifier(key='/db_xref=', value='"taxon:11676"'), Qualifier(key='/country=', value='"Spain"')]
#                         for qualifier in qualifiers:
                                
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
                        
                       # genu pavadinimuose yra nemazai siuksliu, camel case ir t.t.
                    
                       # ['pol', 'env', 'gag',
                       #'vpr', 'vif', 'tat', 'rev', 'vpu', 'nef', 'gag-pol', 'rev1', 'gp120',
                       #'reverse transcriptase', 'HIV-1 protease', 'RT', 'protease', 'PR',
                       #'p24', 'envelope', 'gp160', 'tat/rev', 'vpu*', 'rt', 'ORF',
                       #'v-1 PROTEASE', 'Gag', 'Pol', 'Vif', 'Vpr', 'Tat', 'Rev', 'Vpu', 'Env',
                       #'Nef', 'pro', 'as', 'v-1 reverse transcriptase', 'v-1 protease', 'Pro',
                       #'V-1 protease', 'V-1 reverse transcriptase', 'polyprotein', 'rnv', 'gg',
                       #'env V3', 'gp41', 'GP160', 'GP120', 'GP41', 'env gene', 'RAK alpha',
                       #'V3', 'c2v3', 'vpx']
                        
                        if protein and gene and gene in ['pol', 'env', 'gag', 'vpr', 'vif', 'tat', 'rev', 'vpu', 'nef']:
                            obj[gene] = protein                           
                                
                    
                        
#                     print(feature)
#                     print(dir(feature))
                    
#                     print("---")
                    
                
                data.append(obj)
                
    return data
                


# In[61]:


data = []
for filename in Path(data_path).glob("*.seq"):
    print(f"Parsing {filename}")
    data.extend(parse_file(filename))


# In[62]:


df = pd.DataFrame(data)
df


# In[69]:


df.notnull().sum()

