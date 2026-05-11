"""
1) Input
"""

# Load differential expression results

import pandas as pd
input = pd.read_table('entrada.txt')
input.head()

"""
2) Pre-processing
"""

# Extract GeneID (CICLE_v1000xxxxmg), which is the valid gene identifier for Citrus clementina.
# The previous EntrezGeneID column is no longer usable for UniProt mapping.
lista = input['GeneID'].tolist()
lista = [x for x in lista if str(x) != 'nan']

# Select a subset of genes to query (recommended <100 to avoid server overload)
prots = lista[:25]

"""
UniProt access using REST API
"""

import requests

# Query UniProtKB using the REST API.
# Web scraping and the old 'uploadlists' endpoint are no longer supported.
def get_uniprot(query=''):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{query})",   # Exact match for EnsemblGenome/Phytozome gene IDs
        "format": "txt"                     # Returns plain text compatible with the existing parser
    }

    r = requests.get(url, params=params)
    r.raise_for_status()

    return r.text.splitlines()


"""
3) Processing
"""

table=pd.DataFrame()
full_data = []

# For each gene, parse UniProtKB text output
for index,entry in enumerate(prots):
    # Containers for parsed annotations
    uni_id = []
    pdbs=[]
    functions=[]
    process=[]
    full_name=[]
    kegg = []
    EC = []
    pfam = []
    prosite = []

    # Query UniProtKB for the gene
    data=get_uniprot(query=entry)
    # Append the data to the full_data list
    full_data.append(data)
    # Add the gene ID to the table
    table.loc[index,'GeneID'] = entry
    
#For each line in the output: search for information from each database, replace separator characters, and incorporate into the list.
#Addition to the corresponding column of the information contained in the list.

    for line in data:
        if 'AC   ' in line:
            line=line.strip().replace('AC   ','').replace(';','')
            uni_id.append(line)
            table.loc[index,'Uniprot']=(", ".join(list(set(uni_id))))
            
        if 'DE   RecName: Full=' in line:
            line=line.strip().replace('DE   RecName: Full=','').replace(';','')
            full_name.append(line)
            table.loc[index,'Full Name']=(", ".join(list(set(full_name))))
            
        if 'DR   PDB;' in line:
            line=line.strip().replace('DR   ','').replace(';','')
            pdbs.append((line.split()[1]+':'+line.split()[3]))
            table.loc[index,'PDB:Resol']=(", ".join(list(set(pdbs))))

        if 'DE            EC=' in line:
            line=line.strip().replace('DE            EC=','').replace(';','')
            EC.append(line)
            table.loc[index,'EC']=(", ".join(list(set(EC))))

        if 'DR   KEGG;' in line:
            line=line.strip().replace('DR   KEGG;','').replace(';','')
            kegg.append(line)
            table.loc[index,'KEGG']=(", ".join(list(set(kegg))))

        if 'DR   Pfam' in line:
            line=line.strip().replace('DR   Pfam','').replace(';','')
            pfam.append(line)
            table.loc[index,'Pfam']=(", ".join(list(set(pfam))))
        
        if 'DR   PROSITE;' in line:
            line=line.strip().replace('DR   PROSITE;','').replace(';','')
            prosite.append(line)
            table.loc[index,'PROSITE']=(", ".join(list(set(prosite))))        

        if 'DR   GO; GO:' in line:
            line=line.strip().replace('DR   GO; GO:','').replace(';','').split(':')
            if 'F' in line[0]:
                functions.append(line[1])
                table.loc[index,'GO_function']=(", ".join(list(set(functions))))
            else:
                process.append (line[1])
                table.loc[index,'GO_process']=(", ".join(list(set(process))))

"""
4) Save full raw UniProt output
"""

import json
with open('full_data_uni.json', 'w') as f:
    f.write(json.dumps(full_data))

#Visualization of the table

table


"""
5) Output table
"""

merged = pd.merge(input, table, on='GeneID')
merged

merged.to_csv("output.csv", index=False)
