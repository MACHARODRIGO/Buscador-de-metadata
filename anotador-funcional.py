"""1) Entrada"""

"""Lectura de tabla de resultados de expresion diferencial"""

import pandas as pd
input = pd.read_table('week26_entrezgene.txt')
input.head()

"""Generación de lista con id entrezgene para el parseo de UniProtKB"""

lista = input['Entrez_gene_ID'].tolist()
lista = [x for x in lista if str(x) != 'nan']
lista = [round(x) for x in lista]

""" 2) Procesamiento"""

#Acceso a UniProtKB con urllib y web scraping de pagina con BeautifulSoup (Adaptado 

#1.   Elemento de lista
#2.   Elemento de lista

#de: [https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/](https://chem-workflows.com/articles/2019/10/29/retrieve-uniprot-data-using-python/))

"""Selecciono una lista de proteinas a parsear"""

prots = lista[:25]

import urllib
from bs4 import BeautifulSoup

#Defino una funcion que busca a partir del entrezgeneid, y obtiene el accesion number de UniProtKB
def get_uniprot (query='',query_type='PDB_ID'):
    url = 'https://www.uniprot.org/uploadlists/' #Este es el webser para recuperar los datos de Uniprot
    params = {
    'from':'P_ENTREZGENEID',
    'to':'ACC',
    'format':'txt',
    'query':query
    }

    #Devuelve la salida de la pagina
    data = urllib.parse.urlencode(params)
    data = data.encode('ascii') 
    request = urllib.request.Request(url, data)
    with urllib.request.urlopen(request) as response:
        res = response.read()
        page=BeautifulSoup(res, "lxml").get_text()
        page=page.splitlines()
    return page

#Crea un data frame y una lista donde voy a guardar la salida completa
table=pd.DataFrame()
full_data = []

#Crea para cada proteina: 
for index,entry in enumerate(prots):
    #una lista con la informacion de cada base de datos
    uni_id = []
    pdbs=[]
    functions=[]
    process=[]
    full_name=[]
    kegg = []
    EC = []
    pfam = []
    prosite = []

    #busca en la base de datos 
    data=get_uniprot(query=entry)
    #apendea la salida a una lista
    full_data.append(data)
    #una fila en la table
    table.loc[index,'Entrez_gene_ID']=entry
    
#Para cada linea en la salida, busca la informacion de cada base de datos, reemplaza los caracteres separadores, y agrega esto a la lista correspondiente
#Agrega en la columna correspondiente la informacion contenida en la lista

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

"""Guardado de la salida en un archivo .json como respaldo, para posterior consulta sin necesidad de volver a acceder a UniProtKB"""

import json
with open('full_data_uni.json', 'w') as f:
    f.write(json.dumps(full_data))

"""Visualización de la tabla"""

table

""" 3) Salida  """

merged = pd.merge(input, table, on='Entrez_gene_ID')
merged["Entrez_gene_ID"] = merged["Entrez_gene_ID"].astype(int)
merged

merged.to_csv("output.csv")
