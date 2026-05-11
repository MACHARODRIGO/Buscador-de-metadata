# 🧬 Transcript Metadata Finder – Python Tool 🔎🐍

Tool for retrieving **functional and structural metadata** associated with transcript IDs from *Citrus clementina* using the **UniProt REST API**.  
Originally based on HTML scraping, the workflow was updated because UniProt no longer supports web scraping or the deprecated `uploadlists/` endpoint.

## 🎯 Learning Objectives
- 📥 Load differential expression tables using **pandas**.  
- 🔗 Query UniProtKB using **EnsemblGenome/Phytozome gene IDs** (`CICLE_v1000xxxxmg`).  
- 🧩 Parse UniProt flat-file metadata (AC, DE, GO, KEGG, Pfam, PROSITE, PDB, EC).  
- 🗂️ Merge annotations back into the input dataset.  
- 💾 Export a clean, enriched CSV.

## 🔄 Why This Update?
UniProt disabled:
- HTML-based scraping  
- JavaScript-free metadata access  
- The old `uploadlists/` mapping service  

The tool now uses the **UniProt REST API (`format=txt`)**, which returns a flat-file format compatible with the original parser.

## 📂 Input Format
A tab-delimited file containing at least:

|GeneID	| log2FC |	Padj |	Entrez_gene_ID |
|:-------|:-----------------|:----------------------|:-----------------------------------------|
|CICLE_v10000026mg|	-0.96198402954969|	0.0408352907417974|	18042740
|CICLE_v10000062mg|	0.938230044302119|	0.0146809009781383|	18040165
|CICLE_v10000069mg|	-1.53775125132245|	8.55204037181787e-6|	18042686

Only the **GeneID** column is used for UniProt queries.

## 📊 What the Script Does

1. Loads the input table with `pandas`.  
2. Extracts valid `GeneID` entries.  
3. Queries UniProt using the REST API: https://rest.uniprot.org/uniprotkb/search?query=(gene_exact:ID)&format=txt


4. Parses metadata fields:
   - UniProt accession  
   - Full protein name  
   - GO terms (Function / Process)  
   - KEGG, Pfam, PROSITE  
   - PDB structures  
   - EC numbers  

5. Builds an annotation table.  
6. Merges annotations with the original dataset.  
7. Saves the final output as `output.csv`.  
8. Stores raw UniProt responses in `full_data_uni.json`.

## 📋 Requirements

Python versions used during development:

- **Python 3.12.3**  
- **pandas 2.2.2**  
- **requests 2.32.3**

Install dependencies:


## 🖥️ How to Run

**Clone the repository:**

Default execution:
```
python metadata_finder.py
```

Using a custom input file:

```
python metadata_finder.py --csv path/to/entrada.txt
```

## 📝 Notes

- The script queries UniProt in TXT mode, preserving compatibility with the original parser.

- Querying more than ~100 genes is possible but may slow down execution.

- This version is configured for Citrus clementina gene IDs (CICLE_v1000xxxxmg), but it can be easily adapted to other genomes by adjusting the identifier format and UniProt query parameters.

