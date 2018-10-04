### Welcome to TurkuNLP named entity recognition and normalization systems for BioCreative VI ID assignment shared task

### Requirements
We assume that you have following softwares installed in your system in order to run our tools.
* [Simstring](http://www.chokkan.org/software/simstring/)
* [NERsuite](http://nersuite.nlplab.org/)
* [GENIA Sentence Splitter](http://www.nactem.ac.uk/y-matsu/geniass/)

### Preprocessing 
The gold-standard data from BioCreative shared task has issues in terms of word boundary, the preprocessing was used to resolve such issues. The caption input should be accompanied with full-text documents for the system to collect the correct-boundaries tokens. 

### Named entity recognition and NERsuite models
* [NERsuite models](http://bionlp-www.utu.fi/BC_VI/recognition/models)

### Normalization system
Our normalization system is based on external tools, including Simstring and Solr. We assume that you have installed those mentioned. 
* [pickle files](http://bionlp-www.utu.fi/BC_VI/normalization/pickle)
This folder is needed for both gene/protein and organism normalization systems. It contains the taxonomy tree, scientific name of organisms and lists of gene/proteins for organisms under species taxonomic rank. 
* [mapping files](http://bionlp-www.utu.fi/BC_VI/normalization/map_files)
This folder contains complementary mapping files needed for the organisms normalization systems. They include lists of model organisms, the most studied organisms according to the PubMed Central database and ranks of organisms.
* [Simstring files](http://bionlp-www.utu.fi/BC_VI/normalization/simstring)
The string matching of our normalization system relies on Simstring so we assume that you have it installed together with the python binding. The folder contained pre-compiled simstring database files and the id-symbol mapping. 
* [source data](http://bionlp-www.utu.fi/BC_VI/normalization/src_data)
* [canonical data](http://bionlp-www.utu.fi/BC_VI/normalization/data)
* [solr gene/protein data](http://bionlp-www.utu.fi/BC_VI/normalization/solr)
For gene and proteins, the mapping files are too large and too slow for mapping using the python dictionary as other entity types. So we create Solr core containing the genes/proteins in canonical form, associated taxonomy identifier, symbol type and NCBI Entrez Gene/Uniprot identifiers. This folder contains Entrez Gene and Uniprot mapping files needed for process_solr.py to add and index the entries to solr core. Prior to running the code, you need to create Solr core containing 4 data types: entrezgene_id (int), symbol (text_ws), type (int) and ncbitax_id (int). 

### Citation
If you have used data, models or parts of our systems, please kindly cite our following article.
[Suwisa Kaewphan, Kai Hakala, Niko Miekka, Tapio Salakoski, Filip Ginter; Wide-scope biomedical named entity recognition and normalization with CRFs, fuzzy matching and character level modeling, Database, Volume 2018, 1 January 2018, bay096](https://doi.org/10.1093/database/bay096)

### Authors and Contributors
Department of Future Technologies, University of Turku, Finland
* Suwisa Kaewphan
* Kai Hakala
* Niko Miekko
* Tapio Salakoski
* Filip Ginter

### Support or Contact
Please contact sukaew@utu.fi or kahaka@utu.fi for further information or questions.
