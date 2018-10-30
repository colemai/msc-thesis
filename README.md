# msc-thesis

Master's Thesis: predicting phenotypic impact of chemical exposures on humans

1. Create n-triplets file (RDF Knowledge Graph) from CTDbase.org chemical/gene/disease data --> ctd-to-nt/
2. Create vectors of chemicals and of diseases based on the go functions of their respective positively correlated genes using Opa2Vec --> opa/opa2vec
3. Train a NN to predict from these vectors when a chemical and a disease are positively associated (marker/mechanism relationship) --> opa/opa-nn
