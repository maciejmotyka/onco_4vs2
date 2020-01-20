## Map Ensembl ids to gene symbols
library("org.Hs.eg.db")
my_genes_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = my_genes_ensg, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

