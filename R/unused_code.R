## Map Ensembl ids to gene symbols
library("org.Hs.eg.db")
my_genes_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = my_genes_ensg, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

## Verify the gene symbols using HGNC Multi-symbol checker
# 1. Copy the genes to clipboard
clipr::write_clip(not_mapped)
# 2. Paste into https://www.genenames.org/tools/multi-symbol-checker
# 3. Select "Approved symbol" and "Previous symbol"
# 4. Copy the returned approved symbols
# 5. Read from clipboard
not_mapped_approved <- clipr::read_clip() %>% unique()
