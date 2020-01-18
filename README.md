# onco_4vs2
Analysis of genes changed between immortalized and metastasizing fibroblasts. 

1. Downloaded the read counts for [SRP019968](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP019968) from recount2 website. [link](http://duffel.rail.bio/recount/v2/SRP019968/counts_gene.tsv.gz)
2. Performed DEA with DESeq2 (FDR = 0.01, lfc = 1)
3. Checked which DEGs also show a difference in protein abundance between normal tissue and cancer in the HPA.
4. Performed network enrichment with NEArender using KEGG, Reactome and WikiPathways as FGS and PathwayCommons as network.
5. Summarized the enrichment results with SUMER, performing 1 step using the AGS member genes and 2 step with all of each pathway's genes.
