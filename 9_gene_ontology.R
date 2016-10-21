# 9. Gene set overlap analysis

# Gene Ontology - terms that describe gene function with respect to 
#   molecular function (MF)
#   biological processes (BP) they are involved in
#   cellular component (CC) of the gene's product is part of

# GO (BP) with DESeq2 results
results.deseq2.tested = results.deseq2.LFC1[ !is.na(results.deseq2.LFC1$padj), ]
genelistdown.deseq2 = factor( as.integer( results.deseq2.tested$padj < .1 
                                          & results.deseq2.tested$log2FoldChange < 0))
# now identify those genes that are downregulated, with padj < .1
results.deseq2.tested$gene.id = sapply(strsplit(row.names(results.deseq2.tested), split='.', fixed=T), 
                                       function(x) (x[1])) 
names(genelistdown.deseq2) = results.deseq2.tested$gene.id
GO.deseq2 = new( "topGOdata", ontology = 'BP',   # ontology = select the subontology (MF, BP, CC)
                 allGenes = genelistdown.deseq2,
                 nodeSize = 10, annot = annFUN.org,
                 mapping = "org.Mm.eg.db", ID = 'ensembl')

# Fisher's test (hypergeometric test) to find significant enrichment
GO.test.deseq2 = runTest(GO.deseq2, algorithm='elim', statistic = 'fisher')

GenTable(GO.deseq2, GO.test.deseq2)

# GO (BP) with edgeR results
results.edger.tested = results.edger.LFC1[ !is.na(results.edger.LFC1$FDR), ]
genelistdown.edger = factor( as.integer( results.edger.tested$FDR < .1 
                                         & results.edger.tested$logFC < 0))
results.edger.tested$gene.id = sapply(strsplit(row.names(results.edger.tested), split='.', fixed=T), 
                                       function(x) (x[1])) 
names(genelistdown.edger) = results.edger.tested$gene.id 
GO.edger = new( "topGOdata", ontology = 'BP', 
                allGenes = genelistdown.edger,
                nodeSize = 10, annot = annFUN.org,
                mapping = "org.Mm.eg.db", ID = 'ensembl')
GO.test.edger = runTest(GO.edger, algorithm='elim', statistic = 'fisher')
GenTable(GO.edger, GO.test.edger)

# Repeat for genes up, the other two sub-ontologies
# Can change log fold change (LFC) cut off and false discovery rate
#   use results.deseq2 instead of results.deseq2.LFC1 (same for edger)
#   set padj or FDR cut offs to < .2, instead of < .1