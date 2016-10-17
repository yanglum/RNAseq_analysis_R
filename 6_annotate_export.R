# 6.Annotating and Exporting Results

# columns(Mus.musculus) # this is a list of all available key types

# make a new column (gene.id) with ensembl id without version number
results.deseq2$gene.id = sapply(strsplit(row.names(results.deseq2), split='.', fixed=T), 
                                function(x) (x[1]))
results.edger$gene.id = sapply(strsplit(row.names(results.edger), split='.', fixed=T), 
                               function(x) (x[1]))

# mapIds function adds a column to results table with the key we want
results.deseq2$symbol = mapIds(Mus.musculus, 
                               keys=results.deseq2$gene.id, keytype='ENSEMBL', 
                               column = 'SYMBOL', multiVals='first')
results.deseq2$genename = mapIds(Mus.musculus, 
                                 keys=results.deseq2$gene.id, keytype='ENSEMBL', 
                                 column = 'GENENAME', multiVals='first')

results.edger$symbol = mapIds(Mus.musculus, 
                              keys=results.edger$gene.id, keytype='ENSEMBL', 
                              column = 'SYMBOL', multiVals='first')
results.edger$genename = mapIds(Mus.musculus, 
                                keys=results.edger$gene.id, keytype='ENSEMBL', 
                                column = 'GENENAME', multiVals='first')

write.csv(results.deseq2, file='deseq2_results.txt')
write.csv(results.edger, file='edger_results.txt')