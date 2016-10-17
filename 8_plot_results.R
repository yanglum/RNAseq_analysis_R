# 8. Plotting Results

# Plot Counts
# takes DESeqDataSet as argument, a gene name, and factors/group to plot by
topGene = rownames(results.deseq2)[which.min(results.deseq2$padj)] 
#   this stores the gene name with the smallest padj (in deseq2 resuls)
plotCounts(dds, topGene, 'Factors')

# MA (Minus Average) plot with DESeq2
DESeq2::plotMA(results.deseq2, ylim=c(-5,5), alpha=0.1)
#   red dots are genes with padj < alpha

# MA/Smear plot with edgeR
plotSmear(lrt, de.tags=rownames(tt$table))
#   de.tags defines rownames for genes identified as being differentially expressed

# Heatmap of the most significant genes
no_sig_genes = sum(results.deseq2$padj<.1, na.rm=T) # this is the number of significant genes
mat = assay(vsd)[head(order(results.deseq2$padj), no_sig_genes), ] 
#   significant genes in results.deseq2 (head function with argument for how many to show)
mat = mat - rowMeans(mat)
df = as.data.frame(colData(vsd)[,c('Genotype', 'IP')]) 
#   this is what to plot by; colData(vsd) = sample metadata
pheatmap(mat, annotation_col=df)
