# 7. Comparing DEG results from DESeq2 and edgeR

# make table
table(DESeq2=results.deseq2$padj<0.1, edgeR=results.edger$FDR<0.1)

# make Venn diagram
genes.deseq2 = row.names(results.deseq2)[which(results.deseq2$padj<0.1)]
genes.edger = row.names(results.edger)[which(results.edger$FDR<0.1)]
install.packages('gplots')
library(gplots)
venn(list(DESeq2=genes.deseq2, edgeR=genes.edger))

# Plot Ranks
common = !is.na(results.deseq2$padj)
plot(rank(results.deseq2$padj[common]),
     rank(results.edger$FDR[common]), cex=.1,
     xlab="DESeq2", ylab="edgR")
