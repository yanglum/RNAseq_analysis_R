# 4. Plotting Principal Component Analysis (PCA)

# PCA requires transformations, but DEG needs raw counts
# Transformations should be done in a separate workflow than DEG
dds.trans = DESeqDataSet(se.filter, design = ~ Factors) 
  # this creates a DESeq2 object from the SummarizedExperiment object 

##*Continue here from Part 3alt*##

# To solve problem of heteroscedasticity: 
# DESeq2 transforms count data that stabilize variance across the mean
# Use variance stabilizing transformation method here:
vsd = vst(dds.trans)
# class(vsd) # head(colData(vsd))
# also works on DESeqDataSetFromMatrix ojbect

# The original way of plotting PCA
plotPCA(vsd, 'Factors')
# or
plotPCA(vsd, intgroup = c('IP', 'Genotype'))

# Plotting PCA with sample names
# Run plotPCAWithSampleNames.R script
# https://gist.github.com/yanglum/00f76d65b55520cd4dbe88935bc66e4b
plotPCAWithSampleNames(vsd, 'Factors')
# or
plotPCAWithSampleNames(vsd, c('IP', 'Genotype'))

# Another way to plot PCA is using ggplot; 
#  ask plotPCA to return the data for plot, rather than actually plotting it
PCA.data = plotPCA(vsd, intgroup = c('IP', 'Genotype'), returnData=T)
percentVar = round(100*attr(PCA.data, 'percentVar'))
ggplot(PCA.data, aes(PC1, PC2, color=Genotype, shape=IP)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],'% variance')) +
  ylab(paste0("PC2: ",percentVar[2],'% variance'))

# FYI: Another way to reduce dimensionality is to use multidimensional scaling (MDS)
