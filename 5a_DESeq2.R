# 5. Analysis of Differentially Expressed Genes

# written specifically for analyzing some levels within a 1 factor design
#   i.e. only analyzing levels 'knockout' and 'heterozygous' in the factor 'Genotype', which contains 
#   3 levels: 'wildtype', 'heterozygous', 'knockout'
# for alternative code for comparing every level of 1 factor or comparing every permutation of 2 factors,
#   see workflow_wnotes.R

# 5a. Differential Expression with DESeq2 

# the design formula specifies factors to compare. 
# use it to indicate which column in colData(se) specifies experimental design
# simplest design formula: ~ Factor
# compound design formula: ~ Genotype + IP
# can look at interaction with: ~ Genotype + IP + Genotype:IP
# SummarizedExpression object created earlier can be used directly by edgeR or DESeq2

# Create DESeq2 DataSet
dds = DESeqDataSet(se.filter, design = ~ Factors)

##*Continue here from Part 3alt*##

# Normalize
dds = DESeq(dds)

# Create design matrix
contrast.deseq2 = list('Factorsknockout', 'Factorscontrol')
  # to specify we only want to compare KO vs con; se.ensembl$Factors contains PIP_KO and PIP_con as well

# Run the analysis
results.deseq2 = results(dds, contrast=contrast.deseq2, alpha=0.1)
  # can do this without specifying contrast - it will just compare every factor
  # default alpha is 0.1; if padj < alpha, then significant; set it equal to FDR 
  # default multiple test correction is BH
results.deseq2.LFC1 = results(dds, contrast=contrast.deseq2, lfcThreshold=1, alpha=0.1)
  # this only returns genes with log fold change of > 1 or < -1

# Accessing the results
summary(results.deseq2) 
  # LFC > 0 means upregulated (baseline is alphabetical) 
  # outliers are removed; low counts are removed
  # can do summary(results.deseq2, alpha=0.05)
View(results.deseq2)
  # base mean = average normalize expression across entire dataset
  # lfcSE = log fold change standard error
table(results.deseq2$padj < .1)
  # table(results.deseq2.LFC1$padj < .1)
