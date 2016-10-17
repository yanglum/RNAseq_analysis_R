# 0. Load Packages and Libraries

# If starting from Count Matrix Table, only need DESeq2, edgeR, and everything downstream

# Download packages (only needs to be done once - packages are stored on local drive)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
# biocLite("GenomicFeatures")
# biocLite("GenomicAlignments")
# biocLite("BiocParallel")
# biocLite("GenomeInfoDb")
# biocLite("DESeq2")
# biocLite("edgeR")
# biocLite("ggplot2")
# biocLite("pheatmap")
# biocLite("AnnotationDbi")
# biocLite("Mus.musculus")
# biocLite("topGO")

# Load libraries (needs to be done for each new R session)
library(Rsamtools) # for manipulating BAM files (sorting, indexing)
library(GenomicFeatures) # for making Transcription Database and GRangesList
library(GenomicAlignments) # the summarizeOverlaps method, which produces a SummarizedExperiment object
library(BiocParallel) # this might be necessary
library(GenomeInfoDb) # for seqlevelsStyle fxn to fix chr naming mismatch (UCSC vs NCBI/ensembl) between GRangeList and mapped reads
library(DESeq2)
library(edgeR)
library(ggplot2) # for plotting PCA
library(pheatmap)
library(AnnotationDbi) # for annotating Ensembl gene IDs 
library(Mus.musculus) # for annotating Ensembl gene IDs and for GO analysis
library(topGO) # for Gene Ontology (GO) analysis

sessionInfo() # check what libraries are loaded for this R session