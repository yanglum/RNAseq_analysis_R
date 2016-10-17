# RNAseq_analysis_R

RNAseq pipeline:
Sample prep -> sequencing (.bcl) -> demultiplex (fastq) -> alignment 
(SAM/BAM) -> counting reads (count matrix) -> 
normalization and DEG analysis -> interpret results

Alignment can be done using STAR
Everything downstream can be done using R (Bioconductor package: https://www.bioconductor.org/)

GTAC provides fastq files, sorted and indexed BAM files, and count 
matrix  (excel format included)

###################################################################################################################
If count matrix not included, begin with BAM files

Go through entire R workflow:

0_load_libraries
1_sortindex_bam          *Note IGV requires indexed BAM files
2_gene_model
3_make_count_matrix
4_PCA
5a_DESeq2
5b_edgeR
etc.

###################################################################################################################
If starting from count matrix (excel file or .txt), 

go through alternative workflow:

0_load_libraries 
3alt_start_with_count_table
4_PCA 
5a_DESeq2 
5b_edgeR 
etc.
