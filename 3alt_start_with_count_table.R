# Starting from Count Matrix Table

#*IMPORTANT*: Counts must be raw counts (pre-normalization)

# Save sample info metadata as a csv file named "samples_info_18bswitch.txt"
#   Column1 should be sample names, in same order as appear in count matrix.
#   Other columns should indicate sample info, including factors to compare by, i.e. Genotype or Treatment
# Save count matrix as a csv file named "count_matrix_ensembl.txt"
#   Row1 should be sample names, in same order as appear in sample info metadata
# The csv file must be formmated correctly
# For help converting excel file into the correct csv format, use format_excel_file.R script

# Store the sample metadata
samples.meta = read.csv('samples_info_18bswitch.txt')

# Store the table of counts (don't need to do this if used format_excel_files.R)
counts = read.csv('count_matrix.txt')
counts.mat = data.matrix(counts) # converts the data.frame to a matrix

# Create DESeq2 dataset object
dds.trans = DESeqDataSetFromMatrix(countData = counts.mat, colData = samples.meta, design = ~ Factors)
# Continue with Part 4 PCA
dds = DESeqDataSetFromMatrix(countData = counts.mat, colData = samples.meta, design = ~ Factors)
# Continue with Part 5 DESeq2

# Create DEG object (this is for edgeR)
dge = DGEList(counts = counts.mat, samples = samples.meta)
# Continue with Part 6 edgeR

################################################################################################
# Instructions for manual converstion of excel file into csv for importing to R
# Converting count matrix excel file to csv for importing:
#   Row(1) (Column names) should be sample names, in same order as sample info metadata
#   Column(1) (row names) should contain the gene IDs
#   Cell(1,1) should be empty
#   Column names and row names should be in quotation marks: "name"
#     to add quotation marks to cells, highlight cells -> Format Cells -> Custom -> enter \"@\" into Type field
#   Name file name.txt, type: CSV (comma separated values)
#     check the csv file for extra quotation marks
#     remove extra quotation marks in ms word using search and replace function
#     delete comma after the column names, if present    
#   Refer to count_matrix_ensembl.txt for correct format 