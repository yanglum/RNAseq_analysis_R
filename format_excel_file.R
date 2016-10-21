# save excel file as 'name.txt', type = csv (comma separated values)
# copy file to working directory

tabl = read.csv('name.txt', row.names=1) # row names should be the column containing ensembl ids
# if you saved the file to another directory, store the path to the file in object dir, and replace 'name.txt' with dir
#dir = file.path('H:', 'Lu_Yang', 'Cochlea', 'Fgf20_E14_TRAPseq', 'Sequencing Data', 'htcf.wustl.edu', 'files', 'ynXV2gXO', 'Ornitz_1971_6', 'all.gene_counts.txt')
#tabl = read.csv(dir, row.names=1)

############################################################################################
# if you receive error saying duplicate row names:
tab = read.csv('name.txt')
dups = which(duplicated(tab$ensembl_gene_id))
#   identifies row numbers of duplicated genes
tabl = tab[-dups,]  # subset into complement of rows dups
row.names(tabl) = tabl$ensembl_gene_id
############################################################################################

table = tabl[,9:20] # now just subset the columns with count data 
# (here I'm getting just columns 9 to 20)
counts.mat = data.matrix(table)

# now go on to downstream analyses with counts.mat as your count matrix
