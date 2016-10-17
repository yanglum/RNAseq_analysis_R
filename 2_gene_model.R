# 3. Counting Reads

# 3a. Define gene models 
# need to make Transcription Database (TxDb), and from that, GRangesList
# Three ways to make TxDb:
#   1) make TxDb from GTF file
#   2) make TxDb from UCSC directly
#   3) make TxDb from ensembl directly via BiomaRt (this is the one used in this example)

# Preferred way: 1) Make Transcription Database from GTF file
# download GTF files from Gencode, save in directory (or any directory) with BAM files as gencodeM11.gtf 
#   https://www.gencodegenes.org/mouse_releases/current.html 
#     choose comprehensive gene annotation (content), ALL (region)
gtffile = file.path(dir, 'gencodeM11.gtf')
txdb = makeTxDbFromGFF(gtffile, format='gtf') # from GenomicFeatures
ebg = exonsBy(txdb, by='gene') # produces a GRangesList of all the exons grouped by gene
# look at head(seqlevels(txdb)) or ebg to check how chromosomes are named 

# Save the TxDb file to avoid doing it everytime:
#   saveDb(ensembl, file='ensembl.sqlite')
# To load it again:
#   ensembl = loadDb('ensembl.sqlite')
