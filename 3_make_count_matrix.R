# 3. Counting Reads

# 3b. Counting with summarizeOverlaps 
# run this using computer preferably with 16 GB RAM, 4-core processor (approx 2 hours)
param = SnowParam(workers = 2) # use two workers 
# or define a BiocParallelParam with one worker
se = summarizeOverlaps(features=ebg, reads=bamfiles, mode='Union', 
                       singleEnd=T, ignore.strand=T, fragments=F, BPPARAM = param)
# mode choice see figure 1: http://bioconductor.org/packages/release/bioc/vignettes/GenomicAlignments/inst/doc/summarizeOverlaps.pdf
# fragment = T is only for paired-end reads
#   it means also count unpaired reads in a paired-end experiment
# unless experiment is strand-specific, set ignore.strand = T
# se # dim(se)
# access the counts using assay(se) # head(assay(se))

# Add metadata info about samples to SummarizedExperiment object
# Load the sample metadata
samples = read.csv('samples_info_18bswitch.txt')
samples.dat = read.csv('samples_info_18bswitch.txt', row.names = as.vector(samples$Sample_ID))
# *make sure rownames are sample names*
colData(se) = DataFrame(samples.dat)
colData(se) # now this is filled with sample metadata
# se$Factors # se$colname column name = columns in the colData(se)
# you can relevel the levels in the factor se$Factors so that baseline (reference) is control
# se$Factors = relevel(se.ensembl$Factors, 'control')

# 3c. Writing count matrix table
write.table(assay(se), file='count_matrix.txt', sep=',')

# 3d. Pre-filter genes with small counts 
se.filter = se[rowSums(assay(se)) >= 5, ]
# this only keeps rows where sum of counts across all samples is >= 5 reads