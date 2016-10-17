# References:
#   http://www.bioconductor.org/help/workflows/rnaseqGene/
#   https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#preparing-count-matrices-from-bam-files

# Workflow for RNAseq
#   1. Sequencing: output = .bcl files
#   2. Demultiplexing: output = fastq files
#   3. Sequence Alignment: output = SAM and BAM files
#   4. Count reads, normalize, apply statistics for DEG
#   5. Analyze and interpret results

# The following covers steps 4 and 5:
# RNAseq analysis workflow in R (Beginning with BAM files)

# This workflow uses RNAseq results from Fgf20 KO TRAPseq on E14.5 cochleas
# reads were realigned by Bo Zhang using STAR, giving unsorted BAM files

##################################################################################################
# 0. Load Packages and Libraries

# Download packages (only needs to be done once - packages are stored on local drive) Vignettes on each package can be found on the Bioconductor website
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
# biocLite("GenomicFeatures")
# biocLite("GenomicAlignments")
# biocLite("BiocParallel")
# biocLite("GenomeInfoDb")
# biocLite("biomaRt") # necessary only for make TxDb from BiomaRt
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
library(biomaRt) # this is for the listMarts function for bioMart
library(DESeq2)
library(edgeR)
library(ggplot2) # for plotting PCA
library(pheatmap)
library(AnnotationDbi) # for annotating Ensembl gene IDs 
library(Mus.musculus) # for annotating Ensembl gene IDs and for GO analysis
library(topGO) # for Gene Ontology (GO) analysis

sessionInfo() # check what libraries are loaded for this R session

######################################################################################
#### CREATING COUNT MATRIX FROM BAM FILES ############################################
######################################################################################
# 1. Locate Path to BAM Files

# Directory "H:\Lu_Yang\Cochlea\Fgf20_E14.5_TRAPseq\18b_switch_Bo_realigned\BAM_files"
  # H: is the lab server
dir = file.path('H:', 'Lu_Yang', 'Cochlea', 'Fgf20_E14_TRAPseq', 
                '18b_switch_Bo_realigned', 'BAM_files')

file.names = dir(dir, pattern = '.bam') # *do not reassign after sorting and indexing*
  # file.name is a vector containing names of all 12 BAM files
file.dir = 1:length(file.names)
for(i in 1:length(file.names)) {file.dir[i] = file.path(dir,file.names[i])}
  # file.dir is a vector containing directories to all 12 BAM files

file.short = sapply(strsplit(file.dir, split='.', fixed=T), function(x) (x[1])) 
  # this removes .bam from names in file.dir
sorted.dir = paste0(file.short, '.sorted') # this adds .sorted to names in file.short
sorted.bam = paste0(file.short, '.sorted.bam') # this adds .sorted.bam to names in file.short
# file.exists(sorted.bam) # should return TRUE for each file

################################################################################################
# 2. Sorting and Indexing BAM files (http://biobits.org/samtools_primer.html)

# BAM files must have an index associated in order to be viewed in IGV
# Sorting takes probably 20 minutes per bam file; indexing takes probably 2 minutes per
for(i in 1:length(file.names)) {sortBam(file.dir[i], sorted.dir[i])}
  # produces .sorted.bam files
for(i in 1:length(file.names)) {indexBam(sorted.bam[i])}
  # produces .sorted.bam.bai files

############## only needed if accidenally reassigned file.names ################################
# re-storing names after sorting and indexing
# filenames = dir(dir, pattern = '.sorted.bam$') 
  # file.name is a vector containing names of all 12 BAM files
# filedir = 1:length(filenames)
# for(i in 1:length(filenames)) {filedir[i] = file.path(dir,filenames[i])}
  # filedir is a vector containing directories to all 12 .sorted.bam files
# file.exists(filedir) # should return TRUE for each file
################################################################################################

# Use Rsamtools to indicate in Bioconductor that the files in sorted.bam are BAM files
bamfiles = BamFileList(sorted.bam, yieldSize=2000000) # from Rsamtools
  # use sorted.bam files
  # yield size specifies how many reads to process at a time
  # can also use BamFile(filedir[1]) to do one BAM file at a time
# seqinfo(bamfiles[1]) to check how chromosomes are named (i.e. 1 or chr1); *take note*
  # USCS naming style: chr1, chr2, chr3 etc
  # NCBI/Ensembl naming style: 1, 2, 3 etc

################################################################################################
# 3. Counting Reads

## 3a. Define gene models ##################################################################
  # need to make Transcription Database (TxDb), and from that, GRangesList
# Three ways to make TxDb:
#   1) make TxDb from GTF file
#   2) make TxDb from UCSC directly
#   3) make TxDb from ensembl directly via BiomaRt

######### Method 1 ##########
# Make Transcription Database from GTF file
# download GTF files from Gencode
#   https://www.gencodegenes.org/mouse_releases/current.html 
#     choose comprehensive gene annotation (content), ALL (region)
#   or (not preferred) download them from UCSC or ensembl 
#     UCSC http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=538220331_gAaBXQLWqEipUmzsgHHKbGocFuoW&clade=mammal&org=Mouse&db=0&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=
#     or Ensembl ftp://ftp.ensembl.org/pub/release-86/gtf/mus_musculus
gtffile = file.path(dir, 'gencodeM11.gtf')
txdb = makeTxDbFromGFF(gtffile, format='gtf') # from GenomicFeatures
ebg = exonsBy(txdb, by='gene') # produces a GRangesList of all the exons grouped by gene
# look at head(seqlevels(txdb)) or ebg to check how chromosomes are named 

######### Method 2 ##########
# Make Transcription Database from UCSC
## supportedUCSCtables()[1:4, ]
# ucscmm10 = makeTxDbFromUCSC(genome = "mm10", tablename = "knownGene")
# getChromInfoFromUCSC("mm10") # see how chromosomes are named
  # or check seqlevels(ucscmm10)
# ebg.ucsc = exonsBy(ucscmm10, by='gene') # produces a GRangesList of all the exons grouped by gene

######### Method 3 ##########
# Make Transcription Database from Biomart (Ensembl)
# listMarts() # from BiomaRt
#  ensembl <- makeTxDbFromBiomart(dataset="mmusculus_gene_ensembl") # this is GRCm38.p4 
# getChromInfoFromBiomart(dataset='mmusculus_gene_ensembl')
# seqlevels(ensembl) # look at how chromosomes are named
# ebg.ensembl = exonsBy(ensembl, by='gene') # produces a GRangesList of all the exons grouped by gene

# if naming style of chromosomes do not match between:
#     Reads: seqinfo(bamfiles[1])
#     and
#     GRangesList: seqlevels(ensembl)
# then need to convert one to the other naming scheme
# Use seqlevelsStyle from GenomeInfoDb
  # which only works on an object you can call seqlevels() on, or a character vector
# seqlevelsStyle(ensembl)                 # gets the chr naming style
# seqlevelsStyle(ensembl) = "UCSC"        # sets the chr naming style to UCSC
# now make GRangesList again

# Save the TxDb file to avoid doing it everytime
# saveDb(txdb, file='txdb.sqlite')
# Load it again 
# txdb = loadDb('txdb.sqlite')

## 3b. Counting with summarizeOverlaps ########################################################
# run this using computer preferably with 16 GB RAM, 4-core processor (approx 2 hours)
param = SnowParam(workers = 2) # use two workers 
  # or define a BiocParallelParam with one worker
se = summarizeOverlaps(features=ebg, reads=bamfiles, mode='Union', 
                               singleEnd=T, ignore.strand=T, fragments=F, BPPARAM = param)
# this returns a *SummarizedExperiment object*, which is a count matrix (columns = samples, rows = genes)
  # mode choice see figure 1: http://bioconductor.org/packages/release/bioc/vignettes/GenomicAlignments/inst/doc/summarizeOverlaps.pdf
  # fragment = T is only for paired-end reads
    # it means also count unpaired reads in a paired-end experiment
  # unless experiment is strand-specific, set ignore.strand = T
# se # dim(se)
# access the counts using assay(se)
  # head(assay(se))
# rowRanges(se) # shows the GRangesList used to generate the se object
# str(metadata(rowRanges(se))) shows metadata about construction of the gene model
  # str() is used to display it compactly

# 3b continued. Add metadata info about samples to SummarizedExperiment object
colData(se) # stores metadata about the samples (it's empty now)
  # note it is a dataframe with nrows = # of samples
  # colnames(se) should be in same order as bamfiles
# Load the sample metadata
samples = read.csv('samples_info_18bswitch.txt') # P18b and LY18b are switched here (to correctly named)
samples.dat = read.csv('samples_info_18bswitch.txt', row.names = as.vector(samples$Sample_ID))
  # make sure rownames are sample names
colData(se) = DataFrame(samples.dat)
# can assign meta data of samples as a dataframe to colData(se) (use DataFrame() even though samples is already a data.frame)
colData(se) # now this is filled with sample metadata
# se$Factors
  # se$colname column name = columns in the colData(se)
  # you can relevel the levels in the factor se$Factors so that baseline (reference) is control
    # se$Factors = relevel(se.ensembl$Factors, 'control')

## 3c. Writing count matrix table ##################################################################
# class(assay(se)) # this is a matrix
write.table(assay(se), file='count_matrix.txt', sep=',')

## 3d. Pre-filter genes with small counts ###########################################################
se.filter = se[rowSums(assay(se)) >= 5, ]
  # this only keeps rows where sum of counts across all samples is >= 5 reads

###########################################################################################
############ RNA-SEQ ANALYSIS PIPELINES ####################################################
###########################################################################################
# 4. Plotting Principal Component Analysis (PCA)

# PCA requires transformations, but DEG needs raw counts
# Transformations should be done in a separate workflow than DEG
dds.trans = DESeqDataSet(se.filter, design = ~ Factors) # this creates a DESeq2 object
  # can transform this object for PCA plot. Create another DESeq2 object for actual DEG analysis
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

# Run plotPCAWithSampleNames.R script 
# https://gist.github.com/yanglum/00f76d65b55520cd4dbe88935bc66e4b
plotPCAWithSampleNames(vsd, 'Factors')

# Another way to plot PCA is using ggplot; 
#  ask plotPCA to return the data for plot, rather than actually plotting it
PCA.data = plotPCA(vsd, intgroup = c('IP', 'Genotype'), returnData=T)
percentVar = round(100*attr(PCA.data, 'percentVar'))
ggplot(PCA.data, aes(PC1, PC2, color=Genotype, shape=IP)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],'% variance')) +
  ylab(paste0("PC2: ",percentVar[2],'% variance'))

# FYI: Another way to reduce dimensionality is to use multidimensional scaling (MDS)

################################################################################################
# 5. Analysis of Differentially Expressed Genes

# Need to normalize counts and apply statistics; 2 commonly-used packages to do this:
#   1. DESeq2
#   2. edgeR

## 5a. Differential Expression with DESeq2 ######################################################

# the design formula specifies factors to compare. 
# use it to indicate which column in colData(se) specifies experimental design
# simplest design formula: ~ Factor
# compound design formula: ~ Genotype + IP
# can look at interaction with: ~ Genotype + IP + Genotype:IP
# SummarizedExpression object created earlier can be used directly by edgeR or DESeq2

# Create DESeq2 DataSet
dds = DESeqDataSet(se.filter, design = ~ Factors)
  # note this was stored in a different variable for transformation
dds = DESeq(dds) # normalization

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

## 5b. Differential Expression with edgeR ####################################################

# Create DGE object
countdata.filter = assay(se.filter)
coldata.filter = colData(se.filter)
# for DESeq2, can do dds.mat = DESeqDataSetFromMatrix(countData=countdata.ensfil, colData=coldata.ensfil, design = ~ Factors)
dge = DGEList(counts=countdata.filter, samples=coldata.filter)
# names(dge) # shows the objects contained in the DGEList

# Create design matrix
design = model.matrix(~ 0 + dge$samples$Factors)
colnames(design) = levels(dge$samples$Factors) # this changes colnames from Factorscontrol to just control

# Normalization
dge = calcNormFactors(dge) # uses TMM (upper quartile)
dge = estimateDisp(dge, design)

# GLM model fitting
fit = glmFit(dge, design)
contrast.edger = makeContrasts(knockout - control, levels=design)
  # will not accept Factorscontrol as name; only takes control 
lrt = glmLRT(fit, contrast=contrast.edger)

# Accessing results
res.edger = lrt$table
  # this does not include adjusted pvalue (no multiple test correction)
sig = decideTestsDGE(lrt, adjust.method="BH", p.value = 0.1) 
# View(sig) 1 means upregulated, -1 means downregulated, 0 means no change; gene order is kept the same
genes.edger = row.names(res.edger)[which(sig != 0)] # removes all genes with no change

tt = topTags(lrt, n=nrow(dge), adjust.method='BH', p.value=0.1)
  # this gives you genes meeting p.value<0.1
  # default, if no p.value is specified, is just top 10 genes
tt.all = topTags(lrt, n=nrow(dge), adjust.method='BH', sort.by='none')
  # this give syou all genes
results.edger = tt.all$table
  # a column here shows FDR, which is similar to adjusted p-value

# Set threshold for fold change:
treatres <- glmTreat(fit, contrast=contrast.edger, lfc = 1)
tt.treat <- topTags(treatres, n = nrow(dge), adjust.method='BH', sort.by = "none")
results.edger.LFC1 = tt.treat$table

#################################################################################################
# 6.Annotating and Exporting Results

# Bioconductor's annotation packages mapes gene IDs to each other
# columns(Mus.musculus) # this is a list of all available key types

# make a new column (gene.id) with ensembl id without version number
results.deseq2$gene.id = sapply(strsplit(row.names(results.deseq2), split='.', fixed=T), 
                                function(x) (x[1]))
results.edger$gene.id = sapply(strsplit(row.names(results.edger), split='.', fixed=T), 
                               function(x) (x[1]))

# mapIds function adds a column to results table with the key we want
results.deseq2$symbol = mapIds(Mus.musculus, 
                               keys=results.deseq2$gene.id, keytype='ENSEMBL', 
                               column = 'SYMBOL', multiVals='first')
#   keys indicate where the old keys are; keytype indicates the type of the old keys
#   column indicates from columns(Mus.musculus) the new keytype you want
results.deseq2$genename = mapIds(Mus.musculus, 
                                 keys=results.deseq2$gene.id, keytype='ENSEMBL', 
                                 column = 'GENENAME', multiVals='first')

results.edger$symbol = mapIds(Mus.musculus, 
                              keys=results.edger$gene.id, keytype='ENSEMBL', 
                              column = 'SYMBOL', multiVals='first')
results.edger$genename = mapIds(Mus.musculus, 
                                keys=results.edger$gene.id, keytype='ENSEMBL', 
                                column = 'GENENAME', multiVals='first')

write.csv(results.deseq2, file='deseq2_results.txt')
write.csv(results.edger, file='edger_results.txt')

###########################################################################################
############ ANALYZING AND INTERPRETING RESULTS ####################################################
###########################################################################################
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
common = !is.na(results.deseq2$padj) # this gets all padj's in results.deseq1 that is not NA
plot(rank(results.deseq2$padj[common]),
     rank(results.edger$FDR[common]), cex=.1,
     xlab="DESeq2", ylab="edgR")

##############################################################################################
# 8. Plotting Results

# Plot Counts
# takes DESeqDataSet as argument, a gene name, and factors/group to plot by
topGene = rownames(results.deseq2)[which.min(results.deseq2$padj)] 
#   this stores the gene name with the smallest padj (in deseq2 resuls)
plotCounts(dds, topGene, 'Factors')

# MA (Minus Average) plot with DESeq2
DESeq2::plotMA(results.deseq2, ylim=c(-5,5), alpha=0.1)
#   here, alpha is for thresholding padj
#   red dots are significant genes
#   DESeq2 moderates log2 fold changes from genes with very low counts and highly variable counts
#   this plot should show that only genes with large average normalized count can yield significant padj

# MA/Smear plot with edgeR
plotSmear(lrt, de.tags=rownames(tt$table))
#   de.tags defines rownames for genes identified as being differentially expressed
#   recall tt$table contains only genes with FDR < 0.1

# Heatmap of the most significant genes
no_sig_genes = sum(results.deseq2$padj<.1, na.rm=T)
# this is the number of significant genes
mat = assay(vsd)[head(order(results.deseq2$padj), no_sig_genes), ] 
#   significant genes in results.deseq2 (head function with argument for how many to show: sum(results...) evaluates to 20)
#   recall vsd is the transformed dds object
mat = mat - rowMeans(mat)
df = as.data.frame(colData(vsd)[,c('Genotype', 'IP')]) 
#   this is what to plot by; colData(vsd) = sample metadata
pheatmap(mat, annotation_col=df)

#######################################################################################
# 9. Gene set overlap analysis

# Gene Ontology - terms that describe gene function with respect to 
#   molecular function (MF)
#   biological processes (BP) they are involved in
#   cellular component (CC) of the gene's product is part of
# Mus.musculus package contains mapping of mouse genes to GO terms (org.Mm.eg.db)

# GO (BP) with DESeq2 results
results.deseq2.tested = results.deseq2.LFC1[ !is.na(results.deseq2.LFC1$padj), ]
# this keeps only those genes that have padj values, with greater than 2 fold change
genelistdown.deseq2 = factor( as.integer( results.deseq2.tested$padj < .1 
                                          & results.deseq2.tested$log2FoldChange < 0))
# now identify those genes that are down, with padj < .1
#   these will be assigned 1; rest are 0
names(genelistdown.deseq2) = rownames(results.deseq2.tested) # adds names to the factor

GO.deseq2 = new( "topGOdata", ontology = 'BP', 
                 allGenes = genelistdown.deseq2,
                 nodeSize = 10, annot = annFUN.org,
                 mapping = "org.Mm.eg.db", ID = 'ensembl')
#   nodeSize = 10 means only use sets with at least 10 genes
#   ontology = select the subontology (MF, BP, CC)

# Fisher's test (hypergeometric test) to find significant enrichment
GO.test.deseq2 = runTest(GO.deseq2, algorithm='elim', statistic = 'fisher')
#   elim means eliminate broader terms (e.g. mitochondrion) if a more narrow term 
#   (e.g. mitochondrial membrane) can be found to describe the enrichment

GenTable(GO.deseq2, GO.test.deseq2)

# GO (BP) with edgeR results
results.edger.tested = results.edger.LFC1[ !is.na(results.edger.LFC1$FDR), ]
genelistdown.edger = factor( as.integer( results.edger.tested$FDR < .1 
                                         & results.edger.tested$logFC < 0))
names(genelistdown.edger) = rownames(results.edger.tested) 
GO.edger = new( "topGOdata", ontology = 'BP', 
                allGenes = genelistdown.edger,
                nodeSize = 10, annot = annFUN.org,
                mapping = "org.Mm.eg.db", ID = 'ensembl')
GO.test.edger = runTest(GO.edger, algorithm='elim', statistic = 'fisher')
GenTable(GO.edger, GO.test.edger)

# Repeat for genes up, the other two sub-ontologies
# Can change log fold change (LFC) cut off and false discovery rate
#   use results.deseq2 instead of results.deseq2.LFC1 (same for edger)
#   set padj or FDR cut offs to < .2, instead of < .1

sessionInfo()
