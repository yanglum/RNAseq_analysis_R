# 5. Analysis of Differentially Expressed Genes

# written specifically for analyzing some levels within a 1 factor design
#   i.e. only analyzing levels 'knockout' and 'heterozygous' in the factor 'Genotype', which contains 
#   3 levels: 'wildtype', 'heterozygous', 'knockout'
# for alternative code for comparing every level of 1 factor or comparing every permutation of 2 factors,
#   see workflow_wnotes.R

# 5b. Differential Expression with edgeR 

# Create DGE object
countdata.filter = assay(se.filter)
coldata.filter = colData(se.filter)
dge = DGEList(counts=countdata.filter, samples=coldata.filter)
# names(dge) # shows the objects contained in the DGEList

##*Continue here from Part 3alt*##

# Create design matrix
design = model.matrix(~ 0 + dge$samples$Factors)
colnames(design) = levels(dge$samples$Factors) # this changes colnames from Factorscontrol to just control

# Normalization and dispersion estimate
dge = calcNormFactors(dge) # uses TMM (upper quartile) for normalization
dge = estimateDisp(dge, design) # estimating dispersion

# GLM model fitting
fit = glmFit(dge, design)
contrast.edger = makeContrasts(knockout - control, levels=design)
  # will not accept Factorscontrol as name; only takes control 
lrt = glmLRT(fit, contrast=contrast.edger)

# Accessing results
tt = topTags(lrt, n=nrow(dge), adjust.method='BH', p.value=0.1)
  # this gives you genes meeting p.value<0.1
  # default, if no p.value is specified, is just top 10 genes (lowest p-value)
tt.all = topTags(lrt, n=nrow(dge), adjust.method='BH', sort.by='none')
results.edger = tt.all$table
  # a column here shows FDR, which is the adjusted p-value

# Set threshold for fold change:
treatres = glmTreat(fit, contrast=contrast.edger, lfc = 1)
tt.treat = topTags(treatres, n = nrow(dge), adjust.method='BH', sort.by = "none")
results.edger.LFC1 = tt.treat$table