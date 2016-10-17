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

# Use Rsamtools to indicate in Bioconductor that the files in sorted.bam are BAM files
bamfiles = BamFileList(sorted.bam, yieldSize=2000000) # from Rsamtools
#   use sorted.bam files
#   yield size specifies how many reads to process at a time
#   can also use BamFile(filedir[1]) to do one BAM file at a time
# seqinfo(bamfiles[1]) to check how chromosomes are named (i.e. 1 or chr1); *take note*
#   USCS naming style: chr1, chr2, chr3 etc
#   NCBI/Ensembl naming style: 1, 2, 3 etc
