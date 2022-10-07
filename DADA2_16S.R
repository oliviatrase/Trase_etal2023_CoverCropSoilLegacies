### For detailed instructions, see https://benjjneb.github.io/dada2/tutorial_1_6.html ###
## This script is adapted directly from the link above ##

### 1. Initial Setup ###
### Clear workspace ###

rm(list=ls())

## make sure R knows where Rtools is located
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

### Install packages ###

## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("dada2")

### Load required libraries ###

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(data.table)
#### ASSIGN PATH ####

# CHANGE PATH to the directory containing the fastq files after unzipping.
path <- "16Sfasta" 

#### READ IN FILES ####

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### FILTERING AND TRIMMING #### 

# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter forward and reverse reads
# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F) 

#### LEARN THE ERROR RATES #### 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#### DEREPLICATION (equivalent of unique.seqs) #### 
# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### SAMPLE INFERENCE #### 
# Infer the sequence variants in each sample #
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#### MERGE PAIRED READS #### 
# Merge the denoised forward and reverse reads: #
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#### CONSTRUCT SEQUENCE TABLE #### 
#We can now construct a sequence table of our mouse samples, a #
#higher-resolution version of the OTU table produced by traditional methods.#
seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#### REMOVE CHIMERAS #### 
# Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=TRUE)

# Percent of sequences that are non-chimeric
sum(seqtab.nochim)/sum(seqtab)

#### Export and re-load to save memory ####
#  write.csv(seqtab.nochim, file = "seqtab_nochim.csv")
# seqtab.nochim<-read.csv("seqtab_nochim.csv")
# seqtab.nochim.table<-setDT(seqtab.nochim)
# rownames(seqtab.nochim.table)<-seqtab.nochim.table$X
# seqtab.nochim.table$X<-NULL
# seqtab.nochim.table<-as.matrix(seqtab.nochim.table)

#### ASSIGN TAXONOMY #### 
# Make sure to include PATH info for Silva database
taxa <- assignTaxonomy(seqtab.nochim, "silva_species_assignment_v128.fa.gz", multithread=T, verbose = T)

## Replace column names for ESV table with ESV numbers 
fin.ESV.names<-paste("ESV",sep = "",rep(1:ncol(seqtab.nochim)))
colnames(seqtab.nochim)<-fin.ESV.names

# Write ESV table to file
write.table(seqtab.nochim, file="ESV_Abund_Table.txt")

## Add ESV names to Taxonomy table and then write to file ##
taxa.mod<-cbind(fin.ESV.names,taxa)

# write taxa tavle to txt file
write.table(taxa.mod, file = "ESV_Taxonomy.txt")
