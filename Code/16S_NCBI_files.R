### 1. Initial Setup ###

### Set the working directory; modify to your own ###

### Clear workspace ###
rm(list=ls())

### Install required libraries (only do this once) ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.15")

BiocManager::install(c("phyloseq", "Heatplus"))

### Load required libraries ###
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

library(ade4)
library(vegan)
library(RColorBrewer)
library(phyloseq)
library(gplots)
library(Heatplus)
library(multcompView)
library(dplyr)
library(ggplot2)
library(randomForest)
library(caret)
library(gridExtra)
library(grid)
library(gridtext)
library(Hmisc)
library(lme4)
library(mixlm)
library(plyr)
library(data.table)

### Import Data ###
taxon <- read.table("Data/ESV_Taxonomy.txt")#, sep="\t", header=T, row.names=1)
rownames(taxon)<-taxon$fin.ESV.names
taxon$fin.ESV.names<-NULL
taxon<-taxon[!is.na(taxon$Phylum),] # remove taxa that are not defines at the Phylum level

otus <- read.table("Data/ESV_Abund_Table.txt")
rownames(otus) <- sub('.+-(.+)', '\\1', rownames(otus))

metadat <- read.csv("Data/Metadata.csv")
metadat$SampleNumber<-as.character(metadat$SampleNumber)
metadat<-metadat[order(metadat$SampleNumber),]
metadat <- subset(metadat,Type=='16S')
rownames(metadat)<-metadat$SampleNumber
metadat_rownames<-c(108)
metadat<-metadat[!(row.names(metadat) %in% metadat_rownames),]
metadat<-subset(metadat,SoilType!="ConvFallow")

### Make metadat file for NCBI ###

sample_name <- metadat$SampleName
library_ID <- sample_name
title <- '16S amplicon of soil'
library_strategy <- 'AMPLICON'
library_source <- 'GENOMIC'
library_selection <- 'PCR'
library_layout <- 'paired'
platform <- 'ILLUMINA'
instrument_model <- 'Illumina MiSeq'
design_description <- 'Amplified with 16S universal primers 515F and 806R'
filetype <- 'fastq'
filename <- list.files(pattern = "\\R1_001.fastq.gz$")
filename2 <- list.files(pattern = "\\R2_001.fastq.gz$")

NCBImeta <- data.frame(unlist(sample_name), unlist(library_ID), unlist(title), unlist(library_strategy),
                       unlist(library_source), unlist(library_selection), unlist(library_layout), unlist(platform),
                       unlist(instrument_model), unlist(design_description), unlist(filetype), unlist(filename), unlist(filename2))

names(NCBImeta) <- c("sample_name", "library_ID", "title", "library_strategy", "library_source", "library_selection", 
                     "library_layout", "platform", "instrument_model", "design_description", "filetype", "filename", "filename2")
NCBImeta$filename3<-""
NCBImeta$filename4<-""
NCBImeta$assembly<-""
NCBImeta$fasta_file<-""

write.table(NCBImeta, file = "SRA_metadata16S.txt", sep = "\t",row.names = FALSE, quote=FALSE)
write.csv(NCBImeta, file = "OWT_metadata.csv", row.names = FALSE)
