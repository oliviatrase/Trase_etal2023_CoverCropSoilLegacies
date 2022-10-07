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

### Remove problematic taxa ###

unclass <- subset(taxon, Phylum == "unclassified")
mito <- subset(taxon, Family == "Mitochondria")
arch <- subset(taxon, Kingdom == "Archaea")
chloro <- subset(taxon, Class == "Chloroplast")

un.1<-union(rownames(unclass),rownames(mito))
un.2<-union(rownames(arch),un.1)
un.3<-union(rownames(chloro),un.2)

red.taxon <- taxon[-which(rownames(taxon) %in% un.3),]

## Match ESV table with new taxonomy file ##
otus.filt<-otus[,which(colnames(otus) %in% rownames(red.taxon))]

## Determine minimum available reads per sample ##
min(rowSums(otus.filt))

## Cut samples with too few reads...cutoff here is up to you and depends on dataset ##
n2 <- names(which(rowSums(otus.filt) > 6000))
otus.re<-otus.filt[which(rownames(otus.filt) %in% n2),]

rr <- min(rowSums(otus.re))

### Rarefy to obtain even numbers of reads by sample ###
## set.seed ensures that "random" selection of sequences stays constant each time you run the script #

set.seed(336)
otus.r<-rrarefy(otus.re,rr)

metadat.r<-metadat[match(row.names(otus.r), row.names(metadat)),]
metadat.r<-metadat.r[complete.cases(metadat.r), ]
otus.r2<-otus.r[match(row.names(metadat.r), row.names(otus.r)),]
otus.r2<-otus.r2[complete.cases(otus.r2), ]

write.csv(metadat.r, file = "Data/metadata_clean.csv", row.names = FALSE)
write.csv(otus.r2, file = "Data/ESV_Abund_clean.csv", row.names = FALSE)