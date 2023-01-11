### Clear workspace ###
rm(list=ls())

### Install required libraries (only do this once) ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.15")

BiocManager::install(c("phyloseq", "Heatplus"))

### Load required libraries ###
# tell R where Rtools is
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
library(ecodist)
library(data.table)
library(phyloseq)
library(glmmTMB)
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

# Read in Data#
metadata <- read.csv("Data/16S_metadata_clean.csv")
ESVs <- read.csv("Data/16S_ESV_Abund_clean.csv")
taxon <- read.csv("Data/16S_ESV_Taxonomy_clean.csv",row.names=1)

# Add taxonomy identifier to each taxa
taxon$Kingdom <- paste("K", taxon$Kingdom, sep="_")
taxon$Phylum <- paste("P", taxon$Phylum, sep="_")
taxon$Class <- paste("C", taxon$Class, sep="_")
taxon$Order <- paste("O", taxon$Order, sep="_")
taxon$Family <- paste("F", taxon$Family, sep="_")
taxon$Genus <- paste("G", taxon$Genus, sep="_")

# Subset data by timepoint and variety
metadata1 <- subset(metadata,Timepoint==1)
metadata2 <- subset(metadata,Timepoint==2)
metadataMC <- subset(metadata,Variety=="MastersChoice")
metadataBRO <- subset(metadata,Variety=="BlueRiverOrganic")

## Convert ESV numbers to percentages ##
ESVs_nonzero = ESVs[,colSums(ESVs) >0]
ESV.perc<-ESVs_nonzero/rowSums(ESVs_nonzero)*100

#### Difference between first and second timepoint  ####
# clean data
paired_esv <- ESVs_nonzero
paired_meta <- metadata[c("Sample","SampleName","Timepoint")]
colnames(paired_meta)<-c("subjID","sampID",'time')
rownames(paired_esv)<-paired_meta$sampID

# Ensure dataframe only has data for which there is a timepoint pair
paired_counts<-data.frame(table(paired_meta$subjID))
paired_counts_single<-paired_counts[paired_counts$Freq==1,]
single_list<-paired_counts_single$Var1
sample_list <- paired_meta[paired_meta$subjID %in% single_list, ]$sampID
paired_meta <- paired_meta[ ! paired_meta$subjID %in% single_list, ]
paired_esv <- paired_esv[ ! rownames(paired_esv) %in% sample_list, ]

# run vegdist to find difference between first and second timepoint for each plant sample
BCdist_all<-(vegdist(paired_esv, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)) 
BCdist <- melt(as.matrix(BCdist_all), varnames = c("sample1", "sample2"))
BCdist$matching<-substr(BCdist$sample1, 1, 3) == substr(BCdist$sample2, 1, 3)
BCdist<-BCdist[BCdist$matching=="TRUE",]
BCdist<-BCdist[BCdist$value!=0,]
BCdist$firsttime<-substr(BCdist$sample1, 5, 5)
BCdist<-BCdist[BCdist$firsttime=="1",]
BCdist$Sample<-substr(BCdist$sample1, 1, 3)
BCdist<-BCdist[c("Sample","value")]
# write.csv(BCdist,"Data/16S_distances.csv",row.names = FALSE)

#### 1. Run NMDS analysis ####

# Calculate NMDS
otus.nmds <- metaMDS(ESV.perc,distance="bray",k=5)

### NMDS by time point ###
plot(otus.nmds,type="n",yaxt='n',xaxt='n',main="NMDS of 16S rRNA gene ESVs \nBoth TimePoints \nBoth Varieties",ylim=c(-1,1.1),xlim=c(-1,1.3))
points(otus.nmds,display = "sites",col=c("coral3","cornflowerblue")[metadata$Timepoint],pch=16,cex=2)
axis(2,cex.axis=1.2,lwd=2,tck=-0.01)
axis(1,cex.axis=1.2,lwd=2,tck=-0.01)
text(x=1,y=1,paste("stress=",round(otus.nmds$stress,2)))
text(x=-1,y=1,"A",cex=2)

otus.pca <- prcomp(ESV.perc, center = TRUE,scale. = TRUE)
ggbiplot(otus.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=metadata$SampleName, groups=metadata$Timepoint)

### NMDS by Variety ###
metadata$Variety<-as.factor(metadata$Variety)
plot(otus.nmds,type="n",yaxt='n',xaxt='n',main="NMDS of 16S rRNA gene ESVs \nBoth TimePoints \nBoth Varieties",ylim=c(-1,1.1),xlim=c(-1,1.3))
points(otus.nmds,display = "sites",col=c("cornflowerblue","coral3")[metadata$Variety],pch=16,cex=2)
axis(2,cex.axis=1.2,lwd=2,tck=-0.01)
axis(1,cex.axis=1.2,lwd=2,tck=-0.01)
text(x=1,y=1,paste("stress=",round(otus.nmds$stress,2)))
text(x=-1,y=1,"B",cex=2)

ggbiplot(otus.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=metadata$SampleName, groups=metadata$Variety)


### NMDS by Block ###
plot(otus.nmds,type="n",yaxt='n',xaxt='n',main="NMDS of 16S rRNA gene ESVs \nBoth TimePoints \nBoth Varieties",ylim=c(-1,1.1),xlim=c(-1,1.3))
points(otus.nmds,display = "sites",col=c("black","cornflowerblue","darkorange","darkolivegreen","coral3")[metadata$Block],pch=16,cex=2)
axis(2,cex.axis=1.2,lwd=2,tck=-0.01)
axis(1,cex.axis=1.2,lwd=2,tck=-0.01)
text(x=1,y=1,paste("stress=",round(otus.nmds$stress,2)))
text(x=-1,y=1,"C",cex=2)

ggbiplot(otus.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=metadata$SampleName, groups=metadata$Block)


### NMDS Timepoint 1 ###
# ensure dataframes match by sample
ESV.perc1<-ESV.perc[match(row.names(metadata1), row.names(ESV.perc)),]
ESV.perc1 <- ESV.perc1[complete.cases(ESV.perc1), ]
ESV.perc1 = ESV.perc1[,colSums(ESV.perc1) > 0]
metadata1 <- metadata1[match(row.names(ESV.perc1), row.names(metadata1)),]

# calculate NMDS
otus.nmds1 <- metaMDS(ESV.perc1,distance="bray",k=5)

### NMDS by Cover Crop Timepoint 1 ###
metadata1$SoilType<-as.factor(metadata1$SoilType)
plot(otus.nmds1,type="n",yaxt='n',xaxt='n',main="NMDS of 16S rRNA gene ESVs \nTimepoint 1",ylim=c(-0.5,0.5),xlim=c(-0.4,0.5))
points(otus.nmds1,display = "sites",col=c("cornflowerblue","darkorange","darkolivegreen","coral3","darkorchid3")[metadata1$SoilType],pch=16,cex=2)
axis(2,cex.axis=1.2,lwd=2,tck=-0.01)
axis(1,cex.axis=1.2,lwd=2,tck=-0.01)
text(x=0.4,y=0.4,paste("stress=",round(otus.nmds1$stress,2)))
text(x=-0.4,y=0.4,"A",cex=2)

otus.pca1 <- prcomp(ESV.perc1, center = TRUE,scale. = TRUE)
ggbiplot(otus.pca1,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=metadata1$SampleName, groups=metadata1$SoilType)

### NMDS Timepoint 2 ###
# ensure dataframes match by sample
ESV.perc2<-ESV.perc[match(row.names(metadata2), row.names(ESV.perc)),]
ESV.perc2 <- ESV.perc2[complete.cases(ESV.perc2), ]
ESV.perc2 = ESV.perc2[,colSums(ESV.perc2) > 0]
metadata2 <- metadata2[match(row.names(ESV.perc2), row.names(metadata2)),]

# calculate NMDS
otus.nmds2 <- metaMDS(ESV.perc2,distance="bray",k=5)

### NMDS by Soil Type ###
metadata2$SoilType<-as.factor(metadata2$SoilType)
plot(otus.nmds2,type="n",yaxt='n',xaxt='n',main="NMDS of 16S rRNA gene ESVs \nTimepoint 2",ylim=c(-1.1,0.8),xlim=c(-1.1,1))
points(otus.nmds2,display = "sites",col=c("cornflowerblue","darkorange",
                                         "darkolivegreen","coral3","darkorchid3")[metadata2$SoilType],pch=16,cex=2)
axis(2,cex.axis=1.2,lwd=2,tck=-0.01)
axis(1,cex.axis=1.2,lwd=2,tck=-0.01)
text(x=0.6,y=0.6,paste("stress=",round(otus.nmds2$stress,2)))

otus.pca2 <- prcomp(ESV.perc2, center = TRUE,scale. = TRUE)
ggbiplot(otus.pca2,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE, labels=metadata2$SampleName, groups=metadata2$SoilType)

#### 2. 100% plots ####
### Timepoint 1 ###
# Convert esv table, metadata, and taxonomy into a combined phyloseq class
Workshop_OTU <- otu_table(ESV.perc1, taxa_are_rows=F)
Workshop_metadat <- sample_data(metadata1,errorIfNULL=TRUE)
Workshop_taxo <- tax_table(as.matrix(taxon), errorIfNULL=TRUE)
Workshop.16S <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

# Sum proportions based on class/phylum
Alpha <- subset_taxa(Workshop.16S, Class == "C_Alphaproteobacteria")
Alpha.sum<-rowSums(otu_table(Alpha))
Beta <- subset_taxa(Workshop.16S, Class == "C_Betaproteobacteria")
Beta.sum<-rowSums(otu_table(Beta))
Gamma <- subset_taxa(Workshop.16S, Class == "C_Gammaproteobacteria")
Gamma.sum<-rowSums(otu_table(Gamma))
Firmi <- subset_taxa(Workshop.16S, Phylum == "P_Firmicutes")
Firmi.sum<-rowSums(otu_table(Firmi))
Bacter <- subset_taxa(Workshop.16S, Phylum == "P_Bacteroidetes")
Bacter.sum<-rowSums(otu_table(Bacter))
Acido <- subset_taxa(Workshop.16S, Phylum == "P_Acidobacteria")
Acido.sum<-rowSums(otu_table(Acido))
Actino <- subset_taxa(Workshop.16S, Phylum == "P_Actinobacteria")
Actino.sum<-rowSums(otu_table(Actino))
Other<-100-(Alpha.sum+Beta.sum+Gamma.sum+Firmi.sum+Bacter.sum+Acido.sum+Actino.sum)

# combine new proportions
phyl.mat<-cbind(Alpha.sum,Beta.sum,Gamma.sum,Firmi.sum,Bacter.sum,Acido.sum,Actino.sum,Other)
phyl.ag<-aggregate(phyl.mat~metadata1$SoilType,FUN=mean)
rownames(phyl.ag)<-phyl.ag[,1]
phyl.ag.2<-as.matrix(phyl.ag[,-1])
phyl.ag.names<-rownames(phyl.ag)

# transpose table
phyl.ag.new<-apply(t(phyl.ag.2),2,rev)
phyl.ag.new<-melt(phyl.ag.new)

# plot %100 bar graph for each cover crop
ggplot(data=phyl.ag.new,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity")+
  labs(title="16S Community \nTimepoint 1")+
  scale_fill_manual(name = "Fungal Groups", values = c("gainsboro","midnightblue",
                                                       "mediumturquoise","mediumslateblue",
                                                       "darkmagenta","firebrick2",
                                                       "gold","forestgreen"),
                    labels = rev(c("Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria",
                                   "Firmicutes","Bacteroidetes","Acidobacteria","Actinobacteria","Other"))) + 
  theme(text = element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### Timepoint 2 ###
# Convert esv table, metadata, and taxonomy into a combined phyloseq class
Workshop_OTU <- otu_table(ESV.perc2, taxa_are_rows=F)
Workshop_metadat <- sample_data(metadata2,errorIfNULL=TRUE)
Workshop_taxo <- tax_table(as.matrix(taxon), errorIfNULL=TRUE)
Workshop.16S <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

# Sum proportions based on class/phylum
Alpha <- subset_taxa(Workshop.16S, Class == "C_Alphaproteobacteria")
Alpha.sum<-rowSums(otu_table(Alpha))
Beta <- subset_taxa(Workshop.16S, Class == "C_Betaproteobacteria")
Beta.sum<-rowSums(otu_table(Beta))
Gamma <- subset_taxa(Workshop.16S, Class == "C_Gammaproteobacteria")
Gamma.sum<-rowSums(otu_table(Gamma))
Firmi <- subset_taxa(Workshop.16S, Phylum == "P_Firmicutes")
Firmi.sum<-rowSums(otu_table(Firmi))
Bacter <- subset_taxa(Workshop.16S, Phylum == "P_Bacteroidetes")
Bacter.sum<-rowSums(otu_table(Bacter))
Acido <- subset_taxa(Workshop.16S, Phylum == "P_Acidobacteria")
Acido.sum<-rowSums(otu_table(Acido))
Actino <- subset_taxa(Workshop.16S, Phylum == "P_Actinobacteria")
Actino.sum<-rowSums(otu_table(Actino))
Other<-100-(Alpha.sum+Beta.sum+Gamma.sum+Firmi.sum+Bacter.sum+Acido.sum+Actino.sum)

# combine new proportions
phyl.mat<-cbind(Alpha.sum,Beta.sum,Gamma.sum,Firmi.sum,Bacter.sum,Acido.sum,Actino.sum,Other)
phyl.ag<-aggregate(phyl.mat~metadata2$SoilType,FUN=mean)
rownames(phyl.ag)<-phyl.ag[,1]
phyl.ag.2<-as.matrix(phyl.ag[,-1])
phyl.ag.names<-rownames(phyl.ag)

# transpose table
phyl.ag.new<-apply(t(phyl.ag.2),2,rev)
phyl.ag.new<-melt(phyl.ag.new)

# plot %100 bar graph for each cover crop
ggplot(data=phyl.ag.new,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity")+
  labs(title="16S Community \nTimepoint 2")+
  scale_fill_manual(name = "Fungal Groups", values = c("gainsboro","midnightblue",
                                                       "mediumturquoise","mediumslateblue",
                                                       "darkmagenta","firebrick2",
                                                       "gold","forestgreen"),
                    labels = rev(c("Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria",
                                   "Firmicutes","Bacteroidetes","Acidobacteria","Actinobacteria","Other"))) + 
  theme(text = element_text(size=16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#### How presence of bacterial taxa affect corn growth/insect performance ####
# read in data, clean
otherdf<- read.csv("Data/Master.csv")
rownames(otherdf)<-otherdf$Sample
otherdf$Sample<-NULL



### Mixed models for Plant Data ### 
# Use first timepoint for plant data

# clean data, make sure samples match
rownames(ESV.perc1)<-metadata1$Sample
otherdf.t1<-otherdf[match(row.names(ESV.perc1), row.names(otherdf)),]

# transpose data
t.ESV.perc1<-t(ESV.perc1)
t.ESV.perc1<-data.frame(t.ESV.perc1)

# merge ESV abundance and taxonomy
taxonesv<-merge(taxon,t.ESV.perc1,by=0,all=TRUE)
rownames(taxonesv)<-taxonesv$Row.names
taxonesv$Row.names<-NULL

# Get proportional abundance for each genus, family, order, class, phylum separately
taxongroup.df.list <- list()
num<-length(colnames(taxon))
for (i in 2:num){
  taxongroup <- colnames(taxon)[i]
  tdf<-data.frame(rownames(taxon),taxon[,taxongroup])
  rownames(tdf)<-tdf$rownames.taxon.
  tdf$rownames.taxon.<-NULL
  taxonesv<-merge(tdf,t.ESV.perc1,by=0,all=TRUE)
  rownames(taxonesv)<-taxonesv$Row.names
  taxonesv$Row.names<-NULL
  taxongroup.esv<-taxonesv[!is.na(taxonesv$taxon...taxongroup.),]
  taxongroup.agg <- aggregate(. ~ taxon...taxongroup., taxongroup.esv, function(x) c(sum(x)))
  rownames(taxongroup.agg)<-taxongroup.agg$taxon...taxongroup.
  taxongroup.agg$taxon...taxongroup.<-NULL
  t.taxongroup.agg<-t(taxongroup.agg)
  dfnum <- i-1
  taxongroup.df.list[[dfnum]]<-data.frame(t.taxongroup.agg)
}
phylum<-taxongroup.df.list[[1]]
class<-taxongroup.df.list[[2]]
order<-taxongroup.df.list[[3]]
family<-taxongroup.df.list[[4]]
genus<-taxongroup.df.list[[5]]
# species<-taxongroup.df.list[[6]] # no bacterial species data
taxonlist<-cbind(phylum,class,order,family,genus)

# merge ESV abundances with genus, family, order, class, phylum data
esv.all.perc <- merge(ESV.perc1,taxonlist,by=0,all=TRUE)
rownames(esv.all.perc)<-esv.all.perc$Row.names
esv.all.perc$Row.names<-NULL

# Make esv table binary, presence/absence 
esv.all.perc.bi<-esv.all.perc
esv.all.perc.bi[esv.all.perc.bi>0]<-1

# # remove taxa for which the standard deviation is zero (0% have the microbe or 100% have the microbe)
# esv.all.perc.bi<-esv.all.perc.bi[, !sapply(esv.all.perc.bi, function(x) { sd(if (is.factor(x)) as.integer(x) else x) == 0} )]
# 
# # get the number of taxa left
# num <- length(colnames(esv.all.perc.bi))

# Add in corn and insect data as well as treatment and block effects
esv.all.perc.bi$RootMass <- otherdf.t1$RootMass
esv.all.perc.bi$ShootMass <- otherdf.t1$ShootMass
esv.all.perc.bi$Block <- as.factor(otherdf.t1$DateGerm)
esv.all.perc.bi$Variety <- as.factor(otherdf.t1$Variety)
esv.all.perc.bi$SoilType <- as.factor(otherdf.t1$SoilType)

# remove samples that were used in insect assay
esv.all.perc.bi<-esv.all.perc.bi[!is.na(esv.all.perc.bi$RootMass),] 

# remove taxa for which the standard deviation is zero (0% have the microbe or 100% have the microbe)
esv.all.perc.bi<-Filter(function(x) sd(if (is.factor(x)) as.integer(x) else x,na.rm=T) != 0, esv.all.perc.bi) 

# Loop through taxa and run mixed model for shoot and root mass:
# Is shoot mass significantly greater in the presence of microbe x 
# considering fixed effects cover crop and corn variety and random effect block
num <- length(colnames(esv.all.perc.bi))-5
esvlist<-list()
ratios<-list()
rootp <- list()
shootp <- list()
for (i in 1:num){
  esv <- colnames(esv.all.perc.bi)[i]
  
  esvcol<-esv.all.perc.bi[,esv]
  ratio <- sum(esvcol==0)/sum(esvcol==1) # get ratio to ensure there aren't only one or two samples that have the microbe or don't
  ratios[i]<-ratio
  
  root <- esv.all.perc.bi$RootMass
  shoot <- esv.all.perc.bi$ShootMass
  blocks <-esv.all.perc.bi$Block
  variety <- esv.all.perc.bi$Variety
  soils <- esv.all.perc.bi$SoilType
  
  esvlist[i]<-esv
  
  model <- lmer(root ~ esvcol*variety*soils + (1 | blocks))
  aov<-Anova(model)
  pval <- aov[1,3]
  rootp[i]<-pval
  
  model <- lmer(shoot ~ esvcol*variety*soils + (1 | blocks))
  aov<-Anova(model)
  pval <- aov[1,3]
  shootp[i]<-pval
  
}

# Combine the ESV, the ratio of samples that did not have the microbe to those that did, 
# the pvalue for whether root mass was greater in the presence of the microbe, and
# the pvalue for whether shoot mass was greater in the presence of the microbe
anovadf.t1.corn.anova <- do.call(rbind, Map(data.frame, ESV=esvlist,RootP=rootp, ShootP=shootp,Ratio=ratios))

# limit values by ensuring that there are somewhat equal numbers between samples that do have x microbe
# and those that do not, helps for exluding values where only one or two samples have microbe
anovadf.t1.corn.anova.short<-subset(anovadf.t1.corn.anova, Ratio>0.39 & Ratio<2.51) 

# make separate lists for root and shoot data
rootesvs<-subset(anovadf.t1.corn.anova.short,select=c("ESV","RootP"))
shootesvs<-subset(anovadf.t1.corn.anova.short,select=c("ESV","ShootP"))
implist<-list(rootesvs,shootesvs)
ylist<-c("RootMass","ShootMass")

# loop through roots and shoots dataframes to get more detailed info about differences
plotdata<-list()
for (i in 1:2){ 
  plotdata2<-list()
  
  # dataframe of taxa and p values
  esvs <- data.frame(implist[i])
  
  # label Rootmass or Shootmass
  ycol<- ylist[i] 
  
  # subset ESV proportional abundance table by variable of interest
  df<-esv.all.perc.bi[,which(colnames(esv.all.perc.bi) %in% esvs$ESV)]
  
  # get number of microbial taxa to loop through
  num <- length(colnames(df))
  
  # add in rootmass or shootmass to abundance table
  df$ycol <- esv.all.perc.bi[,ycol] 
  
  # add in cover crop treatment to abundance table
  df$SoilType<-esv.all.perc.bi$SoilType
 
  for(j in 1:num){ # loop through ESVs
    # pvalue for ESV i
    p<-esvs[j,2] 
    
    # name of ESV or microbial taxa
    taxa <- colnames(df)[j] 
   
    # dataframe with root/shootmass values in one column and 0/1 in other column for absence/presence of microbe j
    taxadf <- data.frame(df$ycol,df[,taxa]) 
    colnames(taxadf)<-c("ycol","taxaval")
    
    # change 0 to absent and 1 to present
    taxadf$factor <- 0
    taxadf$factor[taxadf$taxaval==0]<-"absent"
    taxadf$factor[taxadf$taxaval==1]<-"present"
    factorlist <- list(taxadf$factor)
    taxadf$taxaval<-NULL
    taxadf$factor<-NULL
    
    # get mean, se, and counts for difference in shoot/root mass when microbe j is present
    taxadf.agg<-aggregate(taxadf,by=factorlist,
                          FUN = function(x) c(mean = mean(x),se = sd(x)/sqrt(length(x)), count=length(x)))
    
    # change aggregate table to dataframe
    taxadf.agg2<-data.frame(taxadf.agg$ycol)
    
    # add presence/absence labels
    taxadf.agg2$group <- c("absent","present")
    
    diff <- taxadf.agg2$mean[2]-taxadf.agg2$mean[1]
    
    # shoot/root mass when the microbe is present subtracted by shoot/root mass when microbe is absent
    taxadf.agg2$diff[diff<0] <- "worse"
    diff <- taxadf.agg2$mean[2]-taxadf.agg2$mean[1] 
    taxadf.agg2$diff[diff>0] <- "better"
    
    # add column name of microbe
    taxadf.agg2$taxa<-taxa
    
    # add column for taxonomic group of microbe
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"P")]<-"Phylum"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"G")]<-"Genus"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"F")]<-"Family"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"O")]<-"Order"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"C")]<-"Class"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"E")]<-"ESV"
    
    # add column for root or shoot mass
    taxadf.agg2$metric<-ycol
    
    # add column for pvalue
    taxadf.agg2$pval<-p
    
    plotdata2[[j]]<-taxadf.agg2      
  }
  # combine all microbe data for shoot mass then root mass
  plotdata[[i]]<-do.call("rbind",plotdata2)
}
# add shoot and root mass together
plantplotdata<-do.call("rbind",plotdata)



### Mixed models for WCR Data ### 
# Use second timepoint for wcr data

# clean data, make sure samples match
rownames(ESV.perc2)<-metadata2$Sample
otherdf.t2<-otherdf[match(row.names(ESV.perc2), row.names(otherdf)),]

# transpose data
t.ESV.perc2<-t(ESV.perc2)
t.ESV.perc2<-data.frame(t.ESV.perc2)

# merge ESV abundance and taxonomy
taxonesv2<-merge(taxon,t.ESV.perc2,by=0,all=TRUE)
rownames(taxonesv2)<-taxonesv2$Row.names
taxonesv2$Row.names<-NULL

# Get proportional abundance for each genus, family, order, class, phylum separately
taxongroup.df.list2 <- list()
num<-length(colnames(taxon))
for (i in 2:num){
  taxongroup <- colnames(taxon)[i]
  tdf<-data.frame(rownames(taxon),taxon[,taxongroup])
  rownames(tdf)<-tdf$rownames.taxon.
  tdf$rownames.taxon.<-NULL
  taxonesv2<-merge(tdf,t.ESV.perc2,by=0,all=TRUE)
  rownames(taxonesv2)<-taxonesv2$Row.names
  taxonesv2$Row.names<-NULL
  taxongroup.esv<-taxonesv2[!is.na(taxonesv2$taxon...taxongroup.),]
  taxongroup.agg <- aggregate(. ~ taxon...taxongroup., taxongroup.esv, function(x) c(sum(x)))
  rownames(taxongroup.agg)<-taxongroup.agg$taxon...taxongroup.
  taxongroup.agg$taxon...taxongroup.<-NULL
  t.taxongroup.agg<-t(taxongroup.agg)
  dfnum <- i-1
  taxongroup.df.list2[[dfnum]]<-data.frame(t.taxongroup.agg)
}
phylum<-taxongroup.df.list2[[1]]
class<-taxongroup.df.list2[[2]]
order<-taxongroup.df.list2[[3]]
family<-taxongroup.df.list2[[4]]
genus<-taxongroup.df.list2[[5]]
taxonlist2<-cbind(phylum,class,order,family,genus)

# merge ESV abundances with genus, family, order, class, phylum data
esv.all.perc2 <- merge(ESV.perc2,taxonlist2,by=0,all=TRUE)
rownames(esv.all.perc2)<-esv.all.perc2$Row.names
esv.all.perc2$Row.names<-NULL

# Make esv table binary, presence/absence 
esv.all.perc2.bi<-esv.all.perc2
esv.all.perc2.bi[esv.all.perc2.bi>0]<-1

# Add in insect data as well as treatment and block effects
esv.all.perc2.bi$WCRSurvival <- otherdf.t2$WCRSurvival
esv.all.perc2.bi$WCRWeightAvg <- otherdf.t2$WCRWeightAvg
esv.all.perc2.bi$WCR3rdInstar <- otherdf.t2$WCR3rdInstar
esv.all.perc2.bi$Block <- as.factor(otherdf.t2$DateGerm)
esv.all.perc2.bi$Variety <- as.factor(otherdf.t2$Variety)
esv.all.perc2.bi$SoilType <- as.factor(otherdf.t2$SoilType)
esv.all.perc2.bi$Count <- otherdf.t2$WCRRecovered

# remove samples that don't have insect data
esv.all.perc2.bi<-esv.all.perc2.bi[!is.na(esv.all.perc2.bi$WCRWeightAvg),]

# remove taxa for which the standard deviation is zero (0% have the microbe or 100% have the microbe)
esv.all.perc2.bi<-Filter(function(x) sd(if (is.factor(x)) as.integer(x) else x,na.rm=T) != 0, esv.all.perc2.bi)

# Loop through taxa and run mixed model for insect performance:
# Is WCR survival significantly greater in the presence of microbe x 
# considering fixed effects cover crop and corn variety and random effect block
num <- length(colnames(esv.all.perc2.bi))-7
esvlist<-list()
ratios<-list()
survp <- list()
weightp <- list()
devp <- list()

for (i in 1:num){
  esv <- colnames(esv.all.perc2.bi)[i]
  esvlist[i]<-esv
  
  esvcol<-esv.all.perc2.bi[,esv]
  ratio <- sum(esvcol==0)/sum(esvcol==1)
  ratios[i]<-ratio
  
  surv <- esv.all.perc2.bi$WCRSurvival
  weight <- log(esv.all.perc2.bi$WCRWeightAvg) # use log to meet lmer assumptions
  dev <- esv.all.perc2.bi$WCR3rdInstar
  
  blocks <-esv.all.perc2.bi$Block
  variety <- esv.all.perc2.bi$Variety
  soils <- esv.all.perc2.bi$SoilType
  counts <- esv.all.perc2.bi$Count
  
  # WCR survival model
  model <- lmer(surv ~ esvcol*variety*soils + (1 | blocks))
  aov<-Anova(model)
  pval <- aov[1,3]
  survp[i]<-pval
  
  # check assumptions
  # plot<-plot(model)
  # print(plot)
  # qqnorm(resid(model))
  
  # WCR weight gain model, log transformed
  model <- lmer(weight ~ esvcol*variety*soils + (1 | blocks))
  aov<-Anova(model)
  pval <- aov[1,3]
  weightp[i]<-pval
  
  # check assumptions
  # plot<-plot(model)
  # print(plot)
  # qqnorm(resid(model))
  
  # WCR development model, binomial (non normaa and heteroscedastic)
  model <- glmer(dev~esvcol*variety*soils+(1|blocks),weights=counts,family="binomial")
  aov<-Anova(model)
  pval<-aov[1,3]
  devp[i]<-pval
  
  # check assumptions
  # plot<-plot(model)
  # print(plot)
  # qqnorm(resid(model))
}
# Combine the ESV, the ratio of samples that did not have the microbe to those that did, 
# the pvalue for whether root mass was greater in the presence of the microbe, and
# the pvalue for whether shoot mass was greater in the presence of the microbe
anovadf.t2.wcr.anova <- do.call(rbind, Map(data.frame, ESV=esvlist,SurvP=survp, WeightP=weightp, DevP=devp,Ratio=ratios))

# limit values by ensuring that there are somewhat equal numbers between samples that do have x microbe
# and those that do not, helps for exluding values where only one or two samples have microbe
anovadf.t2.wcr.anova.short<-subset(anovadf.t2.wcr.anova, Ratio>0.39 & Ratio<2.51)

# make separate lists for root and shoot data
survesvs<-subset(anovadf.t2.wcr.anova.short,select=c("ESV","SurvP"))
weightesvs<-subset(anovadf.t2.wcr.anova.short,select=c("ESV","WeightP"))
devesvs<-subset(anovadf.t2.wcr.anova.short,select=c("ESV","DevP"))
implist<-list(survesvs,weightesvs,devesvs)
ylist<-c("WCRSurvival","WCRWeightAvg","WCR3rdInstar")

# loop through roots and shoots dataframes to get more detailed info about differences
plotdata<-list()
for (i in 1:3){ 
  plotdata2<-list()
  
  # dataframe of taxa and p values
  esvs <- data.frame(implist[i])
  
  # label Rootmass or Shootmass
  ycol<- ylist[i] 
  
  # subset ESV proportional abundance table by variable of interest
  df<-esv.all.perc2.bi[,which(colnames(esv.all.perc2.bi) %in% esvs$ESV)]
  
  # get number of microbial taxa to loop through
  num <- length(colnames(df))
  
  # add in rootmass or shootmass to abundance table
  df$ycol <- esv.all.perc2.bi[,ycol] 
  
  # add in cover crop treatment to abundance table
  df$SoilType<-esv.all.perc2.bi$SoilType
  
  for(j in 1:num){ # loop through ESVs
    # pvalue for ESV i
    p<-esvs[j,2] 
    
    # name of ESV or microbial taxa
    taxa <- colnames(df)[j] 
    
    # dataframe with root/shootmass values in one column and 0/1 in other column for absence/presence of microbe j
    taxadf <- data.frame(df$ycol,df[,taxa]) 
    colnames(taxadf)<-c("ycol","taxaval")
    
    # change 0 to absent and 1 to present
    taxadf$factor <- 0
    taxadf$factor[taxadf$taxaval==0]<-"absent"
    taxadf$factor[taxadf$taxaval==1]<-"present"
    factorlist <- list(taxadf$factor)
    taxadf$taxaval<-NULL
    taxadf$factor<-NULL
    
    # get mean, se, and counts for difference in shoot/root mass when microbe j is present
    taxadf.agg<-aggregate(taxadf,by=factorlist,
                          FUN = function(x) c(mean = mean(x),se = sd(x)/sqrt(length(x)), count=length(x)))
    
    # change aggregate table to dataframe
    taxadf.agg2<-data.frame(taxadf.agg$ycol)
    
    # add presence/absence labels
    taxadf.agg2$group <- c("absent","present")
    diff <- taxadf.agg2$mean[2]-taxadf.agg2$mean[1]
    # shoot/root mass when the microbe is present subtracted by shoot/root mass when microbe is absent
    taxadf.agg2$diff[diff<0] <- "worse"
    diff <- taxadf.agg2$mean[2]-taxadf.agg2$mean[1] 
    taxadf.agg2$diff[diff>0] <- "better"
    
    # add column name of microbe
    taxadf.agg2$taxa<-taxa
    
    # add column for taxonomic group of microbe
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"P")]<-"Phylum"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"G")]<-"Genus"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"F")]<-"Family"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"O")]<-"Order"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"C")]<-"Class"
    taxadf.agg2$taxagroup[startsWith(taxadf.agg2$taxa,"E")]<-"ESV"
    
    # add column for root or shoot mass
    taxadf.agg2$metric<-ycol
    
    # add column for pvalue
    taxadf.agg2$pval<-p
    
    plotdata2[[j]]<-taxadf.agg2      
  }
  # combine all microbe data for shoot mass then root mass
  plotdata[[i]]<-do.call("rbind",plotdata2)
}
# add shoot and root mass together
wcrplotdata<-do.call("rbind",plotdata)

# Combine plant data and WCR data regarding absence/presence of microbes
allplotdata.list<-(list(plantplotdata,wcrplotdata))
allplotdata<-do.call("rbind",allplotdata.list)

# create list for the taxa that, when present, indicate a significant difference in at least one variable
imptaxa<-unique(subset(allplotdata,pval<=0.05,select=c("taxa")))
allplotdata.sig <- allplotdata[allplotdata$taxa %in% imptaxa$taxa,]
allplotdata.present <- subset(allplotdata.sig,group=="present")
rownames(allplotdata.present)<-NULL

# Remove specific ESv for now
allplotdata.short <- subset(allplotdata.present,taxagroup!='ESV')

# create gradient for pvalues
allplotdata.short$level<-NA
allplotdata.short$level[allplotdata.short$pval>=0.1] <- 0.4
allplotdata.short$level[allplotdata.short$pval<0.1] <- 0.5
allplotdata.short$level[allplotdata.short$pval<0.05] <- .7
allplotdata.short$level[allplotdata.short$pval<0.01] <- .85
allplotdata.short$level[allplotdata.short$pval<0.001] <- 1

# Reordering group factor levels
allplotdata.short$taxagroup <- factor(allplotdata.short$taxagroup,     
                                      levels = c("Class", "Order","Family" ,"Genus"))
# Reordering metric levels
allplotdata.short$metric <-factor(allplotdata.short$metric,
                                    levels = c("ShootMass","RootMass","WCRSurvival","WCRWeightAvg","WCR3rdInstar"))
allplotdata.short$taxa <- substring(allplotdata.short$taxa, 3)
metrics<-c("Shoot \nMass", "Root \nMass", "WCR \nSurvival", "WCR \nWeight", "WCR \nDevelopment")

# Cleanup final graphic
allplotdata.short<-allplotdata.short[- grep("Subgroup", allplotdata.short$taxa),]
allplotdata.short<-allplotdata.short[- grep("[0-9]", allplotdata.short$taxa),]
# Remove duplicate taxa
to_remove<-c("Flavobacteriales","Clostridiales","Pseudonocardiaceae","Caulobacteraceae","Iamia")
allplotdata.short<-allplotdata.short[!(allplotdata.short$taxa %in% to_remove),]

# Plot heatmap
ggplot(allplotdata.short, aes(metric, taxa, fill= diff,alpha=-pval)) + 
  xlab(NULL)+ylab(NULL)+
  geom_tile()+
  facet_grid(taxagroup~., scales = "free",space="free")+
  scale_fill_manual(values=c("darkgreen","mediumpurple4"))+
  scale_x_discrete(labels= metrics,position="top")+
  theme(text = element_text(size=15,family="serif"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_text(size=14,family="serif",face="italic"))
