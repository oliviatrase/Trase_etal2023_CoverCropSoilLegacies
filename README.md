# Cover crop soil legacies alter maize growth, western corn rootworm performance, and microbial composition
Olivia Trase, Terrence Bell, William King, Jared Ali

Code: 20230112_Master.ipynb

Data: Master.csv and Timepoint_Distances.csv

To determine differences in plant growth, larval survival, and larval weight gain between cover crop treatments, we performed an ANOVA using a fitted linear model and a subsequent post-hoc pairwise Nemenyi test which included our block variables (date planted and maize variety) at an alpha value of 0.05. Average weight of WCR larvae was log transformed to meet ANOVA assumptions of normally distributed residuals and homogeneity. Plant growth was assessed on non-herbivore damaged treatments to isolate the influence of cover crop alone. Larval development data, or percentage of larvae which had reached the 3rd instar, were non-normally distributed and semi-continuous so differences between cover crop treatments were assessed separately for each maize variety using Kruskal-Wallis tests. 
To compare Bray-Curtis similarity from pre- to post-planting between treatments, we used Student’s t-test and ANOVA. 


Code: DADA2_16S.R and DADA2_ITS.R

Data: BioProject PRJNA858447

Raw sequence data were processed using the DADA2 pipeline (Callahan et al. 2016 p. 2). Forward and reverse sequences were trimmed to equal lengths and filtered based on overall quality and the number of errors in each read. Sequences were then dereplicated and denoised using DADA2’s parametric error model. Paired reads were then merged and chimeras removed, creating an amplicon sequence variant (ASV) table. The ASVs were assigned taxonomy using the SILVA taxonomy database (Quast et al. 2013). 


Code: 16S_Data_Cleanup.R and ITS_Data_Cleanup.R

Data: 16S_ESV_Abund_Table.txt, 16S_ESV_Taxonomy.txt, ITS_ESV_Abund_Table.txt, ITS_ESV_Taxonomy.txt

Any ASVs which were identified as Archaea, chloroplasts, or mitochondria were removed to ensure that the analysis only included the bacterial and fungal taxa targeted by the 16S rRNA gene and the ITS ribosomal region, respectively. The cleaned bacterial and fungal data were rarefied to 6000 and 5000 sequences, respectively. 


Code: 16S_Data_Analysis.R and ITS_Data_Analysis.R

Data: 16S_ESV_Abund_clean.csv, 16S_ESV_Taxonomy_clean.csv, 16S_metadata_clean.csv
      ITS_ESV_Abund_clean.csv, ITS_ESV_Taxonomy_clean.csv, ITS_metadata_clean.csv

To determine differences in composition between samples and timepoints, we calculated Bray-Curtis similarity using relative abundances. To look at the differences between time points and between treatments, we performed PERMANOVA and non-metric multidimensional scaling (NMDS). We used linear mixed effects models to determine whether each taxon was significantly positively or negatively predictive of root mass, shoot mass, WCR survival, or WCR weight gain. We used generalized linear mixed effects models to determine whether each taxon was significantly positively or negatively predictive of WCR development. For models describing plant biomass, we only used data from the first soil collection pre-planting. For models describing the WCR variables, we only used data from the second soil collection five weeks post-planting. We then controlled for multiple comparisons by using a false discovery rate adjustment.  
