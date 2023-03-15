## Microbiome Tools
library(phyloseq)
library(vegan)
library(microbiome)

#install.packages("devtools")
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R) #for importing QIIME artifacts into phyloseq

library(DESeq2)
library(labdsv)
library(ALDEx2)

##Plotting     
#library(ggplot2)
library(RColorBrewer) #pretty colour pallettes
library(gplots)   #venn diagrams
library(cowplot) #saving graphs
library(pheatmap) #alternative heatmaps

#Data Manipulation
library(tidyverse)
library(MicrobeDS)

#Stats Models
library(lme4)
library(MuMIn)
library(car)

### Build Phyloseq Object --> Cyanobacteria without chloroplast
cyano_v34 <- qza_to_phyloseq(
  "table-V34-cyanobacteria.qza",
  "rooted-tree-V34-cyanobacteria.qza",
  "taxonomy-V34-cyanobacteria.qza",
  "metadata-v1.tsv") # sample metadata must be in TSV format
cyano_v34

# Building phyloseq object without the phylogenetic tree
cyano_v34_no_tree <- qza_to_phyloseq(
  features = "table-V34-cyanobacteria.qza",
  taxonomy = "taxonomy-V34-cyanobacteria.qza", 
  metadata = "metadata-v1.tsv")
cyano_v34_no_tree

# top 20 records from taxonomy table
head(tax_table(cyano_v34),20)

# check how many ASVs are found in different taxa on specific level
table(tax_table(cyano_v34)[,4])

# We can ask what proportion of assignments we have made at each level
apply(tax_table(cyano_v34)[,2:7],2,function(x){1-mean(is.na(x))})

# Removing non-turtle samples
cyano_v34_turtle <- subset_samples(cyano_v34, SampleType=="carapace")
