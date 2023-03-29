# Combined V34 and V6 figures

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
library(patchwork)

#Data Manipulation
library(tidyverse)
library(MicrobeDS)

# Building phyloseq object without the phylogenetic tree
cyano_combined <- qza_to_phyloseq(
  features = "data/table-combined-cyanobacteria.qza",
  taxonomy = "data/taxonomy-combined-cyanobacteria.qza", 
  metadata = "data/metadata-combined.tsv")
cyano_combined

cyano_combined_turtle <- prune_samples(sample_data(cyano_combined)$SampleType =="carapace",cyano_combined)
cyano_combined_turtle

# Taxa bar plots
cyano_order <- cyano_combined_turtle %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional")

cyano_order %>%
  plot_composition(average_by = "Region")

ggsave(filename = "figures/taxa_bar_plot_combo_order.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_combo_order.jpg", width = 6.75, height = 4, dpi = 300)

cyano_family <- cyano_combined_turtle %>%
  aggregate_taxa(level = "Family") %>%  
  microbiome::transform(transform = "compositional")

cyano_family %>%
  plot_composition(average_by = "Region")

ggsave(filename = "figures/taxa_bar_plot_combo_family.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_combo_family.jpg", width = 6.75, height = 4, dpi = 300)


#### Run Again But With Top 'N' Phyla 

#What Are the Names of the most abundant phyla?  
cyano_combined_genus_collapse<- cyano_combined_turtle %>% aggregate_taxa(level="Genus")
cyano_combined_top12genera = names(sort(taxa_sums(cyano_combined_genus_collapse), TRUE)[1:13])

#Subset the phyloseq object to those phyla   
cyano_combined_top12genera_filter<-subset_taxa(cyano_combined_turtle,Genus %in% cyano_combined_top12genera)

#Remake Our Graph  but with grouping by CCL

cyano_combined_top12genera_plot <- cyano_combined_top12genera_filter %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(average_by = "Region")
cyano_combined_top10genera_plot + scale_fill_brewer(palette="Set3")+theme_bw()

ggsave(filename = "figures/taxa_bar_plot_combo_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_combo_genera.jpg", width = 6.75, height = 4, dpi = 300)



