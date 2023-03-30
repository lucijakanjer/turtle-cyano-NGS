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
library(ggsci)

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

order_combined <- cyano_order %>%
  plot_composition(average_by = "Region")+ 
  theme_minimal()+
  scale_fill_manual(values = safe_colorblind_palette_14)+
  labs(title="Cyanobacterial Orders", x="", y="Relative abundance")
order_combined

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


cyano_combined_top12genera_plot <- cyano_combined_top12genera_filter %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(average_by = "Region")
genera_combined <- cyano_combined_top12genera_plot+
  theme_minimal()+ 
  scale_fill_manual(values = safe_colorblind_palette_12, name="Taxa")+
  labs(title= "Top 12 Cyanobacterial Genera",x="", y="Relative abundance")
genera_combined

ggsave(filename = "figures/taxa_bar_plot_combo_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_combo_genera.jpg", width = 6.75, height = 4, dpi = 300)

order_combined / genera_combined

ggsave(filename = "figures/taxa_bar_plot_combo_order_genera.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_combo_order_genera.jpg", width = 6.75, height = 8, dpi = 300)

safe_colorblind_palette_14 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "sienna3", "olivedrab3", "#888888")

safe_colorblind_palette_12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")



