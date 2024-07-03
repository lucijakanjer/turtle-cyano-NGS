
library(phyloseq)
library(vegan)
#install.packages("devtools")
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(patchwork)

library(tidyverse)

### Build Phyloseq Object 
cyano_v6 <- qza_to_phyloseq(
  "data/table-V6-CyanoSeq-cyanobacteria.qza",
  "data/rooted-tree-V6-CyanoSeq-cyanobacteria.qza",
  "data/taxonomy-CyanoSeq-V6.qza",
  "data/metadata-V6-v4.tsv") #sample metadata must be in TSV format
cyano_v6

# rarefaction
cyano_v6_rare <- rarefy_even_depth(cyano_v6,10098,rngseed = 777)
cyano_v6_rare

cyano_v6_turtle <- prune_samples(sample_data(cyano_v6)$SampleType =="carapace",cyano_v6)
cyano_v6_turtle

cyano_v6_rare_turtle <- prune_samples(sample_data(cyano_v6_rare)$SampleType =="carapace",cyano_v6_rare)
cyano_v6_rare_turtle

### Alpha diversity

alpha_diversity_v6 <- estimate_richness(cyano_v6_rare_turtle,
                                        measures=c("Observed","Shannon","InvSimpson"))
head(alpha_diversity_v6)
summary(alpha_diversity_v6)

# Prepare table for Pearson correlation test
cyano_v6_rare_turtle_meta<-data.frame(sample_data(cyano_v6_rare_turtle))
alpha_diversity_v6$Sample<-rownames(cyano_v6_rare_turtle_meta) #making column "Sample"
cyano_v6_rare_turtle_meta$Sample<-rownames(cyano_v6_rare_turtle_meta) #making column "Sample"
alpha_diversity_v6 <-left_join(alpha_diversity_v6,cyano_v6_rare_turtle_meta, by="Sample")

# Pearson correlation test
with(alpha_diversity_v6,cor.test(CCL,Shannon))
with(alpha_diversity_v6,cor.test(CCL,Observed))

# Observed features plot - Age
alpha_v6_observed_age <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  geom_smooth(method =lm, color = "grey30") +
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness", tag = "A") + 
  annotate("text", x = 55, y = 10, label = "R = 0.449, p = 0.0534")+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_observed_age

# Observed features plot - COAST
alpha_v6_observed_coast <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = coast, shape = coast), size = 3) + 
  geom_smooth(method =lm, color = "grey30") +
  scale_color_brewer(palette = "Set2", name ="coast") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness", tag = "A") + 
  annotate("text", x = 55, y = 10, label = "R = 0.449, p = 0.0534")+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_observed_coast

# Shannon diversity plot
alpha_v6_shannon_coast <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(colour = coast, shape = coast), size = 3) + 
  geom_smooth(method =lm, color = "grey30") +
  scale_shape_discrete(name  ="Adriatic coast", breaks = c("north", "south")) +
  scale_color_brewer(palette = "Set2", name ="Adriatic coast", breaks = c("north", "south")) +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index",color = "coast", shape = "coast", tag = "B") +
  annotate("text", x = 55, y = 0.25, label = "R = 0.699, p = 0.0009")+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_shannon_coast

# Combined observed richness and Shannon diversity plot
alpha_v6_observed_age + alpha_v6_shannon_age
ggsave(filename = "figures/alpha_v6_obsv_shannon_age.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_v6_obsv_shannon_age.jpg", width = 6.75, height = 4, dpi = 300)

alpha_v6_observed_coast + alpha_v6_shannon_coast
ggsave(filename = "figures/alpha_v6_obsv_shannon_coast.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_v6_obsv_shannon_coast.jpg", width = 6.75, height = 4, dpi = 300)

### Beta diversity PCA plot - For all samples

# Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu_v6 <- function(cyano_v6) {
  OTU_v6 <- otu_table(cyano_v6)
  if (taxa_are_rows(OTU_v6)) {
    OTU_v6 <- t(OTU_v6)
  }
  return(as(OTU_v6, "matrix"))
}

# CLR Transform  on unrarefied object
cyano_v6_clr <- microbiome::transform(cyano_v6, "clr")
head(otu_table(cyano_v6_clr))

# Extract Matrix and Sample Data  
cyano_v6_clr_vegan<-vegan_otu_v6(cyano_v6_clr)
cyano_v6_clr_sample<-as(sample_data(cyano_v6_clr),"data.frame")

# Principal Components Analysis
cyano_v6_pc<-prcomp(cyano_v6_clr_vegan)

# Scree plot of relative importance of explained by each axis
plot(cyano_v6_pc)

# Variance Explained by First 2 axes  
summary(cyano_v6_pc)$importance[,1:2] # 

# Extract Scores for Plotting        
cyano_v6_pc_scores<-scores(cyano_v6_pc)
cyano_v6_pc_scores_sub<-cyano_v6_pc_scores[,1:2]
#Add Sample Data  
cyano_v6_pc_scores_sub<-cbind(cyano_v6_pc_scores_sub,cyano_v6_clr_sample)

#Plot CLR PCA for all samples
pca_v6<-ggplot(cyano_v6_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level=0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", 
                       breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  scale_color_brewer(palette = "Set2", 
                     name ="Turtle Age", 
                     breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  scale_fill_brewer(palette = "Set2", 
                    name ="Turtle Age", 
                    breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  theme_bw() + 
  labs(x="PC1 (14.7%)",y="PC2 (12.8%)")
pca_v6

ggsave(filename = "figures/pca_v6.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/pca_v6.jpg", width = 6.75, height = 8, dpi = 300)

### Beta diversity PCA plot - For turtle samples

# Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu_v6_turtle <- function(cyano_v6_turtle) {
  OTU_v6_turtle <- otu_table(cyano_v6_turtle)
  if (taxa_are_rows(OTU_v6_turtle)) {
    OTU_v6_turtle <- t(OTU_v6_turtle)
  }
  return(as(OTU_v6_turtle, "matrix"))
}

# CLR Transform  on unrarefied object
cyano_v6_clr_turtle <- microbiome::transform(cyano_v6_turtle, "clr")
head(otu_table(cyano_v6_clr_turtle))

# Extract Matrix and Sample Data  
cyano_v6_clr_vegan_turtle<-vegan_otu_v6_turtle(cyano_v6_clr_turtle)
cyano_v6_clr_sample_turtle<-as(sample_data(cyano_v6_clr_turtle),"data.frame")

# Principal Components Analysis
cyano_v6_pc_turtle<-prcomp(cyano_v6_clr_vegan_turtle)

# Scree plot of relative importance of explained by each axis
plot(cyano_v6_pc_turtle)

# Variance Explained by First 2 axes  
summary(cyano_v6_pc_turtle)$importance[,1:2] # 

# Extract Scores for Plotting        
library(vegan)
cyano_v6_pc_scores_turtle<-scores(cyano_v6_pc_turtle)
cyano_v6_pc_scores_sub_turtle<-cyano_v6_pc_scores_turtle[,1:2]
#Add Sample Data  
cyano_v6_pc_scores_sub_turtle<-cbind(cyano_v6_pc_scores_sub_turtle,cyano_v6_clr_sample_turtle)

# Plot CLR PCA for turtle samples    
pca_v6_turtle<-ggplot(cyano_v6_pc_scores_sub_turtle,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level=0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", 
                       breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", 
                     breaks = c("juvenile", "sub-adult", "adult")) +
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", 
                    breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x="PC1 (17.25%)",y="PC2 (10.52%)", tag = "C")
pca_v6_turtle

# Plot CLR PCA for turtle samples - COAST  
pca_v6_turtle_coast<-ggplot(cyano_v6_pc_scores_sub_turtle,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=coast, fill=coast),level=0.95,alpha=0.2) + 
  geom_point(aes(colour = coast, shape = coast), size = 2) +
  scale_shape_discrete(name  ="Adriatic Coast", 
                       breaks = c("north", "south")) +
  scale_color_brewer(palette = "Set2", name ="Adriatic Coast", 
                     breaks = c("north", "south")) +
  scale_fill_brewer(palette = "Set2", name ="Adriatic Coast", 
                    breaks = c("north", "south")) +
  theme_bw() + 
  labs(x="PC1 (17.25%)",y="PC2 (10.52%)", tag = "C")
pca_v6_turtle_coast

ggsave(filename = "figures/pca_v6_turtle.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/pca_v6_turtle.jpg", width = 6.75, height = 8, dpi = 300)

# Combined alpha and beta diversity plot
(alpha_v6_observed_age + alpha_v6_shannon_age) / pca_v6_turtle 
ggsave(filename = "figures/alpha_pca_v6.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/alpha_pca_v6.jpg", width = 6.75, height = 8, dpi = 300)

(alpha_v6_observed_coast + alpha_v6_shannon_coast) / pca_v6_turtle_coast 
ggsave(filename = "figures/alpha_pca_v6_coast.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/alpha_pca_v6_coast.jpg", width = 6.75, height = 8, dpi = 300)
