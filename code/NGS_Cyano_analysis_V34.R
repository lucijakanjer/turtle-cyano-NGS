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
library(patchwork) #multiple panels

#Data Manipulation
library(tidyverse)
library(MicrobeDS)

#Stats Models
library(lme4)
library(MuMIn)
library(car)

### Build Phyloseq Object --> Cyanobacteria without chloroplast
cyano_v34 <- qza_to_phyloseq(
  "data/table-V34-cyanobacteria.qza",
  "data/rooted-tree-V34-cyanobacteria.qza",
  "data/taxonomy-V34-cyanobacteria.qza",
  "data/metadata-v1.tsv") # sample metadata must be in TSV format
cyano_v34

# Building phyloseq object without the phylogenetic tree
cyano_v34_no_tree <- qza_to_phyloseq(
  features = "data/table-V34-cyanobacteria.qza",
  taxonomy = "data/taxonomy-V34-cyanobacteria.qza", 
  metadata = "data/metadata-v1.tsv")
cyano_v34_no_tree

# top 20 records from taxonomy table
head(tax_table(cyano_v34),20)

# check how many ASVs are found in different taxa on specific level
table(tax_table(cyano_v34)[,4])

# We can ask what proportion of assignments we have made at each level
apply(tax_table(cyano_v34)[,2:7],2,function(x){1-mean(is.na(x))})

# Removing non-turtle samples
cyano_v34_turtle <- subset_samples(cyano_v34, SampleType=="carapace")

# Quality assesment

mean(sample_sums(cyano_v34)); range(sample_sums(cyano_v34))

#Make a data frame of read depths
cyano_v34_reads <- data.frame(reads=sample_sums(cyano_v34))

#Add on the sample ID
cyano_v34_reads$Sample<-rownames(cyano_v34_reads)

#Extract the Metadata from the phyloseq object 
cyano_v34_meta<-data.frame(sample_data(cyano_v34))
cyano_v34_meta$Sample<-rownames(cyano_v34_meta)

#Join on the Metadata
cyano_v34_reads<-left_join(cyano_v34_reads,cyano_v34_meta,"Sample")

#Some Boxplots of Coverage by Population using ggplot2
ggplot(cyano_v34_reads,aes(x=Age,y=reads)) + 
  geom_boxplot(aes(fill=Age)) + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plotopts

min(sample_sums(cyano_v34))
names(sample_sums(cyano_v34))[which(sample_sums(cyano_v34)<1141)]
length(sample_sums(cyano_v34))[which(sample_sums(cyano_v34)<1141)]

# Rarefiying
cyano_v34_rare <- rarefy_even_depth(cyano_v34,1141,rngseed = 777)
#`set.seed(777)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(777); .Random.seed` for the full vector
#...
#16 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
#  
#  TB159TB163TB167TB189TB191
#...
#53OTUs were removed because they are no longer 
#present in any sample after random subsampling

# Removing non-turtle samples
cyano_v34_rare_turtle <- subset_samples(cyano_v34_rare, SampleType=="carapace")
cyano_v34_rare_turtle

# Alpha diversity
alpha_diversity_v34 <- estimate_richness(cyano_v34_rare,measures=c("Observed","Shannon","InvSimpson"))
head(alpha_diversity_v34)

alpha_v34_observed_age <- plot_richness(cyano_v34_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v34_observed_age

alpha_v34_shannon_age <- plot_richness(cyano_v34_rare_turtle, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index",color = "Age", shape = "Age") +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v34_shannon_age


alpha_v34_observed_age + alpha_v34_shannon_age
ggsave(filename = "figures/alpha_v34_obsv_shannon_age.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_obsv_v34_shannon_age.jpg", 
       width = 6.75, height = 4, dpi = 300)

# preparation of data for Pearson correlation test
cyano_v34_rare_turtle_meta<-data.frame(sample_data(cyano_v34_rare_turtle))
alpha_diversity_v34$Sample<-rownames(cyano_v34_rare_turtle_meta) #making column "Sample"
cyano_v34_rare_turtle_meta$Sample<-rownames(cyano_v34_rare_turtle_meta) #making column "Sample"
alpha_diversity_v34 <-left_join(alpha_diversity_v34,cyano_v34_rare_turtle_meta, by="Sample")

# Pearson correlation test
with(alpha_diversity_v34,cor.test(CCL,Shannon))
# Pearson's product-moment correlation
#
# data:  CCL and Shannon
# t = 3.5853, df = 9, p-value = 0.005882
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3094020 0.9361287
# sample estimates:
#      cor 
# 0.7669324 

with(alpha_diversity_v34,cor.test(CCL,Observed))
# data:  CCL and Observed
# t = 5.1017, df = 9, p-value = 0.0006437
# alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.5428244 0.9636051
# sample estimates:
#  cor 
# 0.8620069 

# NMDS
# NMDS - premalo podataka za dobar prikaz !!!
ord_nmds_bray_v34 <- ordinate(cyano_v34_rare_turtle, method="NMDS",k=3, distance="bray",trymax=50)
nmds_bray_v34 <- plot_ordination(cyano_v34_rare_turtle, ord_nmds_bray_v34, color="Age", shape ="Age", title="Bray-Curtis NMDS V34") +
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  stat_ellipse(geom="polygon",type="norm",aes(colour=Age, fill=Age), alpha=0.2) + 
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw()
nmds_bray_v34
# In metaMDS(veganifyOTU(physeq), distance, ...) :
# stress is (nearly) zero: you may have insufficient data
# Too few points to calculate an ellipse

#ggsave(filename = "figures/nmds_v34_bray_age.pdf", width = 6.75, height = 4, device = cairo_pdf)
#ggsave(filename = "figures/nmds_v34_bray_age.jpg", width = 6.75, height = 4, dpi = 300)

# Combined plot of alpha and beta diversity
#(alpha_v6_observed_age + alpha_v34_shannon_age) / nmds_bray_v6
#ggsave(filename = "figures/alpha_nmds_v34.pdf", width = 6.75, height = 8, device = cairo_pdf)
#ggsave(filename = "figures/alpha_nmds_v34.jpg", width = 6.75, height = 8, dpi = 300)

# Heatmaps

cyano_v34_turtle_top20 <- prune_taxa(names(sort(taxa_sums(cyano_v34_turtle),TRUE)[1:20]), cyano_v34_turtle)
plot_heatmap(cyano_v34_turtle_top20,"NMDS",distance = "bray")

ggsave(filename = "figures/heatmap_v34_asv.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v34_asv.jpg", width = 6.75, height = 4, dpi = 300)

plot_heatmap(cyano_v34_turtle_top20,"NMDS",
             distance = "bray",
             sample.label="Age",
             sample.order = "CCL",
             taxa.label = "Genus")

ggsave(filename = "figures/heatmap_v34_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v34_genera.jpg", width = 6.75, height = 4, dpi = 300)

#Pheatmap package - Generate Metadata for Plotting Colours
age_data<-data.frame(Age=sample_data(cyano_v34_turtle_top20)$Age)
rownames(age_data)<-rownames(sample_data(cyano_v34_turtle_top20))

#Plot Heatmap 
pheatmap(otu_table(cyano_v34_turtle_top20),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale="row",
         annotation_col = age_data)

# Taxa bar plots
cyano_order_v34 <- cyano_v34_turtle %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional")

plot_composition(cyano_order_v34)

cyano_order_v34 %>%
  plot_composition(average_by = "Age")+ scale_fill_brewer(palette="Set3")

ggsave(filename = "figures/taxa_bar_plot_v34_order.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v34_order.jpg", width = 6.75, height = 4, dpi = 300)

#### Run Again But With Top 'N' Phyla 

#What Are the Names of the most abundant phyla?  
cyano_v34_genus_collapse<- cyano_v34_turtle %>% aggregate_taxa(level="Genus")
cyano_v34_top10genera = names(sort(taxa_sums(cyano_v34_genus_collapse), TRUE)[1:10])

#Subset the phyloseq object to those phyla   
cyano_v34_top10genera_filter<-subset_taxa(cyano_v34_turtle,Genus %in% cyano_v34_top10genera)

#Remake Our Graph  but with grouping by CCL

cyano_v34_top10genera_plot <- cyano_v34_top10genera_filter %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "CCL")
cyano_v34_top10genera_plot + scale_fill_brewer(palette="Set3")

ggsave(filename = "figures/taxa_bar_plot_v34_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v34_genera.jpg", width = 6.75, height = 4, dpi = 300)



### BETA DIVERISTY PLOT AND TESTS BASED ON CLR-TRANSFORMED DATA OF UNRAREFIED LIBRARY SIZES (XAH)

#### Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu_v34 <- function(cyano_v34) {
  OTU_v34 <- otu_table(cyano_v34)
  if (taxa_are_rows(OTU_v34)) {
    OTU_v34 <- t(OTU_v34)
  }
  return(as(OTU_v34, "matrix"))
}

#CLR Transform  on unrarefied object
cyano_v34_clr <- microbiome::transform(cyano_v34, "clr")
head(otu_table(cyano_v34_clr))

#Extract Matrix and Sample Data  
cyano_v34_clr_vegan<-vegan_otu_v34(cyano_v34_clr)
cyano_v34_clr_sample<-as(sample_data(cyano_v34_clr),"data.frame")

######## Principal Components Analysis
cyano_v34_pc<-prcomp(cyano_v34_clr_vegan)

#Cool Biplot Showing How Diff ASVs affect the primary axes of the ordinatiton
biplot(cyano_v34_pc)

#Scree plot of relative importance of explained by each axis
plot(cyano_v34_pc)

#Variance Explained by FIrst 2 axes  
summary(cyano_v34_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
library(vegan)
cyano_v34_pc_scores<-scores(cyano_v34_pc)
cyano_v34_pc_scores_sub<-cyano_v34_pc_scores[,1:2]
#Add Sample Data  
cyano_v34_pc_scores_sub<-cbind(cyano_v34_pc_scores_sub,cyano_v34_clr_sample)


#Housekeeping + Label Setup    
#cyano_pc_scores_sub$Species<-as.factor(cyano_pc_scores_sub$Full_names)


#Plot    
pca1_v34<-ggplot(cyano_v34_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(type="t",aes(color=Age),level = 0.95,alpha=0.5) + 
  geom_point(aes(colour=Age),size=5) 

pca2_v34<-pca1_v34 + 
  theme_bw() + 
  labs(fill="Transect Name",x="PC1 (11.7%)",y="PC2 (7.3%)") 

pca3_v34<-pca2_v34 +  
  theme(legend.position="top",
        axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text = element_text(size=12),
        legend.title = element_text(size=18))
pca3_v34

ggsave(filename = "figures/pca_clr_v34.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v34.jpg", width = 6.75, height = 8, dpi = 300)

### CLR beta diversity - only turtle samples

vegan_otu_v34_turtle <- function(cyano_v34_turtle) {
  OTU_v34_turtle <- otu_table(cyano_v34_turtle)
  if (taxa_are_rows(OTU_v34_turtle)) {
    OTU_v34_turtle <- t(OTU_v34_turtle)
  }
  return(as(OTU_v34_turtle, "matrix"))
}

#CLR Transform  on unrarefied object
cyano_v34_turtle_clr <- microbiome::transform(cyano_v34_turtle, "clr")
head(otu_table(cyano_v34_turtle_clr))

#Extract Matrix and Sample Data  
cyano_v34_turtle_clr_vegan<-vegan_otu_v34(cyano_v34_turtle_clr)
cyano_v34_turtle_clr_sample<-as(sample_data(cyano_v34_turtle_clr),"data.frame")

######## Principal Components Analysis
cyano_v34_turtle_pc<-prcomp(cyano_v34_turtle_clr_vegan)

#Cool Biplot Showing How Diff ASVs affect the primary axes of the ordinatiton
biplot(cyano_v34_turtle_pc)

#Scree plot of relative importance of explained by each axis
plot(cyano_v34_turtle_pc)

#Variance Explained by FIrst 2 axes  
summary(cyano_v34_turtle_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
library(vegan)
cyano_v34_turtle_pc_scores<-scores(cyano_v34_turtle_pc)
cyano_v34_turtle_pc_scores_sub<-cyano_v34_turtle_pc_scores[,1:2]
#Add Sample Data  
cyano_v34_turtle_pc_scores_sub<-cbind(cyano_v34_turtle_pc_scores_sub,cyano_v34_turtle_clr_sample)


#Housekeeping + Label Setup    
#cyano_pc_scores_sub$Species<-as.factor(cyano_pc_scores_sub$Full_names)


#Plot    
pca1_v34_turtle<-ggplot(cyano_v34_turtle_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(type="t",aes(color=Age),level = 0.95,alpha=0.5) + 
  geom_point(aes(colour=Age),size=5) 

pca2_v34_turtle<-pca1_v34_turtle + 
  theme_bw() + 
  labs(fill="Transect Name",x="PC1 (11.7%)",y="PC2 (7.3%)") 

pca3_v34_turtle<-pca2_v34_turtle +  
  theme(legend.position="top",
        axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text = element_text(size=12),
        legend.title = element_text(size=18))
pca3_v34_turtle

ggsave(filename = "figures/pca_clr_v34_turtles.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v34_turtles.jpg", width = 6.75, height = 8, dpi = 300)

