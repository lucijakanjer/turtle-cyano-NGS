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
bacteria_archea_v34 <- qza_to_phyloseq(
  "data/table-V34.qza",
  "data/rooted-tree-V34.qza",
  "data/taxonomy-V34.qza",
  "data/metadata-v1.tsv") # sample metadata must be in TSV format
bacteria_archea_v34

# Building phyloseq object without the phylogenetic tree
bacteria_archea_v34_no_tree <- qza_to_phyloseq(
  features = "data/table-V34.qza",
  taxonomy = "data/taxonomy-V34.qza", 
  metadata = "data/metadata-v1.tsv")
bacteria_archea_v34_no_tree

# top 20 records from taxonomy table
head(tax_table(bacteria_archea_v34),20)

# check how many ASVs are found in different taxa on specific level
table(tax_table(bacteria_archea_v34)[,1])

bacteria_v34<-prune_taxa(as.logical(tax_table(bacteria_archea_v34)[,1]=="d__Bacteria"),bacteria_archea_v34)
bacteria_v34

# We can ask what proportion of assignments we have made at each level
apply(tax_table(bacteria_v34)[,2:7],2,function(x){1-mean(is.na(x))})

# Removing non-turtle samples
bacteria_v34_turtle <- subset_samples(bacteria_v34, SampleType=="carapace")

# Quality assesment

mean(sample_sums(bacteria_v34)); range(sample_sums(bacteria_v34))

#Make a data frame of read depths
bacteria_v34_reads <- data.frame(reads=sample_sums(bacteria_v34))

#Add on the sample ID
bacteria_v34_reads$Sample<-rownames(bacteria_v34_reads)

#Extract the Metadata from the phyloseq object 
bacteria_v34_meta<-data.frame(sample_data(bacteria_v34))
bacteria_v34_meta$Sample<-rownames(bacteria_v34_meta)

#Join on the Metadata
bacteria_v34_reads<-left_join(bacteria_v34_reads,bacteria_v34_meta,"Sample")

#Some Boxplots of Coverage by Population using ggplot2
ggplot(bacteria_v34_reads,aes(x=Age,y=reads)) + 
  geom_boxplot(aes(fill=Age)) + 
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
#plotopts

min(sample_sums(bacteria_v34))
names(sample_sums(bacteria_v34))[which(sample_sums(bacteria_v34)<11504)]
length(sample_sums(bacteria_v34))[which(sample_sums(bacteria_v34)<11504)]

# Rarefiying
bacteria_v34_rare <- rarefy_even_depth(bacteria_v34,11504,rngseed = 777)
#`set.seed(777)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(777); .Random.seed` for the full vector
#...
#4 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
#  
#  negctrlTB211TB217PTB219P
#...
#3391OTUs were removed because they are no longer 
#present in any sample after random subsampling

# Removing non-turtle samples
bacteria_v34_rare_turtle <- subset_samples(bacteria_v34_rare, SampleType=="carapace")
bacteria_v34_rare_turtle

# Alpha diversity
alpha_diversity_v34_bacteria <- estimate_richness(bacteria_v34_rare,measures=c("Observed","Shannon","InvSimpson"))
head(alpha_diversity_v34_bacteria)

alpha_v34_observed_age_bacteria <- plot_richness(bacteria_v34_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v34_observed_age_bacteria

alpha_v34_shannon_age_bacteria <- plot_richness(bacteria_v34_rare_turtle, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index",color = "Age", shape = "Age") +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v34_shannon_age_bacteria


alpha_v34_observed_age_bacteria + alpha_v34_shannon_age_bacteria
ggsave(filename = "figures/alpha_v34_obsv_shannon_age_bacteria.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_obsv_v34_shannon_age_bacteria.jpg", 
       width = 6.75, height = 4, dpi = 300)

# preparation of data for Pearson correlation test
bacteria_v34_rare_turtle_meta<-data.frame(sample_data(bacteria_v34_rare_turtle))
alpha_diversity_v34_bacteria$Sample<-rownames(bacteria_v34_rare_turtle_meta) #making column "Sample"
bacteria_v34_rare_turtle_meta$Sample<-rownames(bacteria_v34_rare_turtle_meta) #making column "Sample"
alpha_diversity_v34_bacteria <-left_join(alpha_diversity_v34_bacteria,bacteria_v34_rare_turtle_meta, by="Sample")

# Pearson correlation test
with(alpha_diversity_v34_bacteria,cor.test(CCL,Shannon))
 

with(alpha_diversity_v34,cor.test(CCL,Observed))


# NMDS
# NMDS - premalo podataka za dobar prikaz !!!
ord_nmds_bray_v34_bacteria <- ordinate(bacteria_v34_rare_turtle, method="NMDS",k=3, distance="bray",trymax=50)
nmds_bray_v34_bacteria <- plot_ordination(bacteria_v34_rare_turtle, ord_nmds_bray_v34_bacteria, color="Age", shape ="Age", title="Bray-Curtis NMDS V34") +
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  stat_ellipse(geom="polygon",type="norm",aes(colour=Age, fill=Age), alpha=0.2) + 
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw()
nmds_bray_v34_bacteria
# In metaMDS(veganifyOTU(physeq), distance, ...) :
# stress is (nearly) zero: you may have insufficient data
# Too few points to calculate an ellipse

ggsave(filename = "figures/nmds_v34_bray_age_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/nmds_v34_bray_age_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

# Combined plot of alpha and beta diversity
(alpha_v34_observed_age_bacteria + alpha_v34_shannon_age_bacteria) / nmds_bray_v34_bacteria
ggsave(filename = "figures/alpha_nmds_v34_bacteria.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/alpha_nmds_v34_bacteria.jpg", width = 6.75, height = 8, dpi = 300)

# Heatmaps

bacteria_v34_turtle_top20 <- prune_taxa(names(sort(taxa_sums(bacteria_v34_turtle),TRUE)[1:20]), bacteria_v34_turtle)
plot_heatmap(bacteria_v34_turtle_top20,"NMDS",distance = "bray")

ggsave(filename = "figures/heatmap_v34_asv_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v34_asv_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

plot_heatmap(bacteria_v34_turtle_top20,"NMDS",
             distance = "bray",
             sample.label="Age",
             sample.order = "CCL",
             taxa.label = "Genus")

ggsave(filename = "figures/heatmap_v34_genera_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v34_genera_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

#Pheatmap package - Generate Metadata for Plotting Colours
age_data<-data.frame(Age=sample_data(bacteria_v34_turtle_top20)$Age)
rownames(age_data)<-rownames(sample_data(bacteria_v34_turtle_top20))

#Plot Heatmap 
pheatmap(otu_table(bacteria_v34_turtle_top20),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale="row",
         annotation_col = age_data)

# Taxa bar plots
bacteria_phylum_v34 <- bacteria_v34_turtle %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")

plot_composition(bacteria_phylum_v34)

ggsave(filename = "figures/taxa_bar_plot_v34_phyla_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v34_phyla_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

#### Run Again But With Top 'N' Phyla 

#What Are the Names of the most abundant phyla?  
bacteria_v34_phyla_collapse<- bacteria_v34_turtle %>% aggregate_taxa(level="Phylum")
bacteria_v34_top10phyla = names(sort(taxa_sums(bacteria_v34_phyla_collapse), TRUE)[1:10])

#Subset the phyloseq object to those phyla   
bacteria_v34_top10phyla_filter<-subset_taxa(bacteria_v34_turtle,Phylum %in% bacteria_v34_top10phyla)

#Remake Our Graph  but with grouping by CCL

bacteria_v34_top10phyla_plot <- bacteria_v34_top10phyla_filter %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "CCL")
bacteria_v34_top10phyla_plot + scale_fill_brewer(palette="Paired")

ggsave(filename = "figures/taxa_bar_plot_v34_phyla_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v34_phyla_bacteria.jpg", width = 6.75, height = 4, dpi = 300)



### BETA DIVERISTY PLOT AND TESTS BASED ON CLR-TRANSFORMED DATA OF UNRAREFIED LIBRARY SIZES (XAH)

#### Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu_v34_bacteria <- function(bacteria_v34) {
  OTU_v34_bacteria <- otu_table(bacteria_v34)
  if (taxa_are_rows(OTU_v34_bacteria)) {
    OTU_v34_bacteria <- t(OTU_v34_bacteria)
  }
  return(as(OTU_v34_bacteria, "matrix"))
}

#CLR Transform  on unrarefied object
bacteria_v34_clr <- microbiome::transform(bacteria_v34, "clr")
head(otu_table(bacteria_v34_clr))

#Extract Matrix and Sample Data  
bacteria_v34_clr_vegan<-vegan_otu_v34_bacteria(bacteria_v34_clr)
bacteria_v34_clr_sample<-as(sample_data(bacteria_v34_clr),"data.frame")

######## Principal Components Analysis
bacteria_v34_pc<-prcomp(bacteria_v34_clr_vegan)

#Scree plot of relative importance of explained by each axis
plot(bacteria_v34_pc)

#Variance Explained by FIrst 2 axes  
summary(bacteria_v34_pc)$importance[,1:4] # 

#### Extract Scores for Plotting        
bacteria_v34_pc_scores<-scores(bacteria_v34_pc)
bacteria_v34_pc_scores_sub<-bacteria_v34_pc_scores[,1:2]
#Add Sample Data  
bacteria_v34_pc_scores_sub<-cbind(bacteria_v34_pc_scores_sub,bacteria_v34_clr_sample)

#Plot    
pca_v34_bacteria<-ggplot(bacteria_v34_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level=0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks")) +
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks")) +
  theme_bw() + 
  labs(x="PC1 (10.4%)",y="PC2 (7.8%)") 
pca_v34_bacteria

ggsave(filename = "figures/pca_clr_v34_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v34_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

### CLR beta diversity - only turtle samples

vegan_otu_v34_turtle_bacteria <- function(bacteria_v34_turtle) {
  OTU_v34_turtle_bacteria <- otu_table(bacteria_v34_turtle)
  if (taxa_are_rows(OTU_v34_turtle_bacteria)) {
    OTU_v34_turtle_bacteria <- t(OTU_v34_turtle_bacteria)
  }
  return(as(OTU_v34_turtle_bacteria, "matrix"))
}

#CLR Transform  on unrarefied object
bacteria_v34_turtle_clr <- microbiome::transform(bacteria_v34_turtle, "clr")
head(otu_table(bacteria_v34_turtle_clr))

#Extract Matrix and Sample Data  
bacteria_v34_turtle_clr_vegan<-vegan_otu_v34_bacteria(bacteria_v34_turtle_clr)
bacteria_v34_turtle_clr_sample<-as(sample_data(bacteria_v34_turtle_clr),"data.frame")

######## Principal Components Analysis
bacteria_v34_turtle_pc<-prcomp(bacteria_v34_turtle_clr_vegan)

#Scree plot of relative importance of explained by each axis
plot(bacteria_v34_turtle_pc)

#Variance Explained by FIrst 2 axes  
summary(bacteria_v34_turtle_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
bacteria_v34_turtle_pc_scores<-scores(bacteria_v34_turtle_pc)
bacteria_v34_turtle_pc_scores_sub<-bacteria_v34_turtle_pc_scores[,1:2]
#Add Sample Data  
bacteria_v34_turtle_pc_scores_sub<-cbind(bacteria_v34_turtle_pc_scores_sub,bacteria_v34_turtle_clr_sample)

#Plot    
pca_v34_turtle_bacteria<-ggplot(bacteria_v34_turtle_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level = 0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  geom_text(aes(label = AnimalName)) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x="PC1 (11.2%)",y="PC2 (8.0%)") 
pca_v34_turtle_bacteria

ggsave(filename = "figures/pca_clr_v34_turtles_bacteria.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v34_turtles_bacteria.jpg", width = 6.75, height = 4, dpi = 300)

