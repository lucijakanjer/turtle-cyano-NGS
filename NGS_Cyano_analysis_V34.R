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
cyano_v34_rare <- rarefy_even_depth(cyano50,1141,rngseed = 777)
# set.seed(777)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(777); .Random.seed` for the full vector
# ...
# 5 samples removed because they contained fewer reads than `sample.size`.
# "Up to first five removed samples are: 
#  
#  V6-negctrl V6-TB167 V6-TB213 V6-TB217P V6-TB219P
# ...
# 41OTUs were removed because they are no longer 
# present in any sample after random subsampling

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
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v34_shannon_age

#library(patchwork)
alpha_v34_observed_age + alpha_v34_shannon_age
ggsave(filename = "figures/alpha_v34_obsv_shannon_age.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_obsv_v34_shannon_age.jpg", 
       width = 6.75, height = 4, dpi = 300)


alpha_diversity_v34$Sample<-rownames(cyano_v34_meta)
alpha_diversity_v34 <-left_join(alpha_diversity_v34,cyano_v34_meta,"CCL")
with(alpha_diversity_v34,cor.test(CCL,Shannon))

#Make a data frame of read depths
cyano_v34_alpha <- data.frame(alpha_diversity_v34)

#Add on the sample ID
cyano_v34_alpha$Sample<-rownames(cyano_v34_alpha)

#Extract the Metadata from the phyloseq object 
cyano_CCL<-data.frame(sample_data$CCL(cyano_v34))
cyano_meta$Sample<-rownames(cyano_v34_meta)

#Join on the Metadata
cyano_reads<-left_join(cyano_reads,cyano_meta,"Sample")

cyano_v34_alpha <-
  left_join(alpha_diversity_v34,
            cyano_v34_meta %>% dplyr::select(CCL), by ="Observed")
