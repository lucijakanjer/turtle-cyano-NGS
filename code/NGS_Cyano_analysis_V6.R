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

#Stats Models
library(lme4)
library(MuMIn)
library(car)

### Build Phyloseq Object 
cyano_v6 <- qza_to_phyloseq(
  "data/table-V6-cyanobacteria.qza",
  "data/rooted-tree-V6-cyanobacteria.qza",
  "data/taxonomy-V6-cyanobacteria.qza",
  "data/metadata-V6-v1.tsv") #sample metadata must be in TSV format
cyano_v6

# Building phyloseq object without the phylogenetic tree
cyano_v6_no_tree <- qza_to_phyloseq(
  features = "data/table-V6-cyanobacteria.qza",
  taxonomy = "data/taxonomy-V6-cyanobacteria.qza", 
  metadata = "data/metadata-V6-v1.tsv")
cyano_v6_no_tree

# top 20 records from taxonomy table
head(tax_table(cyano_v6),20)

# check how many ASVs are found in different taxa on specific level
table(tax_table(cyano_v6)[,4])

# We can ask what proportion of assignments we have made at each level
apply(tax_table(cyano_v6)[,2:7],2,function(x){1-mean(is.na(x))})


# Filtering based on number of reads

cyano_v6_50 <- prune_taxa(taxa_sums(cyano_v6)>50,cyano_v6)
ntaxa(cyano_v6_50); ntaxa(cyano_v6)
ntaxa(cyano_v6)- ntaxa(cyano_v6_50)

# Quality assesment

mean(sample_sums(cyano_v6_50)); range(sample_sums(cyano_v6_50))

#Make a data frame of read depths
cyano_v6_reads <- data.frame(reads=sample_sums(cyano_v6_50))

#Add on the sample ID
cyano_v6_reads$Sample<-rownames(cyano_v6_reads)

#Extract the Metadata from the phyloseq object 
cyano_v6_meta<-data.frame(sample_data(cyano_v6))
cyano_v6_meta$Sample<-rownames(cyano_v6_meta)

#Join on the Metadata
cyano_v6_reads<-left_join(cyano_v6_reads,cyano_v6_meta,"Sample")

#Some Boxplots of Coverage by Population using ggplot2
#library(ggplot2)
ggplot(cyano_v6_reads,aes(x=Age,y=reads)) + 
  geom_boxplot(aes(fill=Age)) + 
  geom_jitter()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# rarefaction curve
rarecurve(t(otu_table(cyano_v6)), step=50, cex=0.5) 
min(sample_sums(cyano_v6))
names(sample_sums(cyano_v6))[which(sample_sums(cyano_v6)<10000)]
length(sample_sums(cyano_v6))[which(sample_sums(cyano_v6)<10000)]

# Rarefiying
cyano_v6_rare <- rarefy_even_depth(cyano_v6,10000,rngseed = 777)
# set.seed(777)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(777); .Random.seed` for the full vector
# ...
# 10 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#  
#   V6-negctrlV6-TB163V6-TB167V6-TB207V6-TB213
# ...
# 66TUs were removed because they are no longer 
# present in any sample after random subsampling

cyano_v6_turtle <- prune_samples(sample_data(cyano_v6)$SampleType =="carapace",cyano_v6)
cyano_v6_turtle

cyano_v6_rare_turtle <- prune_samples(sample_data(cyano_v6_rare)$SampleType =="carapace",cyano_v6_rare)
cyano_v6_rare_turtle

###################   SHARED AND UNIQUE SEQUENCES BY AGE
#library(gplots) #required for Venn diagram. Online Venn diagram tools like Venny are also available. 

#Empty lists - we have two groups so we need two empty lists
venn.list<-rep(list(NA),3)
poplist<-rep(list(NA),3)
pops<-as.character(unique(sample_data(cyano_v6_turtle)$Age))

#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]



#Name the Groups in the Vennlist
names(venn.list)<-pops

#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:3){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(cyano_v6_rare)$Age %in% poplist[[k]],cyano_v6_turtle)
  
  #Calculate which ASVs are present (non-zero abundance)  
  phy.sub.keep<-apply(otu_table(phy.sub),1,function(x){sum(x)>0})
  
  #Prune down to Inly those ASVs  
  phy.sub.sub<-prune_taxa(phy.sub.keep,phy.sub)
  
  #Extract the names of these ASVs and store them in the appropriate list slot  
  venn.list[[k]]<-rownames(otu_table(phy.sub.sub))
} #loop end


#Draw Venn Diagram, automatically calculate overlap in the list and unique variants
venn(venn.list)

#Alpha diversity

alpha_diversity_v6 <- estimate_richness(cyano_v6_rare_turtle,measures=c("Observed","Shannon","InvSimpson"))
head(alpha_diversity_v6)

plot_richness(cyano_v6_rare,x="CCL",measures=c("Observed"))

alpha_v6_invsimpson_CCL <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("InvSimpson")) + 
  geom_point(aes(color = CCL), size = 5) + 
  scale_color_continuous(type = "viridis") +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Inversed Simpson index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_invsimpson_CCL

alpha_v6_shannon_CCL <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(color = CCL), size = 3) + 
  scale_color_continuous(type = "viridis", name = "CCL (cm)") +
  theme_bw() + 
  labs(x = "Turtle curved carapace length - CCL (cm)",y = "Shannon diversity index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_shannon_CCL



alpha_v6_observed_CCL <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(color = CCL), size = 3) + 
  scale_color_continuous(type = "viridis") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "Turtle curved carapace length - CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_observed_CCL


alpha_v6_observed_CCL + alpha_v6_shannon_CCL
ggsave(filename = "figures/alpha_v6_obsv_shannon_CCL.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_v6_obsv_shannon_CCL.jpg", width = 6.75, height = 4, dpi = 300)

alpha_v6_observed_age <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = Age), size = 3) + 
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_observed_age

alpha_v6_shannon_age <- plot_richness(cyano_v6_rare_turtle, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(colour = Age, shape = Age), size = 3) + 
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index",color = "Age", shape = "Age") +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_v6_shannon_age

alpha_v6_observed_age + alpha_v6_shannon_age
ggsave(filename = "figures/alpha_v6_obsv_shannon_age.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/alpha_v6_obsv_shannon_age.jpg", width = 6.75, height = 4, dpi = 300)



# preparation of data for Pearson correlation test
cyano_v6_rare_turtle_meta<-data.frame(sample_data(cyano_v6_rare_turtle))
alpha_diversity_v6$Sample<-rownames(cyano_v6_rare_turtle_meta) #making column "Sample"
cyano_v6_rare_turtle_meta$Sample<-rownames(cyano_v6_rare_turtle_meta) #making column "Sample"
alpha_diversity_v6 <-left_join(alpha_diversity_v6,cyano_v6_rare_turtle_meta, by="Sample")

# Pearson correlation test
with(alpha_diversity_v6,cor.test(CCL,Shannon))
# data:  CCL and Shannon
# t = 4.0231, df = 17, p-value = 0.0008819
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3575996 0.8750205
# sample estimates:
#  cor 
# 0.6983752 

with(alpha_diversity_v6,cor.test(CCL,Observed))
# data:  CCL and Observed
# t = 2.0791, df = 17, p-value = 0.05306
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.004984767  0.750892055
# sample estimates:
#  cor 
# 0.4502439

# NMDS
ord_nmds_bray_v6 <- ordinate(cyano_v6_rare_turtle, method="NMDS",k=3, distance="bray",trymax=50)
nmds_bray_v6 <- plot_ordination(cyano_v6_rare_turtle, ord_nmds_bray_v6, color="Age", shape ="Age", title="Bray-Curtis NMDS V6") +
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  stat_ellipse(geom="polygon",type="norm",aes(colour=Age, fill=Age), alpha=0.2) + 
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw()
nmds_bray_v6

ggsave(filename = "figures/nmds_v6_bray_age.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/nmds_v6_bray_age.jpg", width = 6.75, height = 4, dpi = 300)

# Combined plot of alpha and beta diversity
(alpha_v6_observed_age + alpha_v6_shannon_age) / nmds_bray_v6
ggsave(filename = "figures/alpha_nmds_v6.pdf", width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/alpha_nmds_v6.jpg", width = 6.75, height = 8, dpi = 300)

# Heatmaps

cyano_v6_turtle_top20 <- prune_taxa(names(sort(taxa_sums(cyano_v6_turtle),TRUE)[1:20]), cyano_v6_turtle)
plot_heatmap(cyano_v6_turtle_top20,"NMDS",distance = "bray")

ggsave(filename = "figures/heatmap_v6_asv.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v6_asv.jpg", width = 6.75, height = 4, dpi = 300)

plot_heatmap(cyano_v6_turtle_top20,"NMDS",
             distance = "bray",
             sample.label="Age",
             sample.order = "CCL",
             taxa.label = "Genus")

ggsave(filename = "figures/heatmap_v6_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/heatmap_v6_genera.jpg", width = 6.75, height = 4, dpi = 300)

#Pheatmap package - Generate Metadata for Plotting Colours
age_data<-data.frame(Age=sample_data(cyano_v6_turtle_top20)$Age)
rownames(age_data)<-rownames(sample_data(cyano_v6_turtle_top20))
#Plot Heatmap 
pheatmap(otu_table(cyano_v6_turtle_top20),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale="row",
         annotation_col = age_data)

# Taxa bar plots
cyano_order_v6 <- cyano_v6_turtle %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional")

plot_composition(cyano_order_v6)

cyano_order_v6 %>%
  plot_composition(average_by = "Age")+ scale_fill_brewer(palette="Set3")

ggsave(filename = "figures/taxa_bar_plot_v6_order.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_order.jpg", width = 6.75, height = 4, dpi = 300)

order_v6 <- cyano_order_v6 %>%
  plot_composition(average_by = "SampleType")+ scale_fill_brewer(palette="Set3")
order_v6



#### Run Again But With Top 'N' Phyla 

#What Are the Names of the most abundant phyla?  
cyano_v6_genus_collapse<- cyano_v6_turtle %>% aggregate_taxa(level="Genus")
cyano_v6_top10genera = names(sort(taxa_sums(cyano_v6_genus_collapse), TRUE)[1:10])

#Subset the phyloseq object to those phyla   
cyano_v6_top10genera_filter<-subset_taxa(cyano_v6_turtle,Genus %in% cyano_v6_top10genera)

#Remake Our Graph  but with grouping by CCL

cyano_v6_top10genera_plot <- cyano_v6_top10genera_filter %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "CCL")
cyano_v6_top10genera_plot + scale_fill_brewer(palette="Set3")

ggsave(filename = "figures/taxa_bar_plot_v6_genera.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_genera.jpg", width = 6.75, height = 4, dpi = 300)


### BETA DIVERISTY PLOT AND TESTS BASED ON CLR-TRANSFORMED DATA OF UNRAREFIED LIBRARY SIZES (XAH)

#### Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu_v6 <- function(cyano_v6) {
  OTU_v6 <- otu_table(cyano_v6)
  if (taxa_are_rows(OTU_v6)) {
    OTU_v6 <- t(OTU_v6)
  }
  return(as(OTU_v6, "matrix"))
}

#CLR Transform  on unrarefied object
cyano_v6_clr <- microbiome::transform(cyano_v6, "clr")
head(otu_table(cyano_v6_clr))

#Extract Matrix and Sample Data  
cyano_v6_clr_vegan<-vegan_otu_v6(cyano_v6_clr)
cyano_v6_clr_sample<-as(sample_data(cyano_v6_clr),"data.frame")

######## Principal Components Analysis
cyano_v6_pc<-prcomp(cyano_v6_clr_vegan)

#Scree plot of relative importance of explained by each axis
plot(cyano_v6_pc)

#Variance Explained by FIrst 2 axes  
summary(cyano_v6_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
library(vegan)
cyano_v6_pc_scores<-scores(cyano_v6_pc)
cyano_v6_pc_scores_sub<-cyano_v6_pc_scores[,1:2]
#Add Sample Data  
cyano_v6_pc_scores_sub<-cbind(cyano_v6_pc_scores_sub,cyano_v6_clr_sample)

#Plot    
pca_v6<-ggplot(cyano_v6_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level=0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult", "rocks", "pool", "negctrl")) +
  theme_bw() + 
  labs(x="PC1 (14.7%)",y="PC2 (12.8%)")
pca_v6

ggsave(filename = "figures/pca_clr_v6.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v6.jpg", width = 6.75, height = 4, dpi = 300)

### CLR beta diversity - only turtle samples

vegan_otu_v6_turtle <- function(cyano_v6_turtle) {
  OTU_v6_turtle <- otu_table(cyano_v6_turtle)
  if (taxa_are_rows(OTU_v6_turtle)) {
    OTU_v6_turtle <- t(OTU_v6_turtle)
  }
  return(as(OTU_v6_turtle, "matrix"))
}

#CLR Transform  on unrarefied object
cyano_v6_turtle_clr <- microbiome::transform(cyano_v6_turtle, "clr")
head(otu_table(cyano_v6_turtle_clr))

#Extract Matrix and Sample Data  
cyano_v6_turtle_clr_vegan<-vegan_otu_v6(cyano_v6_turtle_clr)
cyano_v6_turtle_clr_sample<-as(sample_data(cyano_v6_turtle_clr),"data.frame")

######## Principal Components Analysis
cyano_v6_turtle_pc<-prcomp(cyano_v6_turtle_clr_vegan)

#Scree plot of relative importance of explained by each axis
plot(cyano_v6_turtle_pc)

#Variance Explained by FIrst 2 axes  
summary(cyano_v6_turtle_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
#library(vegan)
cyano_v6_turtle_pc_scores<-scores(cyano_v6_turtle_pc)
cyano_v6_turtle_pc_scores_sub<-cyano_v6_turtle_pc_scores[,1:2]
#Add Sample Data  
cyano_v6_turtle_pc_scores_sub<-cbind(cyano_v6_turtle_pc_scores_sub,cyano_v6_turtle_clr_sample)

#Plot    
pca_v6_turtle<-ggplot(cyano_v6_turtle_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(geom="polygon",type="t",aes(color=Age, fill=Age),level = 0.95,alpha=0.2) + 
  geom_point(aes(colour = Age, shape = Age), size = 2) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_fill_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  theme_bw() + 
  labs(x="PC1 (17.2%)",y="PC2 (10.6%)") 
pca_v6_turtle

ggsave(filename = "figures/pca_clr_v6_turtles.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/pca_clr_v6_turtles.jpg", width = 6.75, height = 4, dpi = 300)
