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

### Build Phyloseq Object 
cyano <- qza_to_phyloseq(
  "table-V6-cyanobacteria.qza",
  "rooted-tree-V6-cyanobacteria.qza",
  "taxonomy-V6-cyanobacteria.qza",
  "metadata-V6-v1.tsv") #sample metadata must be in TSV format
cyano

# Building phyloseq object without the phylogenetic tree
cyano_no_tree <- qza_to_phyloseq(
  features = "table-V6-cyanobacteria.qza",
  taxonomy = "taxonomy-V6-cyanobacteria.qza", 
  metadata = "metadata-V6-v1.tsv")
cyano_no_tree

# top 20 records from taxonomy table
head(tax_table(cyano),20)

# check how many ASVs are found in different taxa on specific level
table(tax_table(cyano)[,4])

# We can ask what proportion of assignments we have made at each level
apply(tax_table(cyano)[,2:7],2,function(x){1-mean(is.na(x))})


# Filtering based on number of reads

cyano50 <-prune_taxa(taxa_sums(cyano)>50,cyano)
ntaxa(cyano50); ntaxa(cyano)
ntaxa(cyano)- ntaxa(cyano50)

# Quality assesment

mean(sample_sums(cyano50)); range(sample_sums(cyano50))

#Make a data frame of read depths
cyano_reads <- data.frame(reads=sample_sums(cyano50))

#Add on the sample ID
cyano_reads$Sample<-rownames(cyano_reads)

#Extract the Metadata from the phyloseq object 
cyano_meta<-data.frame(sample_data(cyano))
cyano_meta$Sample<-rownames(cyano_meta)

#Join on the Metadata
cyano_reads<-left_join(cyano_reads,cyano_meta,"Sample")

#Some Boxplots of COverage by Population using ggplot2
#library(ggplot2)
ggplot(cyano_reads,aes(x=Age,y=reads)) + 
  geom_boxplot(aes(fill=Age)) + 
  geom_jitter()
  #plotopts + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# rarefaction curve
rarecurve(t(otu_table(cyano50)), step=50, cex=0.5) 
# ne radi, javlja gresku: Error in as(x, "matrix")[i, j, drop = FALSE] : 
# (subscript) logical subscript too long

min(sample_sums(cyano50))
names(sample_sums(cyano50))[which(sample_sums(cyano50)<10000)]
length(sample_sums(cyano50))[which(sample_sums(cyano50)<10000)]

# Rarefiying
cyano50_rare <- rarefy_even_depth(cyano50,10000,rngseed = 777)
# set.seed(777)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(777); .Random.seed` for the full vector
# ...
# 10 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#  
#   V6-negctrlV6-TB163V6-TB167V6-TB207V6-TB213
# ...
# 8OTUs were removed because they are no longer 
# present in any sample after random subsampling

cyano50_rare_turtles <- prune_samples(sample_data(cyano50_rare)$SampleType =="carapace",cyano50_rare)
cyano50_rare_turtles

###################   SHARED AND UNIQUE SEQUENCES BY SITE
#library(gplots) #required for Venn diagram. Online Venn diagram tools like Venny are also available. 

#Empty lists - we have two groups so we need two empty lists
venn.list<-rep(list(NA),4)
poplist<-rep(list(NA),4)
pops<-as.character(unique(sample_data(cyano50_rare)$Age))

#Populate List with the Names of Groups we Want to Compare
poplist[[1]]<-pops[1]
poplist[[2]]<-pops[2]
poplist[[3]]<-pops[3]
poplist[[4]]<-pops[4]


#Name the Groups in the Vennlist
names(venn.list)<-pops

#Loop over each population, subset the phyloseq object to just that population and work out which SVs are in that population, store those names in the lisr
for(k in 1:4){
  
  #Subset Phyloseq Object to One Pop at a time  
  phy.sub<-prune_samples(sample_data(cyano50_rare)$Age %in% poplist[[k]],cyano50_rare)
  
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

alpha_diversity <- estimate_richness(cyano50_rare,measures=c("Observed","Shannon","InvSimpson"))
head(alpha_diversity)

plot_richness(cyano50_rare,x="CCL",measures=c("Observed"))

plot_richness(cyano50_rare, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(color = CCL), size = 5) + 
  scale_color_continuous(type = "viridis") +
  theme_bw() + 
  labs(x = "CCL(cm)",y = "Alpha Diversity (Shannon)") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))

#alpha_diversity$Sample<-rownames(alpha_diversity)
#alpha_diversity<-left_join(alpha_diversity,cyano_meta,"Sample")

#Plot the correlation between heterozygosity and richness (shannon diversity)
ggplot(alpha_diversity,aes(x=CCL,y=Shannon)) + 
  geom_point(size=5,fill="white") + 
  theme_bw() + 
  labs(x="CCL(cm)", y="Alpha Diversity (Shannon)") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))

# Pearson test for correlation - but do not work because I can't join tables because 
# of different names in "alpha_diversity" and "cyano_meta"
with(alpha_diversity,cor.test(CCL,Shannon))

# NMDS
ord_nmds_bray <- ordinate(cyano50_rare_turtles, method="NMDS",k=3, distance="bray",trymax=50)
plot_ordination(cyano50_rare_turtles, ord_nmds_bray, color="Age", title="Bray NMDS") + 
  stat_ellipse(geom="polygon",aes(fill=Age),type="norm",alpha=0.3) + 
  theme_bw()

# Heatmaps

cyano50_rare_top20 <- prune_taxa(names(sort(taxa_sums(cyano50_rare_turtles),TRUE)[1:20]), cyano50_rare_turtles)
plot_heatmap(cyano50_rare_top20,"NMDS",distance = "bray")

plot_heatmap(cyano50_rare_top20,"NMDS",
             distance = "bray",
             sample.label="CCL",
             sample.order = "CCL",
             taxa.label = "Genus")

#Pheatmap package - Generate Metadata for Plotting Colours
age_data<-data.frame(Age=sample_data(cyano50_rare_top20)$Age)
rownames(age_data)<-rownames(sample_data(cyano50_rare_top20))


#Plot Heatmap 
pheatmap(otu_table(cyano50_rare_top20),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale="row",
         annotation_col = age_data)

# Taxa bar plots
cyano_order <- cyano50_rare_turtles %>%
  aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional")

plot_composition(cyano_order)

cyano_order %>%
  plot_composition(average_by = "Age")+ scale_fill_brewer(palette="Set3")

#### Run Again But With Top 'N' Phyla 

#What Are the Names of the most abundant phyla?  
cyano_genus_collapse<- cyano50_rare_turtles %>% aggregate_taxa(level="Genus")
cyano_top10genera = names(sort(taxa_sums(cyano_genus_collapse), TRUE)[1:10])

#Subset the phyloseq object to those phyla   
cyano_top10genera_filter<-subset_taxa(cyano50_rare_turtles,Genus %in% cyano_top10genera)

#Remake Our Graph  but with grouping by VEGETATION

cyano_top10genera_plot <- cyano_top10genera_filter %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "CCL")
cyano_top10genera_plot + scale_fill_brewer(palette="Set3")

