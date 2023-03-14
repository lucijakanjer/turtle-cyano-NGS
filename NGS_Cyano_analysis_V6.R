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

#Some Boxplots of Coverage by Population using ggplot2
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

alpha_invsimpson_CCL <- plot_richness(cyano50_rare, x = "CCL", measures = c("InvSimpson")) + 
  geom_point(aes(color = CCL), size = 5) + 
  scale_color_continuous(type = "viridis") +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Inversed Simpson index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_invsimpson_CCL

alpha_shannon_CCL <- plot_richness(cyano50_rare, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(color = CCL), size = 3) + 
  scale_color_continuous(type = "viridis", name = "CCL (cm)") +
  theme_bw() + 
  labs(x = "Turtle curved carapace length - CCL (cm)",y = "Shannon diversity index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_shannon_CCL



alpha_observed_CCL <- plot_richness(cyano50_rare, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(color = CCL), size = 3) + 
  scale_color_continuous(type = "viridis") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "Turtle curved carapace length - CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_observed_CCL

install.packages("patchwork")
library(patchwork)

alpha_observed_CCL + alpha_shannon_CCL
ggsave(filename = "alpha_obsv_shannon_CCL.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "alpha_obsv_shannon_CCL_horizontal.jpg", width = 6.75, height = 4, dpi = 300)

alpha_observed_age <- plot_richness(cyano50_rare, x = "CCL", measures = c("Observed")) + 
  geom_point(aes(colour = Age), size = 3) + 
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = "CCL (cm)",y = "Observed ASV richness") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_observed_age

alpha_shannon_age <- plot_richness(cyano50_rare, x = "CCL", measures = c("Shannon")) + 
  geom_point(aes(colour = Age), size = 3) + 
  scale_color_brewer(palette = "Set2", name ="Age") +
  theme_bw() + 
  labs(x = "CCL (cm)",y = "Shannon diversity index") + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))
alpha_shannon_age

alpha_observed_age + alpha_shannon_age
ggsave(filename = "alpha_obsv_shannon_age.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "alpha_obsv_shannon_age.jpg", width = 6.75, height = 4, dpi = 300)


#alpha_diversity$Sample<-rownames(alpha_diversity)
#alpha_diversity<-left_join(alpha_diversity,cyano_meta,"Sample")

#Plot the correlation between heterozygosity and richness (shannon diversity)
#ggplot(alpha_diversity,aes(x=CCL,y=Shannon)) + 
#  geom_point(size=5,fill="white") + 
#  theme_bw() + 
#  labs(x="CCL(cm)", y="Alpha Diversity (Shannon)") + 
#  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))

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


#################### ORDINATION  USING VEGAN   


#### Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu <- function(cyano) {
  OTU <- otu_table(cyano)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

#Load Libraries
#library(vegan)
#library(RColorBrewer)

#Convert OTU table to abundance matrix
cyano.v<-vegan_otu(cyano50_rare)
cyano_turtles.v<-vegan_otu(cyano50_rare_turtles)

#Convert Sample Data to     
cyano_turtles.s<-as(sample_data(cyano50_rare_turtles),"data.frame")

###### NMDS Ordination
cyano_turtles.nmds<-metaMDS(cyano_turtles.v,k=3,distance="bray",trymax = 50) 

stressplot(cyano.nmds)

##################### NMDS Plot using VEGAN


#Select Some Colours from RColorBrewer  
ordination_cols<-brewer.pal(5,"Set2")
cols<-data.frame(pop=unique(cyano_turtles.s$Age),
                 col=ordination_cols[1:length(unique(cyano_turtles.s$Age))])

#Expand so that each sample in the database has a colour value  
cols.expand<-as.character(cols[,2][match(cyano_turtles.s$Age,cols[,1])])
cols<-cols[order(cols[,1]),]

#Plot Axis Bounds
plot(cyano_turtles.nmds,display="sites",type="n",las=1,ylab="",xaxt="n",yaxt="n",xlab="")

#Add On A Spider linking each point to its population centroid  
ordispider(cyano_turtles.nmds,cyano_turtles.s$Age,col="black",label=F)


#Add an ellipse 
#'ehull' is ellipsoid hull that encloses all points within a group
#'#?ordiellipse for more options under option 'kind'
#ordiellipse(cyano_turtles.nmds,cyano_turtles.s$Age,col="lightgray",label=F,kind="ehull")

#Add the data points  
with(cyano_turtles.s,points(cyano_turtles.nmds,display="sites",pch=21,bg=cols.expand,cex=2))

#Add Axis Labels  
#mtext(c('NMDS 1', 'NMDS 2'), side=c(1,2), adj=1, cex=1.5,font=2, line=c(1, 0))

with(cyano.s,legend("bottomright",legend=cols[,1],pch=21,pt.bg=as.character(cols[,2]),bty="n",pt.cex=2))

############### BETA DIVERISTY PLOT AND TESTS BASED ON CLR-TRANSFORMED DATA OF UNRAREFIED LIBRARY SIZES (XAH)


#### Vegan Function to Convert Phyloseq into Something Vegan can call directly
vegan_otu <- function(cyano) {
  OTU <- otu_table(cyano)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

########### Transform Sample Counts
library(microbiome)

#CLR Transform  on unrarefied object
cyano_clr <- microbiome::transform(cyano50, "clr")
head(otu_table(cyano_clr))

#Extract Matrix and Sample Data  
cyano_clr_v<-vegan_otu(cyano_clr)
cyano_clr_s<-as(sample_data(cyano_clr),"data.frame")

######## Principal Components Analysis
cyano_pc<-prcomp(cyano_clr_v)

#Cool Biplot Showing How Diff ASVs affect the primary axes of the ordinatiton
biplot(cyano_pc)

#Scree plot of relative importance of explained by each axis
plot(cyano_pc)

#Variance Explained by FIrst 2 axes  
summary(cyano_pc)$importance[,1:2] # 

#### Extract Scores for Plotting        
library(vegan)
cyano_pc_scores<-scores(cyano_pc)
cyano_pc_scores_sub<-cyano_pc_scores[,1:2]
#Add Sample Data  
cyano_pc_scores_sub<-cbind(cyano_pc_scores_sub,cyano_clr_s)

######### Plot
library(ggplot2)
library(cowplot)

#Housekeeping + Label Setup    
#cyano_pc_scores_sub$Species<-as.factor(cyano_pc_scores_sub$Full_names)


#Plot    
pca1<-ggplot(cyano_pc_scores_sub,aes(x=PC1,y=PC2)) + 
  stat_ellipse(type="t",aes(color=Age),level = 0.95,alpha=0.5) + 
  geom_point(aes(colour=Age),size=5) 

pca2<-pca1 + 
  theme_bw() + 
  labs(fill="Transect Name",x="PC1 (11.7%)",y="PC2 (7.3%)") 

pca3<-pca2 +  
  theme(legend.position="top",
        axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text = element_text(size=12),
        legend.title = element_text(size=18))
pca3



