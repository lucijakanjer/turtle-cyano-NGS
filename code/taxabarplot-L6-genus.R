# Loading required packages
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(patchwork)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77",  
                             "#332288", "#117733","#AA4499", 
                             "#44AA99", "#999933", "#882255", 
                             "#E69F00", "#661100", "#6699CC", 
                             "#999999")

# Importing dataset 
genus_raw <- read.csv(file = "data/taxabarplots-cyano-V6-level-6-genus.csv", 
                       header = TRUE)

# Creating working table version
genus <- genus_raw

# Visualizing column names
print(colnames(genus))

# Function to keep only the last part of the column name
shorten_name <- function(name) {
  sub(".*__", "", name)
}

# Renaming columns
# creating function to change the name
new_names <- sapply(names(genus), function(name) {
  if (grepl("k__Bacteria.p__Cyanobacteriota.c__Cyanophyceae", name)) {
    shorten_name(name)
  } else {
    name
  }
})

# Assigning new column names to the data frame
names(genus) <- new_names

# Visualizing new column names
print(colnames(genus))

# Renaming unclassified families
# Indices of columns to rename
indices <- c(6, 9, 11, 13, 22, 25, 26, 36, 39, 41, 45, 52, 64, 65)

# Vector with new names for specific columns
unclassified_names <- c("unclassified_Cyanophyceae", #6
                        "unclassified_Caenarcaniphilales", #9
                        "unclassified_Cyanobacteriota", #11
                        "unclassified_Sericytochromatia", #13
                        "unclassified_Nodosilineales", #22
                        "unclassified_Chroococcales", #25
                        "unclassified_Nodosilineaceae", #26
                        "unclassified_Pleurocapsaceae", #36
                        "unclassified_Microcoleaceae", #39
                        "unclassified_Prochlorococcaceae", #41
                        "unclassified_Nostocales", #45
                        "unclassified_Oculatellaceae", #52
                        "unclassified_Rivulariaceae", #64
                        "unclassified_Aphanizomenonaceae" #65
                        )

# Rename columns by their indices
names(genus)[indices] <- unclassified_names

# Print the updated column names
print(colnames(genus))

# Converting the table in long format for plotting
cyano_cols <- names(genus)[2:67]
genus_long <- genus %>%
  pivot_longer(
    cols = all_of(cyano_cols),
    names_to = "genus",
    values_to = "count"
  )

# Changing sample names
genus_long <- genus_long %>%
  mutate(index = str_remove(index, "V6-"))

# Changing names
genus_long <- genus_long %>%
  mutate(Age = case_when(
    Age == "rocks" ~ "R",
    Age == "pool" ~ "P",
    TRUE ~ Age  # Keep other values unchanged
  ))

# Removing negative control sample
genus_long <- genus_long %>%
  filter(index != "negctrl")


# Making stacked bar plot - all families
plot_genus <- ggplot(genus_long, aes(x = index, y = count, fill = genus)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Samples", y = "Relative abundance", tag = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  #scale_fill_brewer(palette = "Dark")+
  facet_grid(.~ factor(Age, 
                       levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
plot_genus
ggplotly(plot_genus)

# Saving last generated plot
ggsave(filename = "figures/taxa_bar_plot_genus_all.pdf", 
       width = 15, height = 5, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_genus_all.jpg", 
       width = 15, height = 5, dpi = 300)

### Filtering and plotting only top 10 orders
# Defining metadata columns
metadata <- genus_long %>% select(sample(1:15)) %>% distinct()

# Calculate total counts of each organism across the entire dataset
total_counts <- genus_long %>%
  group_by(genus) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))

# Identify the top 12 organisms
top_12_genus <- total_counts %>%
  slice_head(n = 12) %>%
  pull(genus)

# View the top 10 families
print(top_12_genus)

# Label the top 10 families and aggregate the rest as "other"
genus_labeled <- genus_long %>%
  mutate(genus = if_else(genus %in% top_12_genus, genus, "other")) %>%
  group_by(index, genus) %>%
  summarise(count = sum(count), .groups = 'drop')

# Ensure the family column is a factor with correct levels
genus_labeled <- genus_labeled %>%
  mutate(genus = factor(genus, levels = c(top_12_genus, "other")))

# Merge Back with Metadata
genus_labeled <- genus_labeled %>%
  left_join(metadata, by = "index")

# Making stacked bar plot - top 12 genera
plot_genus_top12 <- ggplot(genus_labeled, 
                            aes(x = index, y = count, fill = genus)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Samples", y = "Relative abundance", tag = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  #scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, 
                       levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
plot_genus_top12

# Saving last generated plot
ggsave(filename = "figures/taxa_bar_plot_family_top12.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_genus_top12.jpg", 
       width = 6.75, height = 4, dpi = 300)

plot_order_top12 / plot_genus_top12

ggsave(filename = "figures/taxa_bar_plot_comb.pdf", 
       width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_comb.jpg", 
       width = 6.75, height = 8, dpi = 300)
