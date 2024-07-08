# Loading required packages
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77",  
                             "#332288", "#117733","#AA4499", 
                             "#44AA99", "#999933", "#882255", 
                             "#E69F00", "#661100", "#6699CC", 
                             "#999999")

# Importing dataset 
order_raw <- read.csv(file = "data/taxabarplots-cyano-V6-level-4-order.csv", 
                      header = TRUE)

# Creating working table version
order <- order_raw

# Visualizing column names
print(colnames(order_raw))

# Changing column names to keep only short names for Orders
for(i in 1:length(colnames(order))) {
  if(grepl("k__Bacteria.p__Cyanobacteriota.", 
           colnames(order)[i])) {
    colnames(order)[i] <- gsub("k__Bacteria.p__Cyanobacteriota.", "",
                                   colnames(order)[i])
  }
}

for(i in 1:length(colnames(order))) {
  if(grepl("c__Cyanophyceae.o__", 
           colnames(order)[i])) {
    colnames(order)[i] <- gsub("c__Cyanophyceae.o__", "",
                               colnames(order)[i])
  }
}

order <- order %>%
  rename(unclassified_Cyanophyceae = "c__Cyanophyceae.__", 
         Caenarcaniphilales = "c__Vampirivibrionia.o__Caenarcaniphilales",
         unclassified_Cyanobacteriota = "__.__",
         unclassified_Sericytochromatia = "c__Sericytochromatia.o__")

# Checking the changed names
print(colnames(order))

# Converting the table in long format for plotting
cyano_cols <- names(order)[2:20]
order <- order %>%
  pivot_longer(
    cols = all_of(cyano_cols),
    names_to = "order",
    values_to = "count"
  )

# Changing sample names
order <- order %>%
  mutate(index = str_remove(index, "V6-"))

# Changing names
order <- order %>%
  mutate(Age = case_when(
    Age == "rocks" ~ "R",
    Age == "pool" ~ "P",
    TRUE ~ Age  # Keep other values unchanged
  ))

# Removing negative control sample
order <- order %>%
  filter(index != "negctrl")


# Making stacked bar plot - all orders
plot_order <- ggplot(order, aes(x = index, y = count, fill = order)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Samples", y = "Relative abundance", tag = "A") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  #scale_fill_brewer(palette = "Dark")+
  facet_grid(.~ factor(Age, 
                       levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
plot_order

# Saving last generated plot
ggsave(filename = "figures/taxa_bar_plot_order_all.pdf", 
       width = 6.75, height = 6, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_order_all.jpg", 
       width = 6.75, height = 6, dpi = 300)

### Filtering and plotting only top 10 orders
# Defining metadata columns
metadata <- order %>% select(sample(1:15)) %>% distinct()

# Calculate total counts of each organism across the entire dataset
total_counts <- order %>%
  group_by(order) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))

# Identify the top 12 organisms
top_12_order <- total_counts %>%
  slice_head(n = 12) %>%
  pull(order)

# View the top 12 organisms
print(top_12_order)

# Label the top 12 organisms and aggregate the rest as "other"
order_labeled <- order %>%
  mutate(order = if_else(order %in% top_12_order, order, "other")) %>%
  group_by(index, order) %>%
  summarise(count = sum(count), .groups = 'drop')

# Ensure the organism column is a factor with correct levels
order_labeled <- order_labeled %>%
  mutate(order = factor(order, levels = c(top_12_order, "other")))

# Merge Back with Metadata
order_labeled <- order_labeled %>%
   left_join(metadata, by = "index")

# Making stacked bar plot - top 12 orders
plot_order_top12 <- ggplot(order_labeled, 
                           aes(x = index, y = count, fill = order)) +
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
plot_order_top12

# Saving last generated plot
ggsave(filename = "figures/taxa_bar_plot_order_top12.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_order_top12.jpg", 
       width = 6.75, height = 4, dpi = 300)
