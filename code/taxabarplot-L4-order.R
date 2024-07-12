# Loading required packages
library(tidyverse)
library(RColorBrewer)
library(patchwork)

safe_colorblind_palette <- c("#88CCEE","#CC6677","#DDCC77","#332288","#117733",
                             "#AA4499","#44AA99","#999933","#882255","#E69F00",
                             "#661100","#6699CC","#999999")
                             

# Importing dataset (exported from taxa barplot QIIME2 QZV artifact)
order_raw <- read.csv(file = "data/taxabarplots-cyano-V6-level-4-order.csv", 
                      header = TRUE)

# Creating working table version - wide format
order_wide <- order_raw

# Visualizing original column names
print(colnames(order_raw))

# Changing column names to keep only short names for Orders
for(i in 1:length(colnames(order_wide))) {
  if(grepl("k__Bacteria.p__Cyanobacteriota.", 
           colnames(order_wide)[i])) {
    colnames(order_wide)[i] <- gsub("k__Bacteria.p__Cyanobacteriota.", "",
                                   colnames(order_wide)[i])
  }
}

for(i in 1:length(colnames(order_wide))) {
  if(grepl("c__Cyanophyceae.o__", 
           colnames(order_wide)[i])) {
    colnames(order_wide)[i] <- gsub("c__Cyanophyceae.o__", "",
                               colnames(order_wide)[i])
  }
}

order_wide <- order_wide %>%
  rename(unclassified_Cyanophyceae = "c__Cyanophyceae.__", 
         Caenarcaniphilales = "c__Vampirivibrionia.o__Caenarcaniphilales",
         unclassified_Cyanobacteriota = "__.__",
         unclassified_Sericytochromatia = "c__Sericytochromatia.o__")

# Checking the changed names
print(colnames(order_wide))

# Converting the table in long format for plotting
cyano_cols <- names(order_wide)[2:20]
order_long <- order_wide %>%
  pivot_longer(
    cols = all_of(cyano_cols),
    names_to = "order",
    values_to = "count"
  )

# Changing sample names
order_long <- order_long %>%
  mutate(index = str_remove(index, "V6-"))

# Changing names
order_long <- order_long %>%
  mutate(Age = case_when(
    Age == "rocks" ~ "R",
    Age == "pool" ~ "P",
    TRUE ~ Age  # Keep other values unchanged
  ))

# Removing negative control sample
order_long <- order_long %>%
  filter(index != "negctrl")


# Making stacked bar plot - all orders
plot_order <- ggplot(order_long, aes(x = index, y = count, fill = order)) +
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
metadata <- order_long %>% select(sample(1:15)) %>% distinct()

# Calculate total counts of each organism across the entire dataset
total_counts <- order_long %>%
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
order_labeled <- order_long %>%
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

### Order ratios across dataset

# Ensure the first column is treated as row names (if needed)
order_wide <- order_wide %>% column_to_rownames(var = "index")

# Calculate the total abundance for each sample
order_wide <- order_wide %>%
  rowwise() %>%
  mutate(total_abundance = sum(c_across(all_of(cyano_cols))))

# Calculate the ratio of each species' abundance for each sample
order_data_ratio <- order_wide %>%
  mutate(across(all_of(cyano_cols), ~ . / total_abundance))

# Drop the total_abundance column
order_data_ratio <- order_data_ratio %>% select(-total_abundance)

# Convert the rowwise data frame back to a regular data frame for summarization
order_data_ratio <- ungroup(order_data_ratio)

# Calculate the median, maximum, and minimum ratios for each species
order_statistics <- order_data_ratio %>%
  summarise(across(all_of(cyano_cols), list(
    median = ~ median(.x, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE)
  ), .names = "{col}__{fn}"))

# Pivot the data to a long format
order_statistics_long <- order_statistics %>%
  pivot_longer(
    cols = everything(),
    names_to = c("order", ".value"),
    names_sep = "__"
  )

view(order_statistics_long)

# Calculate the total abundance for each order across all samples
total_abundance_order <- colSums(order_wide[ , cyano_cols])

# Calculate the total abundance of all orders across all samples
total_abundance_order_all <- sum(total_abundance_order)

# Calculate the ratio for each order
order_ratio <- total_abundance_order / total_abundance_order_all

# Convert the order ratios to a data frame
order_ratio_df <- data.frame(
  sum = order_ratio)

# Ensure order names are consistent
order_ratio_df <- order_ratio_df %>%
  mutate(Order = rownames(order_ratio_df))

# Reset rownames
rownames(order_ratio_df) <- NULL

# Combine the data frames using left_join
order_combined_df <- order_statistics_long %>%
  left_join(order_ratio_df, by = c("order" = "Order"))

# Arrange the combined data frame by descending values of "sum"
order_combined_df <- order_combined_df %>%
  arrange(desc(sum))

# View the combined data frame
view(order_combined_df)

# Export the combined data frame to a CSV file
write.csv(order_combined_df, "order_combined_data.csv", row.names = FALSE)
