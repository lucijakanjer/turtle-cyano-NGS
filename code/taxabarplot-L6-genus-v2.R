# Loading required packages
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(patchwork)

# Importing dataset (exported from taxa barplot QIIME2 QZV artifact)
genus_raw <- read.csv(file = "data/taxabarplots-cyano-V6-level-6-genus.csv", 
                      header = TRUE)

# Creating working table version
genus_wide <- genus_raw

# Identify columns that start with "k__Bacteria"
bacteria_columns <- grep("^k__Bacteria", names(genus_wide), value = TRUE)

# Compute the sum for each row across the identified columns
genus_wide <- genus_wide %>%
  rowwise() %>%
  mutate(Sum_Bacteria = sum(c_across(all_of(bacteria_columns)), na.rm = TRUE)) %>%
  ungroup()

# Check the DataFrame to ensure the new column is added
colnames(genus_wide)

# Step 1: Identify columns that contain "Nodosilineales"
nodosilineales_columns <- grep("Nodosilineales", names(genus_wide), value = TRUE)

# Step 2: Identify all other bacterial columns starting with "k__Bacteria" but NOT containing "Nodosilineales"
bacteria_columns <- grep("^k__Bacteria", names(genus_wide), value = TRUE)
other_columns <- setdiff(bacteria_columns, nodosilineales_columns)

# Step 3: Sum the other bacteria columns and store the result in a new column named 'other'
genus_wide <- genus_wide %>%
  rowwise() %>%
  mutate(other = sum(c_across(all_of(other_columns)), na.rm = TRUE)) %>%
  ungroup()

# Step 4: Remove the original bacterial columns that don't contain "Nodosilineales"
genus_wide <- genus_wide %>%
  select(-all_of(other_columns))

# Step 5: Verify the result
colnames(genus_wide)

colnames(genus_wide)[2] <- "Salileptolyngbya"  # Column 1 (now index 2) to Salileptolyngbya
colnames(genus_wide)[4] <- "Leptothoe"         # Column 4 to Leptothoe
colnames(genus_wide)[5] <- "Rhodoploca"        # Column 5 to Rhodoploca
colnames(genus_wide)[8] <- "Cymatolege"        # Column 8 to Cymatolege
colnames(genus_wide)[10] <- "Euryhalinema"     # Column 10 to Euryhalinema

# Step 2: Collapse unclassified Nodosilineales columns into a single column
unclassified_columns <- c(3, 6, 7, 9, 11)  # Indices of columns to collapse

# Sum up the values of these columns and create a new "unclassified_Nodosilineales" column
genus_wide <- genus_wide %>%
  rowwise() %>%
  mutate(unclassified_Nodosilineales = sum(c_across(all_of(unclassified_columns)), na.rm = TRUE)) %>%
  ungroup()

# Step 3: Remove the original unclassified Nodosilineales columns
genus_wide <- genus_wide %>%
  select(-all_of(unclassified_columns))

# Step 4: Verify the new column names
colnames(genus_wide)

# Converting the table in long format for plotting
cyano_cols <- names(genus_wide)[c(2,3,4,5,6,22,23)]
genus_long <- genus_wide %>%
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

# Ensure the family column is a factor with correct levels
genus_long <- genus_long %>%
  mutate(genus = factor(genus, levels = c("Rhodoploca", "Leptothoe", "Cymatolege",
                                          "Salileptolyngbya", "Euryhalinema",
                                          "unclassified_Nodosilineales","other")))

plot_genus <- ggplot(genus_long, aes(x = index, y = count, fill = genus)) +
  geom_bar(stat = "identity",, position = "fill") +
  labs(x = "Samples", y = "Relative abundance", tag = "B") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_viridis_d(option = "D", direction = -1) +
  facet_grid(.~ factor(Age, 
                       levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
plot_genus

# Saving last generated plot
ggsave(filename = "figures/taxa_bar_plot_genus_Nodosilineales.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_genus_Nodosilineales.jpg", 
       width = 6.75, height = 4, dpi = 300)
ggsave(filename = "figures/taxa_bar_genus_Nodosilineales.tiff", 
       width = 6.75, height = 4, dpi = 300)

# To make combination plot, first run "taxabarplot-L4-order.R" script 
# object "plot_order_top12" must stay in the environment
plot_order_top6_viridis / plot_genus

ggsave(filename = "figures/taxa_bar_plot_comb_viridis.pdf", 
       width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_comb_viridis.jpg", 
       width = 6.75, height = 8, dpi = 300)
ggsave(filename = "figures/taxa_bar_plot_comb_viridis.tiff", 
       width = 6.75, height = 8, dpi = 300)
