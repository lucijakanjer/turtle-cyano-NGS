library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77",  
                             "#332288", "#117733","#AA4499", 
                             "#44AA99", "#999933", "#882255", 
                             "#E69F00", "#661100", "#6699CC", 
                             "#999999")

# import data from excel
V6_cyano_level6_mod <- read_excel("data/V6-cyano-level6-mod.xlsx")
View(V6_cyano_level6_mod)

# Group the data frame by the group variable and calculate the sum of each group
group_sums <- V6_cyano_level6_mod %>%
  group_by(sample) %>%
  summarise(GroupSum = sum(count))

# Join the original data frame with the calculated sums
df_transformed <- V6_cyano_level6_mod %>%
  left_join(group_sums, by = "sample") %>%
  mutate(count = count / GroupSum)

# Check the transformed data frame
print(df_transformed)

# Filter samples that contain some values in a group
df_filtered <- df_transformed %>%
  group_by(Age, sample) %>%
  filter(sum(count) > 0) %>%
  ungroup()

# Create the stacked bar chart with filtered data
ggplot(df_filtered, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("N", "R", "P", "juvenile", "sub-adult", "adult")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)

ggsave(filename = "figures/taxa_bar_plot_V6_facets.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets.jpg", width = 6.75, height = 4, dpi = 300)

# Define the groups to exclude
groups_to_exclude <- c("P", "R", "N")
groups_to_exclude2 <- c("N")

# Filter the data to include only desired groups
df_plot <- subset(df_filtered, !(Age %in% groups_to_exclude))
df_plot2 <- subset(df_filtered, !(Age %in% groups_to_exclude2))

genera_plot <- ggplot(df_plot, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("juvenile", "sub-adult", "adult")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
genera_plot

ggsave(filename = "figures/taxa_bar_plot_V6_facets_turtles.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_turtles.jpg", 
       width = 6.75, height = 4, dpi = 300)

taxa6 <- ggplot(df_plot2, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
taxa6

ggsave(filename = "figures/taxa_bar_plot_V6_facets_2.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_2.jpg", 
       width = 6.75, height = 4, dpi = 300)

### Order taxa bar plot - level 4

V6_cyano_level_4_mod <- read_excel("data/V6-cyano-level-4-mod.xlsx")
View(V6_cyano_level_4_mod)

# Group the data frame by the group variable and calculate the sum of each group
group_sums4 <- V6_cyano_level_4_mod %>%
  group_by(sample) %>%
  summarise(GroupSum = sum(count))

# Join the original data frame with the calculated sums
df_transformed4 <- V6_cyano_level_4_mod %>%
  left_join(group_sums, by = "sample") %>%
  mutate(count = count / GroupSum)

# Check the transformed data frame
print(df_transformed4)

# Filter samples that contain some values in a group
df_filtered4 <- df_transformed4 %>%
  group_by(Age, sample) %>%
  filter(sum(count) > 0) %>%
  ungroup()

# Create the stacked bar chart with filtered data
order_plot <- ggplot(df_filtered4, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("N", "R", "P", "juvenile", "sub-adult", "adult")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
order_plot

ggsave(filename = "figures/taxa_bar_plot_V6_facets_order.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_order.jpg", 
       width = 6.75, height = 4, dpi = 300)

# Define the groups to exclude
groups_to_exclude4 <- c("P", "R", "N")
groups_to_exclude4_2 <- c("N")

# Filter the data to include only desired groups
df_plot4 <- subset(df_filtered4, !(Age %in% groups_to_exclude))
df_plot4_2 <- subset(df_filtered4, !(Age %in% groups_to_exclude2))

order_plot_turtle <- ggplot(df_plot4, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("juvenile", "sub-adult", "adult")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
order_plot_turtle

ggsave(filename = "figures/taxa_bar_plot_V6_facets_turtles_order.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_turtles_order.jpg", 
       width = 6.75, height = 4, dpi = 300)

taxa4 <- ggplot(df_plot4_2, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #guides(fill = "none") +
  scale_fill_manual(values = safe_colorblind_palette)+
  facet_grid(.~ factor(Age, levels=c("juvenile", "sub-adult", "adult", "P", "R")),
             scales = "free_x", space="free")+
  coord_cartesian(ylim = c(0, NA), expand = FALSE)
taxa4

ggsave(filename = "figures/taxa_bar_plot_V6_facets_2_order.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_2_order.jpg", 
       width = 6.75, height = 4, dpi = 300)

taxa4 / taxa6

ggsave(filename = "figures/taxa_bar_plot_V6_facets_2_combo.pdf", 
       width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_v6_facets_2_combo.jpg", 
       width = 6.75, height = 8, dpi = 300)
