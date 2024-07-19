
# Load required libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(patchwork)
library(dplyr)


# Read the metadata.tsv file
metadata <- read.delim("data/metadata-V6-v4.tsv")

# Subset only samples where SampleType is "carapace"
carapace_samples <- metadata %>%
  filter(SampleType == "carapace")
print(carapace_samples)

# Load world map data
world <- ne_countries(scale = "large", returnclass = "sf")
countries <- c("Croatia", "Slovenia", "Hungary", "Bosnia and Herz.", "Italy")

croatia_map <- world[world$name %in% countries, ]


# Create map plot
map_adriatic <- ggplot() +
  geom_sf(data = croatia_map, fill = "gray", color = NA) +
  coord_sf(xlim = c(13, 18), ylim = c(42, 46), expand = FALSE) + # Set limits for the plot
  theme_minimal() +
  geom_point(data = carapace_samples, 
             aes(x = lon, y = lat, color = Age, shape = Age),
             size = 3) +
  scale_shape_discrete(name  ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  scale_color_brewer(palette = "Set2", name ="Turtle Age", breaks = c("juvenile", "sub-adult", "adult")) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  labs(x = "", y = "", tag = "B") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.85,0.85),
        legend.background = element_rect(fill = "white", color = NA))  

map_adriatic

# Saving plot in JPG and PDF format
ggsave(filename = "figures/map_adriatic.pdf", 
       width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/map_adriatic.jpg", 
       width = 6.75, height = 8, dpi = 300)


# make Europe map

# Filter data for Europe
europe_map <- world[world$continent == "Europe", ]

# Create base plot
map_europe <- ggplot() +
  geom_sf(data = europe_map, fill = "gray", color = NA) +
  coord_sf(xlim = c(-10, 25), ylim = c(35, 60), expand = FALSE) + 
  theme_minimal() +
  geom_rect(
    aes(xmin = 13, xmax = 18, ymin = 42, ymax = 46),
    fill = NA,
    color = "black",
    size = 1) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  labs(x = "", y = "", tag = "A") +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank())  # Remove y-axis text

map_europe

# Saving plot in JPG and PDF format
ggsave(filename = "figures/map_europe.pdf", 
       width = 6.75, height = 8, device = cairo_pdf)
ggsave(filename = "figures/map_europe.jpg", 
       width = 6.75, height = 8, dpi = 300)

map_europe + map_adriatic

# Saving plot in JPG and PDF format
ggsave(filename = "figures/map_combined.pdf", 
       width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/map_combined.jpg", 
       width = 6.75, height = 4, dpi = 300)

