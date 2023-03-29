# Combined V34 and V6 figures


order_v34 + order_v6
ggsave(filename = "figures/taxa_bar_plot_order_combo.pdf", width = 6.75, height = 4, device = cairo_pdf)
ggsave(filename = "figures/taxa_bar_plot_order_combo.jpg", width = 6.75, height = 4, dpi = 300)


cyano_v34_genus_collapse_otu_table <- data.frame(otu_table(cyano_v34_genus_collapse))
cyano_v34_genus_collapse_otu_table$Genus<-rownames(cyano_v34_genus_collapse_otu_table)

cyano_v6_genus_collapse_otu_table <- data.frame(otu_table(cyano_v6_genus_collapse))
cyano_v6_genus_collapse_otu_table$Genus<-rownames(cyano_v6_genus_collapse_otu_table)



cyano_genus_combined <- full_join(cyano_v34_genus_collapse_otu_table,cyano_v6_genus_collapse_otu_table, by = "Genus")


cyano_genus_combined_plot <-ggplot(cyano_genus_combined, aes(fill=Genus), x=TB159, y=1) + 
  geom_bar(position="fill", stat="identity") +  
  theme_classic() 
cyano_genus_combined_plot

test_cyano_merge_taxa <- merge_taxa(cyano_v34, cyano_v6)

test_cyano_merge_taxa <- merge_taxa(cyano_v34_top10genera_filter, cyano_v6_top10genera_filter)
test_cyano_merge_samples <- merge_samples(cyano_v34_top10genera_filter, cyano_v6_top10genera_filter)

test_cyano_merge_phyloseq <- merge_phyloseq(cyano_v34_top10genera_filter, cyano_v6_top10genera_filter)

test_cyano_merge_taxa_plot <- test_cyano_merge_taxa %>%
  aggregate_taxa(level = "Genus") %>%  
  microbiome::transform(transform = "compositional") %>%
  plot_composition(sample.sort = "CCL")
test_cyano_merge_taxa_plot + scale_fill_brewer(palette="Set3")

plot_composition(test_cyano_merge_taxa,
                 sample.sort = CCL,
                 A)
