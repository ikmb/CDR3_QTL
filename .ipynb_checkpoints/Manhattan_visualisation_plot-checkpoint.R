
#manova_df <- fread()
hla_genes <- c('A', 'C', 'B','DRB1','DQA1','DQB1','DPA1', 'DPB1')
bonf <- 0.01/uniqueN(manova_df$pair)
alleles_n_variat_summary <- fread('../data/hla/hla_site_variations_healthy_and_ibd_all_variable_sites_v3.tsv')
alleles_n_variat_summary <- distinct(alleles_n_variat_summary[, .(Site_hla, HLA, n_variat)])




manova_with_variat <- manova_df %>% left_join(alleles_n_variat_summary, by = c('HLA', 'Site_hla')) %>%
    replace_na(list(n_variat = 1))


pl <- ggplot(data = manova_with_variat, aes(x = Site_hla, y = -log10(`Pr(>F)`))) +
    geom_point(alpha=0.6, size=1.3) +
    geom_rug(aes(x = Site_hla, color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    scale_color_gradient(low = "white", high = "black") +  # Gradient for n_variat
    #geom_label(aes(x = independent_signal, y = -log10(`Pr(>F)`), label = independent_signal), 
           #hjust=0.55,vjust=-0.4,size=3, show.legend = FALSE) + 
    labs(x = "HLA sites", y = expression(paste(-log[10], " ", P))) +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", color = "red") +
    theme_cowplot() +
    facet_grid(~factor(HLA, levels = hla_genes), space = "free") +
    theme(
        axis.text.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
    )
pl

#ggsave('../plots/manova_with_variat_per_site.jpg', pl)
