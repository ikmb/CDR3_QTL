alleles_n_variat <- fread('../../Emerson_DeWitt/hla_variation_summary.tsv') %>% rename(HLA = gene, Site_hla = site)
alleles_n_variat_summary <- alleles_n_variat %>% dplyr::select(HLA, Site_hla, n_variat) %>% distinct()

ishigaki <- read.table('../../cdr3qtl_MVML_amino_acid_geno_all_results.rm_gl.txt', sep = '\t', header = TRUE)
bonf <- 0.01/nrow(ishigaki)
ishigaki_manova_with_variat <- ishigaki %>% left_join(alleles_n_variat_summary, by = c('HLA', 'Site_hla')) %>%
    replace_na(list(n_variat = 1))
independent_sites <-list(DRB1 = c(13, 71, 32, 74, 86, 30))
ishigaki_manova_variat_condition <- ishigaki_manova_with_variat %>% group_by(HLA, Site_hla) %>%  mutate(min_Pvalue = min(Pvalue)) %>% 
    rowwise() %>%
    mutate(independent_signal = ifelse(Pvalue == min_Pvalue & Site_hla %in% independent_sites[[HLA]], Site_hla, NA)) %>% ungroup()

pl_ishigaki <- ggplot(data = ishigaki_manova_variat_condition, aes(x = Site_hla, y = -log10(Pvalue))) +
    geom_point(alpha=0.6, size=1.3) +
    geom_rug(aes(color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    scale_color_gradient(low = "white", high = "black") +  # Gradient for n_variat
    geom_label(aes(x = independent_signal, y = -log10(Pvalue), label = independent_signal), 
           hjust=0.55,vjust=-0.4,size=3, show.legend = FALSE) + 
    labs(x = "HLA sites", y = expression(paste(-log[10], " ", P))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_cowplot() +
    facet_grid(~factor(HLA, levels = hla_genes), space = "free") +
    theme(
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
    )
pl_ishigaki

ggsave('../../ishigaki_manova_with_variat_per_site.jpg')

###-------just for DRB1------------
drb1_ishigaki <- ggplot(data = ishigaki_manova_variat_condition %>% filter(HLA == 'DRB1'), aes(x = Site_hla, y = -log10(Pvalue))) +
    geom_point(alpha=0.6, size=1.3) +
    geom_rug(aes(color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    scale_color_gradient(low = "white", high = "black") +  # Gradient for n_variat
    geom_label(aes(x = independent_signal, y = -log10(Pvalue), label = independent_signal), 
           hjust=0.55,vjust=-0.4,size=3, show.legend = FALSE) + 
    labs(x = "HLA DRB1 sites", y = expression(paste(-log[10], " ", P))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_cowplot() +
    theme(
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
    )
drb1_ishigaki
ggsave('../../ishigaki_drb1_with_variat_per_site.jpg')