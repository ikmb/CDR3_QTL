hla_genes = c('A','C','B','DRB1','DQA1', 'DQB1', 'DPA1','DPB1')
smallest_number <- .Machine$double.xmin
our_manova <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/HEALTHY/manova_results_wo_correlation_in_PCA.tsv')
bonf <- 0.01/uniqueN(our_manova$pair)

our_manova <- na.omit(our_manova) %>%
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>%
    rename('Pvalue' = grep('Pr', names(our_manova), value = TRUE)) %>%
    mutate(Site_hla = as.integer(Site_hla), Pvalue = ifelse(Pvalue == 0, smallest_number, Pvalue)) 
bonf <- 0.01/nrow(our_manova)
ishigaki <- read.table('/work_beegfs/sukmb667/projects/cdr3-qtl/cdr3qtl_MVML_amino_acid_geno_all_results.rm_gl.txt', 
    sep = '\t', header = TRUE)

both_datasets <- our_manova %>% dplyr::select(Pvalue,HLA,Site_hla,Length_cdr3,Position_cdr3,variance_explained) %>% 
    mutate(Dataset = 'in_house') %>% 
    rename(Variance_explained = variance_explained) %>%
    rbind(ishigaki %>% mutate(Dataset = 'Ishigaki')) %>% mutate(Site_hla = as.integer(Site_hla))



pl_ours <- ggplot(both_datasets %>% filter(Dataset == 'in_house')) +
    geom_point(aes(x = Site_hla, y = -log10(Pvalue), color = ifelse(Pvalue < bonf, 'significant', ' ')), 
        alpha=0.5, size=1.3, show.legend = FALSE) +
    scale_color_manual(values = c('significant' = '#44AA99', ' ' = 'grey')) +
    #scale_x_discrete(breaks = as.integer(seq(0,300, length.out = 5))) +
    #geom_rug(aes(x = Site_hla, color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    #scale_color_gradient(low = "white", high = "black") +
    labs(x = "HLA sites", y = expression(paste(-log[10], " ", P, "  (Our study)"))) +
    geom_hline(yintercept= - log10(bonf), linetype="dashed", color = "red") +
    geom_hline(yintercept= - log10(smallest_number), linetype="dashed", color = "lightgrey") +
    theme(legend.position = "none") +
    theme_cowplot() +
    facet_grid(~factor(HLA, levels = hla_genes), space = "free") +
    theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(vjust = -1, 
            margin = margin(t = 15)),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12, face = 'bold'),
        axis.ticks.x=element_blank()) 
    #geom_segment(data = hla_annotation, aes(x = start, xend = end, y=-0.01, yend=-0.01), color = 'blue', label = NA)
#ggsave('../plots/manhattan_plot.jpg',pl, width = 10, height = 4)


pl <- ggplot(both_datasets %>% filter(Dataset == 'Ishigaki')) +
    geom_point(aes(x = Site_hla, y = -log10(Pvalue), color=ifelse(Pvalue < (0.01/24360), 'significant', ' ')), 
        alpha=0.5, size=1.3, show.legend = FALSE) +
    scale_color_manual(values = c('significant' = '#88CCEE', ' ' = 'grey')) +
    scale_x_continuous(breaks = as.integer(seq(0, 300, length.out = 5)), position = 'top') +
    #geom_rug(aes(x = Site_hla, color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    #scale_color_gradient(low = "white", high = "black") +
    labs(x = " ", y = expression(paste( -log[10], " ", P, "  (Ishigaki et al.)"))) +
    geom_hline(yintercept= -log10(0.01/24360), linetype="dashed", color = "red") +
    scale_y_reverse() +
    theme(legend.position = "none") +
    theme_cowplot() +
    ylim(300, 0) +
    facet_grid(~factor(HLA, levels = hla_genes), space = "free") +
    theme(axis.text.x = element_blank(), 
        axis.ticks.x=element_blank(), #strip.background = element_rect(fill = "white"), 
        strip.background = element_blank(),
        strip.text.x = element_blank()) #+ element_text(angle = 60, hjust = 1
        
    #geom_segment(data = hla_annotation, aes(x = start, xend = end, y=-0.01, yend=-0.01), color = 'blue', label = NA)


#ggsave('../plots//manhattan_plots_both_one_axis_inverted.jpg', plots_both_one_axis)