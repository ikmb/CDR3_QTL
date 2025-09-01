
hla_genes <- c('A','C','B','DRB1','DQA1', 'DQB1', 'DPA1','DPB1')
smallest_number <- .Machine$double.xmin
#our_manova <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/HEALTHY/manova_results_wo_correlation_in_PCA.tsv')
our_manova <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/main_manova_as_in_Ishigaki_removing_90_corr.tsv')
bonf <- 0.01/uniqueN(our_manova$pair)

our_manova <- na.omit(our_manova) %>%
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = '_', remove = FALSE) %>%
    rename('Pvalue' = grep('Pr', names(our_manova), value = TRUE)) %>%
    mutate(Site_hla = as.integer(Site_hla), Pvalue = ifelse(Pvalue == 0, smallest_number, Pvalue)) %>% 
    as.data.table()
bonf <- 0.01/nrow(our_manova)
ishigaki <- read.table('/work_ikmb/sukmb667/projects/cdr3-qtl/cdr3qtl_MVML_amino_acid_geno_all_results.rm_gl.txt', 
    sep = '\t', header = TRUE) %>% as.data.table()

both_datasets <- our_manova %>% dplyr::select(Pvalue,HLA,Site_hla,Length_cdr3,Position_cdr3,variance_explained) %>% 
    mutate(Dataset = 'in_house') %>% 
    rename(Variance_explained = variance_explained) %>%
    rbind(ishigaki %>% mutate(Dataset = 'Ishigaki')) %>% mutate(Site_hla = as.integer(Site_hla))

unique_sites_from_both_analyses <- both_datasets[Pvalue < bonf][, .(all_sites = list(unique(Site_hla))), by = .(HLA, Dataset)]
unique_sites_from_both_analyses <- unique(both_datasets[Pvalue < bonf][, .(HLA, Site_hla, Dataset)])
first <- TRUE
for (hla in hla_genes){
    dt <- unique_sites_from_both_analyses[HLA == hla]
    unique_sites <- fsetdiff(dt[Dataset == 'in_house', .(Site_hla)], dt[Dataset == 'Ishigaki', .(Site_hla)] )
    unique_sites[, HLA := hla]
    if (first){
        unique_sites_all <- unique_sites
        first <- FALSE
    } else {
        unique_sites_all <- rbind(unique_sites_all, unique_sites)
    }
}

unique_sites_all_annotated <- hla_pocket(unique_sites_all)
unique_sites_all_annotated[, novelty := "novel"][, Dataset := "in_house"]
both_datasets_annotated <- merge(unique_sites_all_annotated, both_datasets, by = c('HLA', 'Site_hla', "Dataset"), all.y = TRUE)

novel_with_pval <- both_datasets_annotated %>% filter(Dataset == 'in_house', Pvalue < bonf, novelty == 'novel') %>% filter(!is.na(Pocket)) %>% 
    group_by(HLA, Site_hla) %>% filter(Pvalue == min(Pvalue)) %>% slice(1) %>% ungroup() %>% mutate(start = ifelse(HLA == "B", 150, 80))

pl_ours <- ggplot(both_datasets %>% filter(Dataset == 'in_house')) +
    geom_point(aes(x = Site_hla, y = -log10(Pvalue), color = ifelse(Pvalue < bonf, 'significant', ' ')), 
        alpha=0.5, size=1.3, show.legend = FALSE) +
    scale_color_manual(values = c('significant' = 'orange', ' ' = 'grey')) +
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

pl_our <- ggplot(both_datasets_annotated %>% filter(Dataset == 'in_house')) +
    geom_point(aes(x = Site_hla, y = -log10(Pvalue), color = ifelse(Pvalue < bonf, 'significant', ' ')), 
        alpha=0.5, size=1.5, show.legend = FALSE) +
    scale_color_manual(values = c('significant' = '#fc8961', ' ' = 'grey')) +
    new_scale_color() +
    geom_point(data = both_datasets_annotated %>% filter(Dataset == 'in_house', Pvalue < bonf, novelty == 'novel'), 
        aes(x = Site_hla, y = -log10(Pvalue), color = novety, shape = ifelse(is.na(Pocket), ' ', "ABP")), 
        alpha=0.5, size=1.5, show.legend = FALSE, colour = "violet") +
    #scale_color_manual(values = c('novel' = 'violet', ' ' = 'grey')) +
    #scale_x_discrete(breaks = as.integer(seq(0,300, length.out = 5))) +
    #geom_rug(aes(x = Site_hla, color = n_variat), sides = "b", show.legend = TRUE) +  # Use n_variat for the color aesthetic
    #scale_color_gradient(low = "white", high = "black") +
    labs(x = "HLA allelic variants (sites)", y = expression(paste(-log[10], " ", P, "  (Current study)"))) +
    geom_hline(yintercept= - log10(bonf), linetype="dashed", color = "salmon") +
    geom_hline(yintercept= - log10(smallest_number), linetype="dashed", color = "lightgrey") +
    theme(legend.position = "none") +
    theme_cowplot() +
    facet_grid(~factor(HLA, levels = hla_genes), space = "free") +
    theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(vjust = -1, 
            margin = margin(t = 15)),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 12, face = 'bold'),
        axis.ticks.x=element_blank()) +
    geom_segment(data = novel_with_pval, 
        aes(x = Site_hla, xend = Site_hla, y= start, yend= start - 20), color = '#b73779', alpha = 0.5,
        arrow = arrow(length = unit(0.2, "cm"), type = "closed"), 
        size = 0.5, show.legend = FALSE) #+ geom_text(data = na.omit(unique_sites_all_annotated), aes(x = Site_hla, y=150, label = Pocket))
#ggsave('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/plots/manhattan_ours_with_novel.png', pl_our, dpi = 600, width = 10, height = 5)
#ggsave('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/plots/manhattan_ours_with_novel.pdf', pl_our, dpi = 600, width = 10, height = 5)

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