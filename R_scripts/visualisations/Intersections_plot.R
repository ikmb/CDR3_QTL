#path_plot_out <- '../plots/dots_across_all_with_lengths_updated_december_colored.jpg'

expanded_hla_df_all <- data.frame()
intersection_df_all <- data.frame()
sizes_df_all <- data.frame()

hla_genes <- c('A', 'C', 'B','DRB1','DQA1','DQB1','DPA1', 'DPB1')

for (g in hla_genes){
    hla_gene <- list(#dIBD_gene_wise = c(sort(unlist(lapply(list.files(paste0('../conditional_analysis/IBD/IBD_downsampled/gene_wise/',g,'/')), 
                                                #function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    dIBD_across_all = c(sort(unlist(lapply(grep('manova',grep(paste0(g, '_'), list.files('..//conditional_analysis/IBD/IBD_downsampled//all_sites_in_downsampled/'), 
                                                                              value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    dIBD_across_all_9PCs = c(sort(unlist(lapply(grep('manova',grep(paste0(g, '_'), list.files('../conditional_analysis/IBD/IBD_downsampled//all_9PCs/'), 
                                                                                   value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    dIBD_across_all_3PCs = c(sort(unlist(lapply(grep('manova',grep(paste0(g, '_'), list.files('../conditional_analysis/IBD/IBD_downsampled//all_3PCs/'), 
                                                                                   value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    #HLH_gene_wise = c(sort(unlist(lapply(list.files(paste0('../conditional_analysis/HEALTHY/gene_wise/',g,'/')), 
                                                   #function(x)as.integer(unlist(str_split(x, '_'))[2]))))),
                    HLH_across_all = c(sort(unlist(lapply(grep('manova',grep(g, list.files('../conditional_analysis/HEALTHY/across_all/'), 
                                                                             value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    HLH_across_all_9PCs = c(sort(unlist(lapply(grep('manova', grep(paste0(g, '_'), list.files('../conditional_analysis/HEALTHY/all_9PCs/'), 
                                                                                   value = TRUE), value = TRUE, invert = TRUE), 
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    HLH_across_all_3PCs = c(sort(unlist(lapply(grep('manova', grep(paste0(g, '_'), list.files('../conditional_analysis/HEALTHY/all_3PCs/'), 
                                                                                   value = TRUE), value = TRUE, invert = TRUE), 
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    #BOTH_gene_wise = c(sort(unlist(lapply(list.files(paste0('../conditional_analysis/v2/',g,'/')), 
                                                #function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    BOTH_across_all = c(sort(unlist(lapply(grep('manova', grep(paste0(g, '_'), list.files('../conditional_analysis/all/'), 
                                                                               value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    BOTH_across_all_3PCs = c(sort(unlist(lapply(grep('manova', grep(paste0(g, '_'), list.files('../conditional_analysis/all_3PCs/'), 
                                                                                    value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))),
                    BOTH_across_all_9PCs = c(sort(unlist(lapply(grep('manova', grep(paste0(g, '_'), list.files('..//conditional_analysis/all_9PCs/'), 
                                                                                    value = TRUE), value = TRUE, invert = TRUE),
                                                function(x) as.integer(unlist(str_split(x, '_'))[2]))))) )
    hla_gene <- lapply(hla_gene, function(x) if (is.null(x)) NA else x)
    melted_hla_df <- na.omit(melt(hla_gene))
    
    expanded_hla_df <- melted_hla_df %>% complete(value, L1) %>% rowwise() %>% mutate(status = ifelse(value %in% hla_gene[[L1]], "yes", "no"), gene = g)
    if (nrow(expanded_hla_df_all) == 0){
      expanded_hla_df_all <- expanded_hla_df  
    }else{
       expanded_hla_df_all <- rbind(expanded_hla_df_all, expanded_hla_df) 
    }
    
    
    intersection_df <- melted_hla_df %>% group_by(value) %>% summarise(intersections = uniqueN(L1)) %>% mutate(gene = g)
    if (nrow(intersection_df_all) == 0){
      intersection_df_all <- intersection_df  
    }else{
       intersection_df_all <- rbind(intersection_df_all, intersection_df) 
    }
    
    sizes_df <- data.frame(L1 = names(hla_gene), number_sites = sapply(hla_gene, function(x) length(x)), gene = g)
    if (nrow(sizes_df_all) == 0){
      sizes_df_all <- sizes_df  
    }else{
       sizes_df_all <- rbind(sizes_df_all, sizes_df) 
    }
}

expanded_hla_df_all_with_intersections  <- expanded_hla_df_all %>%
  group_by(gene, value) %>%
  mutate(intersections = sum(status == "yes")) %>%
  ungroup() %>% 
  filter(intersections > 1)

intersection_df_all <- intersection_df_all %>% filter(intersections > 1)
expanded_hla_df_all_with_intersections <- expanded_hla_df_all_with_intersections %>% 
    filter(status == 'yes') %>%
    mutate(dataset = case_when(
    grepl("^BOTH_", L1) ~ "both",   
    grepl("^HLH_", L1) ~ "healthy", 
    grepl("^dIBD_", L1) ~ "IBD"))

dots <- ggplot()+  
    geom_point(data = expanded_hla_df_all_with_intersections,
                 aes(x = factor(value), y = factor(L1, levels = names(hla_gene)), color = factor(dataset)), 
                 show.legend = FALSE, size = 3) + 
    scale_x_discrete(expand = c(0.2, 0.2)) + 
    scale_color_viridis_d() + 
    new_scale_color() +
    geom_col(data = intersection_df_all,
         aes(y = intersections, x = factor(value)), 
         width = ifelse(intersection_df_all$gene == 'DPA1', 0.1, 0.4), fill = 'grey', alpha = 0.5) +
    #geom_text(data = expanded_hla_df_all_with_intersections,
              #aes(y = intersections, x = factor(value), label = intersections), vjust = 1, size = 5) +
    theme_cowplot()  + 
    labs(y = ' ', x = ' ') +
    facet_wrap(~factor(gene, levels = hla_genes), ncol = 8, scale = 'free_x') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1.2, size = 7, margin = margin(t = 5)),
         strip.background = element_rect(fill = "white", color = 'black'))
   
ggsave(path_plot_out, dots, width = 25, height = 5)
                                                                       
                                                                       
                                                                       