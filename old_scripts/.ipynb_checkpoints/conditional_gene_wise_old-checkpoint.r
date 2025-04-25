main_manova <- na.omit(fread('../IBD/manova_results_downsampled.tsv')) %>% 
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
    arrange(., Pr..F., desc(variance_explained))
bonf <- 0.01 / nrow(main_manova)
hla_genes <- c('DRB1','A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1','DQB1')
cond_manova_template <- read_tsv('../conditional_analysis//v2//manova_template.tsv')[0,] %>% 
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE)

for (gene in hla_genes){
    
    manova_df_all <- cond_manova_template
    
    gene_manova <- main_manova %>% filter(HLA == gene) 
    hit <- gene_manova$pair[1]
    g_s_l_p <- unlist(strsplit(hit, ':'))
    condition_site <- paste(g_s_l_p[1:2], collapse= '_')
    condition_site_df <- fread(paste0('../hla_matrices/all_variable_sites/', grep(paste0(condition_site, '_'), list.files('../hla_matrices/all_variable_sites/'), value = TRUE))) %>%
        dplyr::select(-G_S)
    colnames(condition_site_df) <- gsub('allele', condition_site, colnames(condition_site_df))
    sites <- colnames(condition_site_df)[-1]
    hit_pattern <- paste0(gsub('_', ':', condition_site),':')
    pairs_to_test <- grep(hit_pattern, paste0('../IBD/cdr3_hla_pairs_downsampled/',grep(paste0(gene,':'), list.files('../IBD/cdr3_hla_pairs_downsampled/'), value = TRUE)), value = TRUE, invert = TRUE)
    dir_matrices <- paste0('../conditional_analysis//IBD_downsampled/',gene,'/', condition_site,'/cdr3_hla_pairs/')
    dir_results <- paste0('../conditional_analysis//IBD_downsampled/',gene,'/', condition_site,'/')
        
    if (!file.exists(dir_matrices)) {
    dir.create(dir_matrices, recursive = TRUE)
    }
    if (!file.exists(dir_results)) {
    dir.create(dir_results, recursive = TRUE)
    }
    manova_df_site <- cond_manova_template
    
    for (f in pairs_to_test){
        cdr3_hla_condition_matrix <- fread(f) %>% left_join(condition_site_df, by = 'patient_id')
        cdr3_hla_pair_name <- cdr3_hla_condition_matrix$pair[1]
        fwrite(cdr3_hla_condition_matrix, paste0(dir_matrices,cdr3_hla_pair_name,':',condition_site,'.tsv'), sep = '\t')
        manova_df <- mlm_fun_conditional(cdr3_hla_condition_matrix, dir_results, sites) %>% 
            separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
            mutate(condition = condition_site)
        manova_df_site <- rbind(manova_df_site, manova_df)
        
    }
    fwrite(manova_df_site, paste0(dir_results, 'manova_',condition_site,'.tsv'))
    manova_df_all <- rbind(manova_df_all, manova_df_site)

    while (nrow(manova_df_site %>% filter(Pr..F.<bonf)) > 0){
        manova_results <- na.omit(fread(paste0('../conditional_analysis/IBD_downsampled/',gene,'/',condition_site,'/manova_',condition_site,'.tsv'))) %>% 
            arrange(., Pr..F.,desc(variance_explained))
         
        hit <- manova_results$pair[1]
        g_s_l_p <- unlist(strsplit(hit, ':'))

        condition_site <- paste(g_s_l_p[1:2], collapse= '_')
        hit_pattern <- paste0(gsub('_', ':', condition_site),':')
        condition_site_df <- fread(paste0('../hla_matrices/all_variable_sites/', grep(paste0(condition_site, '_'), list.files('../hla_matrices/all_variable_sites/'), value = TRUE))) %>% 
            dplyr::select(-G_S)
        colnames(condition_site_df) <- gsub('allele', condition_site, colnames(condition_site_df))
        site <- colnames(condition_site_df)[-1]
        sites <- c(sites, site)
        pairs_survived <- unique(lapply(unique(manova_results %>% filter(Pr..F. < bonf) %>% pull(pair)), function(x) paste0(paste(unlist(strsplit(x, ':'))[1:2], collapse = ':'), ':') ))
        pairs_to_test <- grep(hit_pattern, paste0(dir_matrices, grep(paste(pairs_survived, collapse = '|'), list.files(dir_matrices), value = TRUE)), value = TRUE, invert = TRUE)
    
        dir_matrices <- paste0('../conditional_analysis//IBD_downsampled/',gene,'/', condition_site,'/cdr3_hla_pairs/')
        dir_results <- paste0('../conditional_analysis//IBD_downsampled/',gene,'/', condition_site,'/')
            
        if (!file.exists(dir_matrices)) {
        dir.create(dir_matrices, recursive = TRUE)
        }
        if (!file.exists(dir_results)) {
        dir.create(dir_results, recursive = TRUE)
        }
        manova_df_site <- cond_manova_template
        
        for (f in pairs_to_test){
            cdr3_hla_condition_matrix <- fread(f) %>% left_join(condition_site_df, by = 'patient_id')
            cdr3_hla_pair_name <- cdr3_hla_condition_matrix$pair[1]
            fwrite(cdr3_hla_condition_matrix, paste0(dir_matrices,cdr3_hla_pair_name,':',condition_site,'.tsv'), sep = '\t')
            manova_df <- mlm_fun_conditional(cdr3_hla_condition_matrix, dir_results, sites) %>% 
                separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
                mutate(condition = condition_site)
            manova_df_site <- rbind(manova_df_site, manova_df)
            
        }
        fwrite(manova_df_site, paste0(dir_results, 'manova_',condition_site,'.tsv'), sep = '\t')
        manova_df_all <- rbind(manova_df_all, manova_df_site)
            
    }
    fwrite(manova_df_all, paste0('../conditional_analysis//IBD_downsampled/manova_',gene,'.tsv'), sep = '\t')
    }
