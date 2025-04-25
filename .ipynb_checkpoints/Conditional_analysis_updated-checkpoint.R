
# hla_matrices <- 
#dir_analysis <- '../conditional_analysis/all_9PCs/' 
#pairs_to_test <- '../cdr3_hla_pairs_9PCs/'
#pcs <- 9
# manova_cond_all <- fread('../manova_results_9PCs.tsv')

first <- TRUE
while (nrow(manova_cond_all %>% filter(`Pr(>F)` < bonf)) > 0){
    manova_results <- manova_cond_all %>% 
        filter(model == 'full') %>%
        arrange(., `Pr(>F)`,desc(variance_explained))
     
    hit <- manova_results$pair[1]
    g_s_l_p <- unlist(strsplit(hit, '_'))
    condition_site <- paste(g_s_l_p[1:2], collapse= '_')
    condition_site_df <- fread(paste0(hla_matrices, grep(paste0(condition_site, '_'), list.files(hla_matrices), value = TRUE))) %>%
        dplyr::select(-G_S) %>% dplyr::select(!contains('PC'))
    colnames(condition_site_df) <- gsub('allele', condition_site, colnames(condition_site_df))
    if(first){
        covariates_df <- condition_site_df
        first <- FALSE
        } else {
        covariates_df <- inner_join(covariates_df, condition_site_df)
        }
    
    sites <- colnames(covariates_df)[-1] 
    
    hit_pattern <- paste0(condition_site,'_')
    #################
    pairs_to_test <- grep(hit_pattern, pairs_to_test, value = TRUE, invert = TRUE)
    
    dir_results <- paste0(dir_analysis, condition_site,'/')
        
    if (!file.exists(dir_results)) {
    dir.create(dir_results, recursive = TRUE)
    }
    
    first <- TRUE
    for (f in pairs_to_test){
        
        cdr3_hla_matrix <- na.omit(fread(f)) %>% left_join(covariates_df, by = 'patient_id')
        pair_name <- cdr3_hla_matrix$pair[1]
        
        manova_df <- mlm_fun_conditional_PCs(cdr3_hla_matrix, dir_results, sites, n_pcs = pcs) 
        if (first){
            manova_cond_all <- manova_df[0,]
            first <- FALSE
        }
        
        if (is.null(manova_df)){
            next
        } else {

            manova_cond_all <- rbind(manova_cond_all, manova_df, fill=TRUE)
            
        }
    }
    manova_cond_all <- manova_cond_all %>% separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = '_', remove = FALSE) %>% 
        mutate(condition = condition_site)
    fwrite(manova_cond_all, paste0(dir_analysis, 'manova_',condition_site,'.tsv'), sep = '\t')
        
}


                                    ###-------combine manova_df_site together--------
manova_files <- list.files(path = dir_analysis, pattern = "*.tsv", full.names = TRUE)
manova_df_all <- rbindlist(lapply(manova_files, fread))
fwrite(manova_df_all, paste0(dir_analysis, 'manova_all_',pcs,'PCs.tsv'), sep = '\t')