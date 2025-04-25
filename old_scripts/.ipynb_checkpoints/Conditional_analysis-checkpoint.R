source('cdr3-QTL_functions.R')
source('libraries.R')

#hla_features <- fread('../data/hla_features_wo_duplicates_healthy_and_ibd_those_with_cdr3.tsv')
#main_manova <- fread('../manova_results_9PCs.tsv')
#dir_analysis <- '../conditional_analysis/all_9PCs/' 
#dir_cdr3_hla_pairs <- '../cdr3_hla_pairs_9PCs/'
#pcs <- 9

main_manova <- na.omit(main_manova) %>% 
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
    arrange(., Pr..F., desc(variance_explained))

bonf <- 0.01 / nrow(main_manova)

cond_manova_template <- read_tsv('../conditional_analysis//v2//manova_template.tsv')[0,] %>% 
    separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE)

  
manova_df_all <- cond_manova_template

hit <- main_manova$pair[1]
g_s_l_p <- unlist(strsplit(hit, ':'))
condition_site <- paste(g_s_l_p[1:2], collapse= '_')
condition_site_df <- fread(paste0('../hla_matrices_v3/wo_correlation//', grep(paste0(condition_site, '_'), list.files('../hla_matrices_v3/wo_correlation//'), value = TRUE))) %>%
    dplyr::select(-G_S)
colnames(condition_site_df) <- gsub('allele', condition_site, colnames(condition_site_df))
sites <- colnames(condition_site_df)[-1]
hit_pattern <- paste0(gsub('_', ':', condition_site),':')

pairs_to_test <- grep(hit_pattern,paste0(dir_cdr3_hla_pairs, 
                                         list.files(dir_cdr3_hla_pairs)), value = TRUE, invert = TRUE)

dir_matrices <- paste0(dir_analysis, condition_site,'/cdr3_hla_pairs/')
dir_results <- paste0(dir_analysis, condition_site,'/')
    
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
    manova_df <- mlm_fun_conditional_PCs(cdr3_hla_condition_matrix, dir_results, sites, n_pcs = pcs) %>% 
        separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
        mutate(condition = condition_site)
    manova_df_site <- rbind(manova_df_site, manova_df)
    
}
fwrite(manova_df_site, paste0(dir_results, 'manova_',condition_site,'.tsv'))
#backup
fwrite(manova_df_site, paste0(dir_analysis, 'manova_',condition_site,'.tsv'))

while (nrow(manova_df_site %>% filter(Pr..F.<bonf)) > 0){
    manova_results <- na.omit(fread(paste0(dir_analysis,condition_site,'/manova_',condition_site,'.tsv'))) %>% 
        arrange(., Pr..F.,desc(variance_explained))
     
    hit <- manova_results$pair[1]
    g_s_l_p <- unlist(strsplit(hit, ':'))

    condition_site <- paste(g_s_l_p[1:2], collapse= '_')
    hit_pattern <- paste0(gsub('_', ':', condition_site),':')
    condition_site_df <- fread(paste0('../hla_matrices_v3//wo_correlation/', grep(paste0(condition_site, '_'), 
                                                                                  list.files('../hla_matrices_v3//wo_correlation/'), value = TRUE))) %>% 
        dplyr::select(-G_S)
    colnames(condition_site_df) <- gsub('allele', condition_site, colnames(condition_site_df))
    site <- colnames(condition_site_df)[-1]
    sites <- c(sites, site)
    pairs_survived <- unique(lapply(unique(manova_results %>% 
                                           filter(Pr..F. < bonf) %>% 
                                           pull(pair)), 
                                    function(x) paste0(paste(unlist(strsplit(x, ':'))[1:2], collapse = ':'), ':') ))
    pairs_to_test <- grep(hit_pattern, paste0(dir_matrices, grep(paste(pairs_survived, collapse = '|'), 
                                                                 list.files(dir_matrices), value = TRUE)), value = TRUE, invert = TRUE)

    dir_matrices <- paste0(dir_analysis, condition_site,'/cdr3_hla_pairs/')
    dir_results <- paste0(dir_analysis, condition_site,'/')
        
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
        manova_df <- mlm_fun_conditional_PCs(cdr3_hla_condition_matrix, dir_results, sites, n_pcs = pcs) %>% 
            separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = ':', remove = FALSE) %>% 
            mutate(condition = condition_site)
        manova_df_site <- rbind(manova_df_site, manova_df)
        
    }
    fwrite(manova_df_site, paste0(dir_results, 'manova_',condition_site,'.tsv'))
    fwrite(manova_df_site, paste0(dir_analysis, 'manova_',condition_site,'.tsv'), sep = '\t')
        
}


                                    ###-------combine manova_df_site together--------
manova_files <- list.files(path = dir_analysis, pattern = "*.tsv", full.names = TRUE)
manova_df_all <- rbindlist(lapply(manova_files, fread))
fwrite(manova_df_all, paste0(dir_analysis, 'manova_all_',pcs,'PCs.tsv'), sep = '\t')