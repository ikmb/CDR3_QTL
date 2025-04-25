source('libraries.R')
source('cdr3-QTL_functions.R')


#hla_features <- read_tsv('../data/hla_features_wo_duplicates_healthy_and_ibd_those_with_cdr3.tsv')
#pairs <- '../cdr3_hla_pairs_v3/'
#dir_pair_path <- paste0('../cdr3_hla_pairs_9PCs/')
#dir_results <- paste0('../manova_results_9PCs/')
#pcs <- 9

pairs_paths <- paste0(pairs, list.files(pairs))

manova_df_all <- fread('../manova_results_v3.tsv')[0,]

manova_df_all_3PCs <- fread('../manova_results_v3.tsv')[0,]
    
if (!file.exists(dir_pair_path)) {
dir.create(dir_pair_path, recursive = TRUE)
}
if (!file.exists(dir_results)) {
dir.create(dir_results, recursive = TRUE)
}
if (!file.exists(dir_results_3PCs)) {
dir.create(dir_results_3PCs, recursive = TRUE)
}


for (f in pairs_paths){
    
    cdr3_hla_matrix <- read_tsv(f) 
    ind_in_pair <- unique(cdr3_hla_matrix_prev$patient_id)
    
    pca_hla <- pca_hla_fun(hla_features %>% filter(patient_id %in% ind_in_pair))
    
    cdr3_hla_matrix <- cdr3_hla_matrix_prev %>%
        inner_join(pca_hla, by = 'patient_id')
    name_pair <- cdr3_hla_matrix$pair[1]
    
    write_tsv(cdr3_hla_matrix, paste0(dir_pair_path,name_pair,'_matrix.tsv'))
    
    manova_results <- mlm_fun(cdr3_hla_matrix, dir_results, n_pcs = pcs)
    manova_df_all <- rbind(manova_df_all, manova_results)
        
        
    }


fwrite(manova_df_all, paste0('../manova_results_',pcs,'PCs.tsv', sep = '\t')