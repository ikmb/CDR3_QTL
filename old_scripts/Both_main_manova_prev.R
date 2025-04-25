source('libraries.R')
source('cdr3-QTL_functions.R')


hla_features <- read_tsv('../data/hla_features_wo_duplicates_healthy_and_ibd_those_with_cdr3.tsv')


cdr3_freq <- fread('../data/cdr3_freq_healthy_public_v2.tsv')
cdr3_freq_L_P <- split(cdr3_freq, by = c("length", "IMGT"))

dir_pair_path <- paste0('../HEALTHY/cdr3_hla_pairs_wo_correlation_in_alleles_using_9_PCs/')
hla_gene = c('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1','DRB1')
hla_matrices <- paste0('../hla_matrices_v3/wo_correlation/',list.files('../hla_matrices_v3/wo_correlation/'))

dir_results <- paste0('../HEALTHY/manova_results_wo_correlation_in_alleles_using_9_PCs/')

manova_df_all <- fread('../manova_results_v3.tsv')[0,]
    
if (!file.exists(dir_pair_path)) {
dir.create(dir_pair_path, recursive = TRUE)
}
if (!file.exists(dir_results)) {
dir.create(dir_results, recursive = TRUE)
}


for (L_P_group in seq(1,length(cdr3_freq_L_P))){

    cdr3_L_P <- cdr3_freq_L_P[[L_P_group]]
    L <- unique(cdr3_L_P$length)
    P <- unique(cdr3_L_P$IMGT)
    
    cdr3_matrix <- cdr3_matrix_fun(cdr3_L_P)
    cdr3_matrix$L_P <- paste0(L,':',P)
    
    ind_in_cdr3 <- unique(cdr3_L_P$patient_id)

    for (f in hla_matrices){
        
        hla_matrix <- read_tsv(f) %>% filter(patient_id %in% ind_in_cdr3)
        pca_hla <- pca_hla_fun(hla_features %>% filter(patient_id %in% ind_in_cdr3))
        
        cdr3_hla_matrix <- hla_matrix %>%
            inner_join(cdr3_matrix, by='patient_id') %>% unite('pair', G_S,L_P, sep = ':') %>%
            inner_join(pca_hla, by = 'patient_id')
        name_pair <- cdr3_hla_matrix$pair[1]
        
        write_tsv(cdr3_hla_matrix, paste0(dir_pair_path,name_pair,'_matrix.tsv'))
        manova_results <- mlm_fun(cdr3_hla_matrix, dir_results)
        manova_df_all <- rbind(manova_df_all, manova_results)
        
    }
    }


fwrite(manova_df_all, '../HEALTHY/manova_results_wo_correlation_in_alleles_using_9_PCs.tsv')


#---------------
dir_pair_path <- paste0('../HEALTHY/cdr3_hla_pairs_wo_correlation_in_alleles_using_9_PCs/')
dir_results <- paste0('../HEALTHY/manova_results_wo_correlation_in_alleles_using_9_PCs/')
manova_df_all <- fread('../manova_results_v3.tsv')[0,]

for (f in paste0(dir_pair_path, list.files(dir_pair_path))){
    manova_results <- mlm_fun(fread(f), dir_results)
    manova_df_all <- rbind(manova_df_all, manova_results)   
    }
    

fwrite(manova_df_all, '../HEALTHY/manova_results_wo_correlation_in_alleles_using_9_PCs.tsv')
#---------------



cdr3_freq <- cdr3_freq %>% filter(patient_id %in% downsampling_ibd_ids)
cdr3_freq_L_P <- split(cdr3_freq, by = c("length", "IMGT"))

dir_pair_path <- paste0('../IBD/cdr3_hla_pairs_downsampled_9PCs/')
hla_gene = c('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1','DRB1')
hla_matrices <- paste0('../hla_matrices_v3/wo_correlation/',list.files('../hla_matrices_v3/wo_correlation/'))

dir_results <- paste0('../IBD/manova_results_downsampled_9PCs/')

manova_df_all <- fread('../manova_results_v3.tsv')[0,]
    
if (!file.exists(dir_pair_path)) {
dir.create(dir_pair_path, recursive = TRUE)
}
if (!file.exists(dir_results)) {
dir.create(dir_results, recursive = TRUE)
}


for (L_P_group in seq(1,length(cdr3_freq_L_P))){

    cdr3_L_P <- cdr3_freq_L_P[[L_P_group]]
    L <- unique(cdr3_L_P$length)
    P <- unique(cdr3_L_P$IMGT)
    
    cdr3_matrix <- cdr3_matrix_fun(cdr3_L_P)
    cdr3_matrix$L_P <- paste0(L,':',P)
    
    ind_in_cdr3 <- unique(cdr3_L_P$patient_id)

    for (f in hla_matrices){
        
        hla_matrix <- read_tsv(f) %>% filter(patient_id %in% ind_in_cdr3)
        pca_hla <- pca_hla_fun(hla_features %>% filter(patient_id %in% ind_in_cdr3))
        
        cdr3_hla_matrix <- hla_matrix %>%
            inner_join(cdr3_matrix, by='patient_id') %>% unite('pair', G_S,L_P, sep = ':') %>%
            inner_join(pca_hla, by = 'patient_id')
        name_pair <- cdr3_hla_matrix$pair[1]
        
        write_tsv(cdr3_hla_matrix, paste0(dir_pair_path,name_pair,'_matrix.tsv'))
        manova_results <- mlm_fun(cdr3_hla_matrix, dir_results)
        manova_df_all <- rbind(manova_df_all, manova_results)
        
    }
    }


fwrite(manova_df_all, '../IBD/manova_results_downsampled_9PCs.tsv')