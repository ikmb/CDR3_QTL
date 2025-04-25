source('cdr3_functions.R')
source('libraries.R')

#----------with IRT---------

dir_pair_path <- paste0('../cdr3_hla_pairs/')
    
if (!file.exists(dir_pair_path)) {
dir.create(dir_pair_path, recursive = TRUE)
}

hla_matrices <- paste0('../hla_matrices/with_pcs/',list.files('../hla_matrices/with_pcs/'))
    prep_mode <- 'irt_freq_expand'

cdr3_freq <- fread('../data/cdr3_all_freq_with_IRT.tsv')
cdr3_freq <- cdr3_freq[(total_imgt_count > 20 & n_carriers > 200)]
cdr3_freq_split_l_imgt <- split(cdr3_freq, by = c('length_seq', 'IMGT'))



for (L_P_group in seq(1,length(cdr3_freq_split_l_imgt))){

    cdr3_L_P <- cdr3_freq_split_l_imgt[[L_P_group]]
    ind_in_cdr3 <- unique(cdr3_L_P$patient_id)

    L <- unique(cdr3_L_P$length_seq)
    P <- unique(cdr3_L_P$IMGT)
    
    threshold_carriers <- uniqueN(cdr3_L_P$patient_id)*0.9

    
    #-----cdr3_matrix_fun() takes the l_imgt group, filters AA with at least 'threshold' carriers, and the name of measurment for widening into matrix
    cdr3_group_long <- cdr3_prep_for_matrix_fun(cdr3_L_P, threshold_carriers)    
    
    path_group<- paste0(dir_pair_path, prep_mode, '/')
    fwrite(cdr3_group_long, paste0(path_group,L,'_',P,'.tsv'), sep = '\t')
    
    cdr3_matrix <- dcast(cdr3_group_long, patient_id ~ AA, value.var = prep_mode)
    
    cdr3_matrix$L_P <- paste0(L,'_',P)
    

    for (f in hla_matrices){
        
        hla_matrix <- fread(f)
        
        cdr3_hla_matrix <- merge(hla_matrix, cdr3_matrix, by = 'patient_id', all.x = TRUE)
        cdr3_hla_matrix[is.na(norm_freq_expand), norm_freq_expand := 0]
        final_dt[is.na(norm_freq_unique), norm_freq_unique := 0]
        cdr3_hla_matrix$pair <- paste0(hla_matrix$G_S[1],'_', cdr3_matrix$L_P[1])
        name_pair <- cdr3_hla_matrix$pair[1]

        path_out <- paste0(dir_pair_path, prep_mode, '/')
        if (!file.exists(path_out)) {
            dir.create(path_out, recursive = TRUE)
            }

        fwrite(cdr3_hla_matrix, paste0(path_out,name_pair,'_matrix.tsv'), sep = '\t')
        
    }

    
}

irt_groups <- paste0('../cdr3_hla_pairs/CDR3_GROUPS_IRT/', list.files('../cdr3_hla_pairs/CDR3_GROUPS_IRT'))
prep_mode <- 'irt_freq_unique'

for (group in irt_groups){
    cdr3_group_long <- fread(group)
    L <- unique(cdr3_group_long$length_seq)
    P <- unique(cdr3_group_long$IMGT)
    ind_in_cdr3 <- unique(cdr3_group_long$patient_id)
    cdr3_matrix <- dcast(cdr3_group_long, patient_id ~ AA, value.var = prep_mode)

    cdr3_matrix$L_P <- paste0(L,'_',P)
    

    for (f in hla_matrices){
        
        hla_matrix <- fread(f)
        
        hla_matrix <- hla_matrix[patient_id %in% ind_in_cdr3]
        
        cdr3_hla_matrix <- merge(hla_matrix, cdr3_matrix, by = 'patient_id')
        cdr3_hla_matrix$pair <- paste0(hla_matrix$G_S[1],'_', cdr3_matrix$L_P[1])
        name_pair <- cdr3_hla_matrix$pair[1]

        path_out <- paste0(dir_pair_path, prep_mode, '/')
        if (!file.exists(path_out)) {
            dir.create(path_out, recursive = TRUE)
            }

        fwrite(cdr3_hla_matrix, paste0(path_out,name_pair,'_matrix.tsv'), sep = '\t')
        
    }

    
}



#--------------with norm frequencies, log and binary state ---------
cdr3_freq <- fread('../data/cdr3_all_freq_with_IRT.tsv')
cdr3_freq <- cdr3_freq[(total_imgt_count > 20 & n_carriers > 200)]
cdr3_freq_split_l_imgt <- split(cdr3_freq, by = c('length_seq', 'IMGT'))

dir_pair_path <- paste0('../cdr3_hla_pairs/')
hla_matrices <- paste0('../hla_matrices_v3/with_pcs/',list.files('../hla_matrices_v3/with_pcs/'))

for (L_P_group in seq(1,length(cdr3_freq_split_l_imgt))){

    cdr3_L_P <- cdr3_freq_split_l_imgt[[L_P_group]]
   
    L <- unique(cdr3_L_P$length_seq)
    P <- unique(cdr3_L_P$IMGT)
        
    #-----cdr3_matrix_fun() takes the l_imgt group, filters AA with at least 'threshold' carriers, and the name of measurment for widening into matrix
    cdr3_group_long <- cdr3_prep_for_log_matrix_fun(cdr3_L_P)    
    
    path_group <- '../cdr3_hla_pairs/CDR3_GROUPS_FREQ/'
    fwrite(cdr3_group_long, paste0(path_group,L,'_',P,'.tsv'), sep = '\t')
    
    prep_modes <- c('norm_freq_expand', 'norm_freq_unique', 'log_c_freq_expand', 'log_c_freq_unique', 'binary_state' )
    for (prep_mode in prep_modes){
        cdr3_matrix <- dcast(cdr3_group_long, patient_id ~ AA, value.var = prep_mode)
            
        for (f in hla_matrices){
            
            hla_matrix <- fread(f)
            
            cdr3_hla_matrix <- merge(hla_matrix, cdr3_matrix, by = 'patient_id', all.x = TRUE)
            cdr3_hla_matrix$pair <- paste0(hla_matrix$G_S[1],'_', L, '_', P)
            name_pair <- cdr3_hla_matrix$pair[1]
    
            path_out <- paste0(dir_pair_path,prep_mode, '/')
            if (!file.exists(path_out)) {
                dir.create(path_out, recursive = TRUE)
                }
            cdr3_hla_matrix[, G_S := NULL]
            fwrite(cdr3_hla_matrix, paste0(path_out,name_pair,'_matrix.tsv'), sep = '\t')
            
        }
        }

    
}




#------with UNIQUE COUNTS------

cdr3_freq <- fread('../data/cdr3_all_freq_with_IRT.tsv')
cdr3_freq <- cdr3_freq[(total_imgt_count > 20 & n_carriers > 200)]
cdr3_freq_split_l_imgt <- split(cdr3_freq, by = c('length_seq', 'IMGT'))

dir_pair_path <- paste0('../cdr3_hla_pairs/')
hla_matrices <- paste0('../hla_matrices_v3/with_pcs/',list.files('../hla_matrices_v3/with_pcs/'))

for (L_P_group in seq(1,length(cdr3_freq_split_l_imgt))){

    cdr3_L_P <- cdr3_freq_split_l_imgt[[L_P_group]]
   
    L <- unique(cdr3_L_P$length_seq)
    P <- unique(cdr3_L_P$IMGT)
        
    #-----cdr3_matrix_fun() takes the l_imgt group, filters AA with at least 'threshold' carriers, and the name of measurment for widening into matrix
    cdr3_group_long <- cdr3_prep_for_count_matrix_fun(cdr3_L_P)    
    
    path_group <- '../cdr3_hla_pairs/CDR3_GROUPS_UNIQUE_COUNT/'
    fwrite(cdr3_group_long, paste0(path_group,L,'_',P,'.tsv'), sep = '\t')
    
    prep_modes <- c('unique_count')
    for (prep_mode in prep_modes){
        cdr3_matrix <- dcast(cdr3_group_long, patient_id ~ AA, value.var = prep_mode)
            
        for (f in hla_matrices){
            
            hla_matrix <- fread(f)
            
            cdr3_hla_matrix <- merge(hla_matrix, cdr3_matrix, by = 'patient_id', all.x = TRUE)
            cdr3_hla_matrix$pair <- paste0(hla_matrix$G_S[1],'_', L, '_', P)
            name_pair <- cdr3_hla_matrix$pair[1]
    
            path_out <- paste0(dir_pair_path,prep_mode, '/')
            if (!file.exists(path_out)) {
                dir.create(path_out, recursive = TRUE)
                }
            cdr3_hla_matrix[, G_S := NULL]
            fwrite(cdr3_hla_matrix, paste0(path_out,name_pair,'_matrix.tsv'), sep = '\t')
            
        }
        }

    
}