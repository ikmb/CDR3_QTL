phenotypes <- fread('../phenotypes.tsv')
#hla_matrices <- '../hla_matrices/with_pcs/'
#path_in <- 

cdr3_freq <- fread(paste0(path_in,'/cdr3_all_irt_excluded_germ_long.tsv')) 

cdr3_freq_split <- split(cdr3_freq, by = c("length_seq", "IMGT"))
pcs <- 3
prep_mode <- c('irt_freq_unique')
hla_matrices_files <- paste0(hla_matrices, list.files(hla_matrices))

first <- TRUE
for (cdr3_L_P in cdr3_freq_split){
    if (uniqueN(cdr3_L_P$patient_id) < 200){
        next
    }
   
    L <- unique(cdr3_L_P$length_seq)
    P <- unique(cdr3_L_P$IMGT)
        
    cdr3_matrix <- dcast(cdr3_L_P, patient_id ~ AA, value.var = prep_mode, fill = 0)
    
    dir_results <- paste0(path_in, '/manova_results/', prep_mode, '/')
    if (!file.exists(dir_results)) {
            dir.create(dir_results, recursive = TRUE)
            }
    
    for (f in hla_matrices_files){
        
        hla_matrix <- fread(f)
        
        pair_matrix <- merge(hla_matrix, cdr3_matrix, by = 'patient_id', all.y = TRUE)
        pair_matrix <- merge(pair_matrix, phenotypes, by = 'patient_id', all.x = TRUE)
        pair_matrix <- na.omit(pair_matrix)
        pair_matrix$pair <- paste0(hla_matrix$G_S[1],'_', L, '_', P)
        name_pair <- pair_matrix$pair[1]

        pair_matrix[, G_S := NULL]
        path_pairs <- paste0(dir_results, '../cdr3_hla_pairs/', prep_mode, '/')
        if (!file.exists(path_pairs)) {
            dir.create(path_pairs, recursive = TRUE)
            }
        name_pair <- pair_matrix$pair[1]
        
        fwrite(pair_matrix, paste0(path_pairs, name_pair,'.tsv'), sep = '\t')
        
        
        manova_df <- mlm_fun(pair_matrix, dir_results, n_pcs = pcs) %>% 
            separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = '_', remove = FALSE)
        if (first){
            manova_df_all <- manova_df[0,]
            first <- FALSE
        }
        if (is.null(manova_df)){
            next
        } else {
            manova_df_all <- rbind(manova_df_all, manova_df) 
            }
        }
    
        
    manova_df_all <- manova_df_all %>% 
        separate(pair, into = c('HLA', 'Site_hla', 'Length_cdr3', 'Position_cdr3'), sep = '_', remove = FALSE)
    fwrite(manova_df_all, paste0(dir_results, L, '_',P, '_main_manova.tsv'), sep = '\t')
                    
        }


