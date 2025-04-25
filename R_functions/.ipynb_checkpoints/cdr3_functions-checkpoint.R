
cdr3_dt_preproccesing <- function(path_to_parquet_file, sample_id, path_out){
    
    cdr3_dt <- cdr3_dt[sequenceStatus == 'In',                             
        .(amino_acid = aminoAcid,
          v_gene = vGeneName, 
          j_gene = jGeneName, 
          v_family = vFamilyName, 
          j_family = jFamilyName, 
          expansion_freq = frequencyNormalized)][
         , patient_id := sample_id][
         , unique_count := 1][
         , length_seq := nchar(amino_acid)][  
         , expansion_freq := as.numeric(expansion_freq) / 100][
         , c("v_gene", "v_gene_2") := tstrsplit(v_gene, "/")][  
         , v_gene := fifelse(is.na(v_gene), v_family, v_gene)][  
         , j_gene := fifelse(is.na(j_gene), j_family, j_gene)][  
         , .(unique_count = sum(unique_count), cum_expansion_freq = sum(expansion_freq)),                             
         by = .(patient_id, amino_acid, v_gene, j_gene, length_seq)][
         (length_seq >= 12) & (length_seq <=18)] 
    write_parquet(cdr3_dt, paste0(path_out,'preprocessed_',sample_id,'.parquet'))
    }
    


#-------all IMGT positions--------
all_imgt_pos <- c(
  "P104","P105","P106","P107","P108","P109","P110","P111",
  "P111.1","P112.2","P112.1",
  "P112","P113","P114","P115","P116","P117","P118")


#------alignment in matrix----
align_imgt <- function(CDR3){
    AA <- unlist(strsplit(CDR3, ""));
    N_AA_fow <- nchar(CDR3) %/% 2 + nchar(CDR3) %% 2;
    N_AA_na <- 18 - nchar(CDR3);
    AAmod <- c(
        AA[1:(N_AA_fow)],
        rep("NA",N_AA_na),
        AA[(N_AA_fow + 1):nchar(CDR3)]
    );
    names(AAmod) <- all_imgt_pos
    return(t(AAmod))
    #------another version----
    #AAmod_matrix <- matrix(AAmod, ncol=18)
    #colnames(AAmod_matrix) <- all_imgt_pos
    #return(AAmod_matrix)
}
#------------------------------


#-------------V and J germ-line database-----

#AAs encoded in V gene germline sequences
V_germ_line <- as.data.table(read.table("../../reference_data/references/V_exclude_AA.info",
                          header=T,colClasses="character"))
V_germ_line$IMGT <- paste0("P",V_germ_line$IMGT)
setnames(V_germ_line, "Vgene", "v_gene")
V_exclude <- paste0(V_germ_line$Vgene,":",V_germ_line$IMGT,":",V_germ_line$AA)

#AAs encoded in J gene germline sequences (we need to consider CDR3 length for J genes)
J_germ_line <- as.data.table(read.table("../../reference_data/references/J_exclude_AA.info",
                          header=T,colClasses="character"))
J_germ_line$IMGT <- paste0("P",J_germ_line$IMGT)
setnames(J_germ_line, "Jgene", "j_gene")
setnames(J_germ_line, "Length", "length_seq")
J_exclude <- paste0(J_germ_line$Jgene,':',J_germ_line$Length,":",J_germ_line$IMGT,":",J_germ_line$AA)




###--------updated version-------24.06.2024-------
v_j_exclude <- function(cdr3_group, path_out){
    
    cdr3_aligned <- as.data.frame(t(sapply(cdr3_group$amino_acids, align_imgt)))
    colnames(cdr3_aligned) <- all_imgt_pos
    cdr3_group_aligned <- cdr3_group %>% bind_cols(cdr3_aligned)
    
    v <- cdr3_group$v_gene[1]
    j <- cdr3_group$j_gene[1]
    l <- cdr3_group$length_seq[1]
    
    V_exclude <- V_germ_line %>% filter(Vgene == v) %>% dplyr::select(-Vgene)
    if (nrow(V_exclude) == 0){
        V_exclude <- V_germ_line[grep(v,V_germ_line$Vgene),] %>% dplyr::select(IMGT, AA)
    }
    
    J_exclude <- J_germ_line %>% filter(Jgene == j & Length == l) %>% dplyr::select(-c(Jgene, Length))
    if (nrow(J_exclude) == 0){
        J_exclude <- J_germ_line[grep(j, J_germ_line$Jgene),] %>% filter(Length == l) %>% 
            dplyr::select(IMGT, AA)
    }
    
    VJ <- full_join(V_exclude,J_exclude)
    VJ_list <- setNames(VJ$AA, VJ$IMGT)
    
    cdr3_group_aligned_excluded_germ <- cdr3_group_aligned %>%
      mutate(across(names(VJ_list), ~if_else(. == VJ_list[cur_column()], NA_character_, .)))
    
    return(cdr3_group_aligned)

    }



v_j_exclude_old <- function(cdr3_group){
    
    cdr3_aligned <- t(sapply(cdr3_group$amino_acid, function(x){align_imgt(x)}))
    colnames(cdr3_aligned) <- all_imgt_pos
    cdr3_group_aligned <- cdr3_group %>% bind_cols(cdr3_aligned)
    v <- cdr3_group$v_gene[1]
    j <- cdr3_group$j_gene[1]
    l <- cdr3_group$length[1]
    
    V_exclude <- V_germ_line_aligned %>% filter(Vgene == v) %>% select(-Vgene)
    J_exclude <- J_germ_line_aligned %>% filter(Jgene == j & Length == l) %>% select(-c(Jgene, Length))
    V_exclude <- V_exclude[,colSums(is.na(V_exclude)) == 0]
    J_exclude <- J_exclude[,colSums(is.na(J_exclude)) == 0]
  
    cdr3_group_aligned_excluded_germ <- cdr3_group_aligned
    
    positionsV <- colnames(V_exclude)
    for (pos in positionsV){
        matching <- which(cdr3_group_aligned_excluded_germ[pos] == unlist(V_exclude[pos]))
        cdr3_group_aligned_excluded_germ[matching,pos] <- 'NA'                                                             
    }
    positionsJ <- colnames(J_exclude)
    for (pos in positionsJ){
        matching <- which(cdr3_group_aligned_excluded_germ[pos] == unlist(J_exclude[pos]))
        cdr3_group_aligned_excluded_germ[matching,pos] <- 'NA'                                                             
    }
    
    cdr3_group_aligned_excluded_germ_long <- cdr3_group_aligned_excluded_germ %>% select(patient_id, count, length, contains('P')) %>% 
        pivot_longer(names_to = 'IMGT', values_to = 'AA', cols = 4:21) %>% 
        filter(AA != 'NA') %>%
        group_by(patient_id,IMGT, AA) %>%
        mutate(count = sum(count)) %>% 
        distinct() %>%
        ungroup()
    
    write.table(cdr3_group_aligned_excluded_germ_long, '../cdr3/cdr3_excluded_germ_long.tsv', sep = '\t', quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
    return(cdr3_group_aligned_excluded_germ_long)
    }


cdr3_freq_with_clonal_expansion <- function(cdr3_aligned){
    cdr3_long <- cdr3_aligned %>% select(patient_id, count, length, contains('P')) %>% 
        pivot_longer(names_to = 'IMGT', values_to = 'AA', cols = 4:21) %>% 
        filter(AA != 'NA') %>%
        group_by(patient_id,length, IMGT, AA) %>%
        summarise(clonal_count = sum(count), .groups = 'keep') %>% 
        ungroup() %>% 
        group_by(patient_id, length, IMGT) %>%
        mutate(total_IMGT_count = sum(clonal_count)) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(relativ_freq = total_AA_count/total_IMGT_count) %>%
        ungroup() %>% 
        group_by(length, IMGT, AA) %>%
        mutate(n_carriers = length(unique(patient_id)), n_ind_high_freq = sum(relativ_freq > 0.6)) %>%
        ungroup()
    
    write_tsv(cdr3_long, '../data/cdr3_freq_healthy_and_ibd_clonal_expansion.tsv')
    return(cdr3_long)
}


cdr3_long_fun <- function(cdr3_aligned){
    dt <- cdr3_aligned[, c("v_gene", "j_gene", "count") := NULL]
    processed_dt <- melt(dt, id.vars = c("patient_id", "amino_acids","length"), 
                        measure.vars = all_imgt_pos, 
                        variable.name = "IMGT", 
                        value.name = "AA")
    processed_dt <- na.omit(processed_dt)
    write_tsv(processed_dt, '../data/cdr3_freq_healthy_and_ibd_long.tsv', append = TRUE, col_names = FALSE)
}


cdr3_freq_fun <- function(cdr3_long){
    

    dt_summary <- cdr3_long[, .(total_AA_count = .N), by = .(patient_id, length, IMGT, AA)]
    
    # Filter the rows where total_AA_count is greater than 5
    dt_summary <- dt_summary[total_AA_count > 5]
    
    dt_imgt <- dt_summary[, .(total_IMGT_count = sum(total_AA_count)), by = .(patient_id, length, IMGT)]
    dt_summary <- merge(dt_summary, dt_imgt, by = c("patient_id","length","IMGT"), all = TRUE)                            
    
    # Filter the rows where total_IMGT_count is greater than 100
    dt_filtered <- dt_summary[total_IMGT_count > 100]
    
    dt_freq <- dt_filtered[,relativ_freq := total_AA_count/total_IMGT_count]
    
    cdr3_frequencies <- dt_freq[, `:=`(
      n_carriers = uniqueN(patient_id),
      n_ind_high_freq = sum(relativ_freq > 0.6)
    ), by = .(length, IMGT, AA)]
    
    fwrite(cdr3_frequencies, '../data/cdr3_freq_healthy_and_ibd_updated_expansion.tsv', sep = '\t')
    return(cdr3_frequencies)
}


cdr3_prep_for_irt_matrix_fun <- function(dt){ 
    
    dt <- dt[, .(patient_id, length_seq, IMGT, AA, norm_freq_expand, norm_freq_unique)]
    dt <- dt[, `:=`(
        irt_freq_unique = qnorm((rank(norm_freq_unique, na.last="keep") - 0.5) / sum(!is.na(norm_freq_unique))),
        irt_freq_expand = qnorm((rank(norm_freq_expand, na.last="keep") - 0.5) / sum(!is.na(norm_freq_expand)))
    ), by = .(length_seq, IMGT, AA)]
    
    # Step 1: Define all possible amino acids (AA)
    all_AAs <- unique(dt$AA)
    
    # Step 2: Use CJ to cross-join patient_id, IMGT, length_seq, and all_AAs to expand the data
    expanded_dt <- CJ(patient_id = dt$patient_id,
       length_seq = dt$length_seq,
        IMGT = dt$IMGT,
        AA = all_AAs, unique = TRUE)
       
    # Step 3: Merge with original data to fill in the frequencies
    # This will merge on 'patient_id', 'length_seq', 'IMGT', and 'AA'
    final_dt <- merge(expanded_dt, dt, by = c("patient_id", "length_seq", "IMGT", "AA"), all.x = TRUE)
    
    # Step 4: Replace NA frequencies with 0 or with min(freq) 9.40581138056186e-09
    final_dt[is.na(norm_freq_expand), norm_freq_expand := 0]
    final_dt[is.na(norm_freq_unique), norm_freq_unique := 0]
    final_dt[is.na(irt_freq_unique), irt_freq_unique := 0]
    final_dt[is.na(irt_freq_expand), irt_freq_expand := 0]
    final_dt[, log_c_freq_expand := log(norm_freq_expand + 0.5)]
    final_dt[, log_c_freq_unique := log(norm_freq_unique + 0.5)]
    final_dt <- final_dt[, `:=`(
        irt_postproc_freq_unique = qnorm((rank(norm_freq_unique, na.last="keep") - 0.5) / sum(!is.na(norm_freq_unique))),
        irt_postproc_freq_expand = qnorm((rank(norm_freq_expand, na.last="keep") - 0.5) / sum(!is.na(norm_freq_expand)))
    ), by = .(length_seq, IMGT, AA)]
    



    return(final_dt)
}

cdr3_prep_for_log_matrix_fun <- function(dt){ 
    
    dt <- dt[, .(patient_id, length_seq, IMGT, AA, n_carriers, norm_freq_expand, norm_freq_unique)]
    
    # Step 1: Define all possible amino acids (AA)
    all_AAs <- unique(dt$AA)

    ids <- fread('../data/ids_in_cdr3.txt', header = FALSE)
    
    # Step 2: Use CJ to cross-join patient_id, IMGT, length_seq, and all_AAs to expand the data
    expanded_dt <- CJ(patient_id = ids$V1,
       length_seq = dt$length_seq,
        IMGT = dt$IMGT,
        AA = all_AAs, unique = TRUE)
    
    # Step 3: Merge with original data to fill in the frequencies
    # This will merge on 'patient_id', 'length_seq', 'IMGT', and 'AA'
    final_dt <- merge(expanded_dt, dt, by = c("patient_id", "length_seq", "IMGT", "AA"), all.x = TRUE)
    final_dt[, n_carriers := nafill(n_carriers, type = "nocb"), by = AA]
    
    # Step 4: Replace NA frequencies with 0 or with min(freq) 9.40581138056186e-09
    final_dt[is.na(norm_freq_expand), norm_freq_expand := 0]
    final_dt[is.na(norm_freq_unique), norm_freq_unique := 0]

    final_dt[, log_c_freq_expand := log(norm_freq_expand + 0.5)]
    final_dt[, log_c_freq_unique := log(norm_freq_unique + 0.5)]
    final_dt[, binary_state := ifelse(norm_freq_expand > 0, 1, 0)]

    
    return(final_dt)
}

cdr3_prep_for_count_matrix_fun <- function(dt){ 
    
    dt <- dt[, .(patient_id, length_seq, IMGT, AA, unique_count)]
    
    # Step 1: Define all possible amino acids (AA)
    all_AAs <- unique(dt$AA)

    ids <- fread('../data/ids_in_cdr3.txt', header = FALSE)
    
    # Step 2: Use CJ to cross-join patient_id, IMGT, length_seq, and all_AAs to expand the data
    expanded_dt <- CJ(patient_id = ids$V1,
       length_seq = dt$length_seq,
        IMGT = dt$IMGT,
        AA = all_AAs, unique = TRUE)
    
    # Step 3: Merge with original data to fill in the frequencies
    # This will merge on 'patient_id', 'length_seq', 'IMGT', and 'AA'
    final_dt <- merge(expanded_dt, dt, by = c("patient_id", "length_seq", "IMGT", "AA"), all.x = TRUE)
     
    final_dt[is.na(unique_count), unique_count := 0]

  
    return(final_dt)
}



downsampling_based_on_expansion <- function(cdr3_data, repertoire_size){
    #total_frequency <- sum(cdr3_data$unique_count)
    #cdr3_data$prob <- cdr3_data$unique_count / total_frequency
    sampled_clonotypes <- data.table(sample(
      x = cdr3_data$tags,
      size = repertoire_size,
      replace = TRUE, # Sampling with replacement to reflect relative probabilities
      prob = cdr3_data$cum_expansion_freq))
    cdr3_downsampled_dt <- sampled_clonotypes[, c("amino_acid", "v_gene", "j_gene", "length_seq") := tstrsplit(V1, split = ":", fixed = TRUE)][, V1 := NULL]
    return(cdr3_downsampled_dt)
    }

downsampling_based_on_unique_freq <- function(cdr3_data, repertoire_size){
    total_frequency <- sum(cdr3_data$unique_count)
    cdr3_data$prob <- cdr3_data$unique_count / total_frequency
    sampled_clonotypes <- data.table(sample(
      x = cdr3_data$tags,
      size = repertoire_size,
      replace = TRUE, # Sampling with replacement to reflect relative probabilities
      prob = cdr3_data$prob))
    cdr3_downsampled_dt <- sampled_clonotypes[, c("amino_acid", "v_gene", "j_gene", "length_seq") := tstrsplit(V1, split = ":", fixed = TRUE)][, V1 := NULL]
    return(cdr3_downsampled_dt)
    }



