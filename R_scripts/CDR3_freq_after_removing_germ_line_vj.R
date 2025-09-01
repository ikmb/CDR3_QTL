

#path_in <- '../data/preprocessed_cdr3/'
#path_out <- '../data/cdr3_freq_with_expansion/'
#cdr3_freq_output_file <- '../data/cdr3_all_freq_with_expansion_excluded_germ_long.tsv'
#cdr3_freq_output_file_irt <- '../data/cdr3_all_irt_with_expansion_excluded_germ_long.tsv'

first_file <- TRUE
for (f in list.files(path_in)){
    cdr3 <- as.data.table(read_parquet(paste0(path_in,f)))
    cdr3 <- cdr3[(length_seq >= 12) & (length_seq <=18)]
    id <- cdr3$patient_id[1]
    cdr3_aligned <- as.data.table(t(sapply(cdr3$amino_acid, align_imgt)))
    setnames(cdr3_aligned, all_imgt_pos)
    cdr3_aligned <- cbind(cdr3,cdr3_aligned)
    rm(cdr3)
    fwrite(cdr3_aligned, paste0(path_in, '../cdr3_wide_aligned_imgt.tsv'), sep = '\t',
           quote = FALSE, row.names = FALSE, append = TRUE, 
           col.names = first_file)
    first_file <- FALSE
    
    processed_cdr3 <- melt(cdr3_aligned, id.vars = c("patient_id", "amino_acid", "v_gene", "j_gene", "length_seq", "cum_expansion_freq", "unique_count"), 
                    measure.vars = all_imgt_pos, 
                    variable.name = "IMGT", 
                    value.name = "AA")
    processed_cdr3 <- na.omit(processed_cdr3)
    rm(cdr3_aligned)
    setDT(processed_cdr3)
    
    V_germ_line[, v_gene_family := sub("-.*", "", v_gene)]
    V_germ_line[grepl('TCRBV12', v_gene)]
    
    J_germ_line[, j_gene_family := sub("-.*", "", j_gene)][
        ,length_seq := as.integer(length_seq)]
    
    # Perform an anti-join to remove rows from main_dt where v_gene and j_gene or matches v_gene_family, j_gene_family and IMGT, AA match
    processed_cdr3_removing_v_j <- processed_cdr3[!V_germ_line, 
                                            on = c("v_gene", "IMGT", "AA")][!V_germ_line, 
                                            on = c("v_gene" = "v_gene_family", "IMGT", "AA")][!J_germ_line, 
                                            on = c("j_gene", "length_seq", "IMGT", "AA")][!J_germ_line, 
                                            on = c("j_gene" = "j_gene_family", "length_seq", "IMGT", "AA")][!(AA == 'NA')]
    
    rm(processed_cdr3)
    
    cdr3_long_freq <- processed_cdr3_removing_v_j[, .(patient_id, length_seq, IMGT, AA, cum_expansion_freq, unique_count )][, .(
        cum_expansion_freq = sum(cum_expansion_freq),  # Summing cum_expansion_freq
        unique_count = sum(unique_count)               # Summing unique_count
        ), by = .(patient_id, length_seq, IMGT, AA)]
    cdr3_long_freq <- cdr3_long_freq[, `:=`(total_imgt_count = sum(unique_count),
                                            norm_freq_expand = cum_expansion_freq/sum(cum_expansion_freq),
                                            norm_freq_unique = unique_count/total_imgt_count), 
                                     by = .(patient_id, length_seq, IMGT)]
    fwrite(cdr3_long_freq, paste0(path_out, id, '_cdr3_freq_excluded_germ_long.tsv'), sep = '\t')
    }



first_file = TRUE
for (f in list.files(path_out)){
    cdr3_freq <- fread(paste0(path_out,f))
    fwrite(cdr3_freq, cdr3_freq_output_file, 
           sep = '\t', quote = FALSE, row.names = FALSE, append = TRUE, 
           col.names = first_file)
    first_file <- FALSE
}


cdr3_freq <- fread(cdr3_freq_output_file)

cdr3_freq <- cdr3_freq[(total_imgt_count > 20)]
cdr3_freq <- cdr3_freq[, `:=`(
    n_carriers = uniqueN(patient_id),
    n_ind_high_expand_freq = sum(norm_freq_expand > 0.6),
    irt_freq_unique = qnorm((rank(norm_freq_unique, na.last="keep") - 0.5) / sum(!is.na(norm_freq_unique))),
    irt_freq_expand = qnorm((rank(norm_freq_expand, na.last="keep") - 0.5) / sum(!is.na(norm_freq_expand)))
), by = .(length_seq, IMGT, AA)]

fwrite(cdr3_freq, cdr3_freq_output_file_irt, sep = '\t', quote = FALSE, row.names = FALSE)





