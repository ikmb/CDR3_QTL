
### GETTING HLA INFO---------

#--------getting allele variation part of HLA genes--------

pivot_and_filter <- function(hla_split, alleles_in_dataset){
  #transform the hla dataframe to a long format

    # pivoting df starting from the 3rd column, because the first two are 'gene' and 'allele', and -2 due to non existent status of the last two columns
  pivot_hla <- pivot_longer(hla_split, cols = 3:(length(hla_split)-2),  names_to = 'site', values_to = 'aa')
    
  hla_position_aa <- pivot_hla %>%
    filter(!aa == '*') %>%
    separate(allele, into = c('gene_group', 'protein', 'synon_substitution', 'non_coding'), sep = ':', remove = FALSE)  %>%
    select(gene, gene_group, protein, site, aa) %>%
    unite(gene_group, protein, col = 'allele', sep = ':') %>%
    filter(allele %in% alleles_in_dataset) %>%
    distinct()
  
  hla_reference <- hla_position_aa %>%
    filter(allele == first(allele)) %>% 
    rename(ref_aa = aa) %>% select(site,ref_aa)
  
  hla_position_aa <- merge(hla_position_aa,hla_reference, by = 'site') %>% 
    mutate(change = ifelse(aa == ref_aa, 0,1)) 
  
  return(hla_position_aa)
}

#--------getting the allele variations of hla genes per positions--------
patient_site_variation_fun <- function(hla_features, gene_name){
  
  hla_patients <- hla_features[gene == gene_name]
  hla_reference_long <- fread(paste0(getwd(), "/reference_data/hla_msa/",gene_name, "_long.tsv"), sep = '\t')
  hla_reference_long <- hla_reference_long[allele %in% unique(hla_patients$allele)][
      , n_variat := uniqueN(AA), by = 'site'][n_variat >= 2]
  threshold <- uniqueN(hla_patients$patient_id) * 2
  
  hla_patients <- hla_patients %>% 
      inner_join(hla_reference_long, relationship = "many-to-many", by = join_by(allele))  %>%
      group_by(patient_id, gene, site, AA) %>%
      summarise(homo_hetero = sum(homo_hetero), .groups = 'keep') %>%
      ungroup() 
  
  hla_patients_proper_sites <- hla_patients %>%
      group_by(site, AA) %>%
      summarise(n_allele_carriers = sum(homo_hetero), .groups = 'keep') %>% 
      right_join(hla_patients, by = c('site', 'AA')) %>% 
      ungroup() %>%
      group_by(site) %>%
      mutate(site_to_remove = ifelse(n_allele_carriers >= round(threshold*0.95), 'yes', 'no')) %>%
      filter(site_to_remove == 'no') %>%
      filter(n_allele_carriers >= round(threshold*0.05)) %>%
      ungroup()

  return(hla_patients_proper_sites)
}




pca_hla_fun <- function(hla_genotypes_matrix, cor_threshold, n_pcs){
    #ids_in_cdr3 <- fread('../data/ids_in_cdr3.txt', header = FALSE)
    #hla_genotypes_matrix <- hla_genotypes_matrix[patient_id %in% ids_in_cdr3$V1]
    hla_genotypes_matrix <- remove_correlated_columns(hla_genotypes_matrix, cor_threshold)  
    hla_patients_id <- hla_genotypes_matrix$patient_id
    hla_matrix <- hla_genotypes_matrix %>%
        dplyr::select_if(is.numeric) %>%
        t() %>%
        as.matrix()
    colnames(hla_matrix) <- hla_patients_id
    pca <- prcomp(hla_matrix, center = TRUE, scale. = TRUE)
    pca_hla <- as.data.frame(pca$rotation) %>% 
        dplyr::select(c(paste0('PC',seq(1,n_pcs)) )) %>% 
        mutate(patient_id = rownames(pca$rotation)) %>%
        relocate(patient_id, 1)
    rownames(pca_hla) <- NULL
    return(as.data.table(pca_hla))
    
}


#--------annotating HLA sequences with locations of antigen-binding pockets-------

hla_pocket <- function(manova_results){
    A_pocket <- c(5, 7, 59,163, 167, 171)
    B_pocket <- c(24, 34, 45, 63, 66, 67)
    C_pocket <- c(9, 70, 73, 74, 97)
    D_pocket <- c(99, 155, 159, 160)
    E_pocket <- c(114, 147, 152, 156)
    F_pocket <- c(77, 80, 81, 84, 95, 123, 143, 146)
    class_I_pockets_dict <- data.frame(Pocket = c(rep('A_pocket',length(A_pocket)), rep('B_pocket',length(B_pocket)), rep('C_pocket',length(C_pocket)), 
                                           rep('D_pocket',length(D_pocket)), rep('E_pocket',length(E_pocket)), rep('F_pocket',length(F_pocket))),
                               Site_hla = c(A_pocket,B_pocket,C_pocket,D_pocket,E_pocket,F_pocket))
    Class_I <- c('A', 'B', 'C')
    pockets_df <- tidyr::crossing(class_I_pockets_dict, HLA = Class_I)
  
    P1 <- c(7, 9, 24, 31, 32, 43, 52)
    P6 <- c(11, 62, 65, 66, 69)
    P9 <- c( 68, 69, 72, 73, 76)
    class_II_alpha_pockets_dict <- data.frame(Pocket = c(rep('P1',length(P1)), rep('P6',length(P6)), rep('P9',length(P9))), 
                             Site_hla = c(P1,P6,P9))
    Class_II_alpha <- c('DQA1', 'DPA1')
    pockets_df_II_alpha <- tidyr::crossing(class_II_alpha_pockets_dict, HLA = Class_II_alpha)
    
    P1 <- c(85, 86, 89, 90)
    P4 <- c(13, 26, 28, 70, 74, 78)
    P6 <- c(11, 30)
    P7 <- c(28, 47, 61, 67, 71)  
    P9 <- c(9, 37, 57)
    class_II_beta_pockets_dict <- data.frame(Pocket = c(rep('P1',length(P1)), rep('P4',length(P4)),rep('P6',length(P6)), 
                                             rep('P7',length(P7)), rep('P9',length(P9))), 
                                             Site_hla = c(P1,P4,P6,P7,P9))
    Class_II_beta <- c('DRB1','DQB1', 'DPB1')
    pockets_df_II_beta <- tidyr::crossing(class_II_beta_pockets_dict, HLA = Class_II_beta)
    pockets_df <- rbind(pockets_df, pockets_df_II_alpha, pockets_df_II_beta)
    
    hla_annotated <- left_join(manova_results, pockets_df, by = c("HLA", 'Site_hla'))
    return(hla_annotated)
}



hla_correlation_fun <- function(hla_matrix){

    cor_matrix <- list(cor(hla_matrix %>% select(contains('allele'))))
    names(cor_matrix) <- name_pair
    cor_list <- append(cor_list, cor_matrix)
    file_path <- "../list_of_hla_sites_correlation_v2.RData"
    save(cor_list, file = file_path)
    corrplot(cor_list[[1]], method = "circle")
    
}


remove_correlated_columns <- function(hla_matrix, cutoff_value){

    cor_matrix <- cor(hla_matrix[,-1])
    columns_to_remove <- findCorrelation(cor_matrix, cutoff = cutoff_value, exact = FALSE, names = TRUE)
    df_wo_correlations <- hla_matrix[, c(columns_to_remove) := NULL]
    
    return(df_wo_correlations)
}



