#pca_hla <- fread('../data/hla/pca_hla_9PCs.tsv')
#hla_features_all <- hla_features[0,]

hla_genes = c('A','C','B','DRB1','DQA1', 'DQB1', 'DPA1','DPB1')

for (gene_name in hla_genes){   
 
    hla_patients_proper_sites <- patient_site_variation_fun(hla_features, gene_name)
    
    site_split <- hla_patients_proper_sites %>% dplyr::select(patient_id, homo_hetero, site, AA) %>% 
        group_split(site)
    
    for (site in seq(1, length(site_split))){
        site_location <- site_split[[site]]$site[1]
        site_matrix <- site_split[[site]] %>% dplyr::select(-site) %>% 
            pivot_wider(names_from=AA, values_from = homo_hetero, values_fn = sum, values_fill = 0, names_prefix = 'allele_') %>%
            inner_join(pca_hla, by = 'patient_id') 
        
        site_matrix$G_S <- paste0(gene_name,'_',site_location)
        write_tsv(site_matrix, paste0("../hla_matrices_v3/with_pcs/",gene_name, "_",site_location, ".tsv"))
    }
    

}
                                                         
                                                         