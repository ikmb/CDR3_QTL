hla_variation_Emerson <- fread('../../Emerson_DeWitt//HLA_DeWitt_variations.tsv') %>% filter(site > 0, n_variat>1) %>% dplyr::select(gene, site, allele, aa)
hla_patients_Emerson <- fread('../../Emerson_DeWitt//hla_features_DeWitt.tsv') %>% mutate(patient_id = paste0('E', patient_id)) %>% dplyr::select(-gene)

alleles_summary <- hla_patients_Emerson %>% group_by(allele) %>% 
    summarise(sum_alleles = sum(homo_hetero), .groups = 'keep') 
alleles_summary_with_sites <- hla_variation_Emerson %>%
    left_join(alleles_summary, by = "allele") %>%
    dplyr::select(-allele) %>%
    group_by(gene, site, aa) %>% 
    summarise(sum_alleles_per_aa = sum(sum_alleles), .groups = 'keep') %>% ungroup()
alleles_n_variat <- alleles_summary_with_sites %>% group_by(gene, site) %>% 
    mutate(site_to_remove = ifelse(any(sum_alleles_per_aa > 1193), 'yes', 'no'), n_variat = ifelse(any(site_to_remove == 'yes'), 1, uniqueN(aa))) %>% ungroup()


