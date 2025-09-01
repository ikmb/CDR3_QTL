
n_rounds <- 4


first <- TRUE
for (phenotype in c('I', 'H', 'Both')){
    if (phenotype == 'Both'){
        path_to_phenotype <- '/with_groups/both/'
    } else {
        path_to_phenotype <- paste0(phenotype,'/')
    }
    for (hla_gene in hla_genes){

        n_tests <- 0
        for (i in seq(1,n_rounds)){
            file_path <- paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/',path_to_phenotype,hla_gene,'_conditional_round_',i,'.tsv')
            n_df <- nrow(fread(file_path))
            n_tests <- n_tests + n_df/2
        }
        for (i in seq(1,n_rounds)){
            file_path <- paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/',path_to_phenotype,hla_gene,'_conditional_round_',i,'.tsv')
            df_rounds <- na.omit(fread(file_path)) %>% 
                filter(Pvalue <= 0.05/n_tests) %>%
                group_by(Length_cdr3) %>%
                filter(Pvalue == min(Pvalue)) %>%
                filter(variance_explained == max(variance_explained)) %>%
                rowwise() %>%
                mutate(Site_hla = as.numeric(unlist(strsplit(condition, '_'))[[i+1]]), Cond_round = i,
                    Dataset = ifelse(phenotype == 'I', 'IBD', ifelse(phenotype == 'H', 'Healthy', 'Combined')), 
                    HLA = hla_gene, Pvalue = ifelse(Pvalue == 0, smallest_number, Pvalue)) %>%
                ungroup() %>%
                dplyr::select(Dataset, HLA, Length_cdr3, IMGT, Pvalue, Omnibus, condition, Site_hla, Cond_round) %>% 
                rename('variance_explained' = 'Omnibus') %>% 
                unique()
            if (first){
                df_rounds_all <- df_rounds
                first <- FALSE
            } else {
                df_rounds_all <- rbind(df_rounds_all, df_rounds)
            }
        }
        
}}

df_combined <- na.omit(fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/main_manova_as_in_Ishigaki.tsv')[, 
    Length_cdr3 := paste0('L', Length_cdr3)][HLA %in% c('DQA1', 'DQB1', 'DRB1')])  %>% 
    rename('Pvalue' = grep('Pr', names(.), value = TRUE), 'IMGT' = 'Position_cdr3') %>% 
    group_by(HLA, Length_cdr3) %>%  
    filter(Pvalue == min(Pvalue)) %>% 
    filter(variance_explained == max(variance_explained))%>%
    mutate(Dataset = 'Combined', Cond_round = 0, condition = unique(Site_hla)) %>%
    ungroup() %>%
    dplyr::select(HLA, Dataset, Length_cdr3, IMGT, Pvalue, variance_explained, Site_hla, Cond_round, condition) %>%
    unique()


df_healthy <- na.omit(fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/HEALTHY/manova_results_wo_correlation_in_PCA.tsv'))
df_healthy <- df_healthy[,
    c('HLA', 'Site_hla', 'Length_cdr3', 'IMGT'):= tstrsplit(pair,':')][, 
    Length_cdr3 := paste0('L', Length_cdr3)][HLA %in% c('DQA1', 'DQB1', 'DRB1')] %>% 
    rename('Pvalue' = grep('Pr', names(df_healthy), value = TRUE)) %>% 
    group_by(HLA, Length_cdr3) %>% unique() %>% 
    filter(Pvalue == min(Pvalue)) %>% 
    filter(variance_explained == max(variance_explained))%>%
    mutate(Dataset = 'Healthy', Cond_round = 0, condition = unique(Site_hla)) %>%
    ungroup() %>%
    dplyr::select(HLA, Dataset, Length_cdr3, IMGT, Pvalue, variance_explained, Site_hla, Cond_round, condition) %>% 
     unique()

df_ibd <- na.omit(fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/IBD/manova_results_downsampled_wo_correlation_in_alleles_using_3_PCs.tsv'))
df_ibd <- df_ibd[,
    c('HLA', 'Site_hla', 'Length_cdr3', 'IMGT'):= tstrsplit(pair,':')][, 
    Length_cdr3 := paste0('L', Length_cdr3)][HLA %in% c('DQA1', 'DQB1', 'DRB1')] %>% 
    rename('Pvalue' = grep('Pr', names(df_ibd), value = TRUE)) %>% 
    group_by(HLA, Length_cdr3) %>% unique() %>% 
    filter(Pvalue == min(Pvalue)) %>% 
    filter(variance_explained == max(variance_explained)) %>%
    mutate(Dataset = 'IBD', Cond_round = 0, condition = unique(Site_hla)) %>%
    ungroup() %>%
    dplyr::select(HLA, Dataset, Length_cdr3, IMGT, Pvalue, variance_explained, Site_hla, Cond_round, condition) %>% 
    unique()

df_all_combined <- rbind(df_rounds_all, df_combined, df_healthy, df_ibd)

fwrite(, df_all_combined, '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/conditional_cdr3qtls_results.tsv', sep = '\t')

cdr3_qtl_results <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/conditional_cdr3qtls_results.tsv')

cdr3_qtl_results <- df_all_combined %>% 
    group_by(Dataset, HLA, Site_hla) %>% 
    mutate(n_occur = n(), Site_hla = as.numeric(Site_hla), Pvalue = ifelse(Pvalue == 0, smallest_number, Pvalue)) %>% 
    ungroup() %>% 
    filter(n_occur > 1) 