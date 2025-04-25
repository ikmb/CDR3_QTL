source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')


cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/cdr3_all_irt_with_expansion_excluded_germ_long.tsv')
n_ind <- uniqueN(cdr3_freq$patient_id)
imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)]
cdr3_freq$length_seq <- paste0('L', cdr3_freq$length_seq)

cdr3_freq_split_length <- split(cdr3_freq, cdr3_freq$length_seq)
cdr3_freq_split_length_wide <- lapply(cdr3_freq_split_length, function(x){
    l <- unique(x$length_seq)
    z <- split(x, x$IMGT)
    lapply(z, function(y){
        p <- unique(y$IMGT)
        df <- dcast(y, patient_id ~ AA, value.var = 'irt_freq_unique', fill = 0)
        as.data.table(df)[, pair := paste0(l, '_', p)]
    })
})
cdr3_lengths <- names(cdr3_freq_split_length_wide)


hla_gene <- 'DQA1'
cond_round <- 3

hla_alleles_long <- fread(paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/reference_data/hla_msa/',hla_gene,'_long.tsv'))
hla_alleles_patients <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/hla/hla_features_matched_DeWitt.tsv')[gene == hla_gene]
alleles_in_dataset <- unique(hla_alleles_patients$allele)

pca_dt <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/pca_hla_v2.tsv')[,
    patient_id_Emerson := paste0('P', patient_id_Emerson)]
setnames(pca_dt, new = 'patient_id', old = 'patient_id_Emerson')

hla_var_sites <- sapply(grep(paste0(hla_gene,'_'), list.files('/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/hla_matrices/'), value = TRUE), function(x) {
    unlist(strsplit(x, '_'), use.names = FALSE)[[2]]
})
hla_var_sites <- unname(hla_var_sites)

dir_results_cond <- '/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/conditional/'

if (cond_round == 0){
    manova_cond_all <- na.omit(fread('/work_beegfs/sukmb667/projects/cdr3-qtl/Emerson_DeWitt/manova_results/manova_results_3PCs.tsv'))[, 
        Length_cdr3 := paste0('L', Length_cdr3)]
    setnames(manova_cond_all, old = grep('Pr', colnames(manova_cond_all), value = TRUE), new = 'Pvalue')
    significant_hits_with_length <- define_cond_hits(manova_cond_all[HLA == hla_gene], 1)
} else {
    manova_cond_all <- fread(paste0(dir_results_cond, hla_gene, '_conditional_round_',cond_round, '.tsv'))
    significant_hits_with_length <- define_cond_hits(manova_cond_all)
}

bonf <- 0.05 / nrow(na.omit(manova_cond_all))


while (min(na.omit(manova_cond_all$Pvalue)) <= bonf){
    
    first <- TRUE
    for (i in seq_along(significant_hits_with_length)){
        l <- names(significant_hits_with_length[i])
        significant_hits <- significant_hits_with_length[[l]]
        hla_reduced <- sites_recategorate(hla_alleles_patients, site_matrix(significant_hits))
        hla_reduced_with_pca <- merge(hla_reduced, pca_dt, by = 'patient_id')
        x_reduced_alleles <- colnames(hla_reduced)[-1]

        hla_var_sites_to_test <- hla_var_sites[!(hla_var_sites %in% significant_hits)]
        combinations_to_test <- lapply(hla_var_sites_to_test, function(x) c(significant_hits,x))
        for (j in seq_along(combinations_to_test)){
            sites <- combinations_to_test[[j]]
            cond <- paste(sites, collapse = '_')
            site_comb <- site_matrix(sites)
            conditional_matrix <- sites_recategorate(hla_alleles_patients, site_comb)
            cond_covariates <- colnames(conditional_matrix)[-1]
            for (dt in cdr3_freq_split_length_wide[[l]]){
                cdr3_matrix <- merge(dt, hla_reduced_with_pca, by = 'patient_id')
                test <- merge(cdr3_matrix, conditional_matrix, by = 'patient_id')
                try({
                    manova_cond <- conditional_fun(test, x_reduced_alleles, cond_covariates)
                    manova_cond$condition <- cond
                    if (first){
                        manova_cond_all <- manova_cond
                        first <- FALSE
                    } else {
                        manova_cond_all <- rbind(manova_cond_all, manova_cond)
                    }
                }, silent = TRUE)
                
            }
        }
        }
    setnames(manova_cond_all, old = grep('Pr', colnames(manova_cond_all), value = TRUE), new = 'Pvalue')
    manova_cond_all <- manova_cond_all %>% separate(pair, into = c('Length_cdr3', 'IMGT'), sep = '_')
    significant_hits_with_length <- define_cond_hits(manova_cond_all)    
    fwrite(manova_cond_all, paste0(dir_results_cond, hla_gene, '_conditional_round_',length(significant_hits), '.tsv'), sep = '\t')

}