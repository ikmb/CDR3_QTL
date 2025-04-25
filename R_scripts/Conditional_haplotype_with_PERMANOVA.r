# inputs:
# phenotype 'both', 'I', 'H'
# hla_gene from hla_genes <- c('DRB1', 'DQB1', 'DQA1', 'B', 'DPB1', 'DPA1', 'A', 'C')
# cond_round
# using_groups TRUE, FALSE
# analyses_pairs <- 

source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')


dir_results_cond <- paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/permanova/', hla_gene, '/')
analysed_pairs <- strsplit(list.files(dir_results_cond), '_')
analysed_pairs <- unlist(lapply(analysed_pairs, function(x) paste0(x[2:3], collapse = '_')))


if (phenotype == 'both'){
    cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')
    pca_dt <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')
} else {
    cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')[patient_id %like% phenotype]
    pca_dt <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')[patient_id %like% phenotype]
}

n_ind <- uniqueN(cdr3_freq$patient_id)
imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)]
cdr3_freq$pair <- paste0(cdr3_freq$length_seq, '_', cdr3_freq$IMGT)
cdr3_freq <- cdr3_freq[!(pair %in% analysed_pairs)]

cdr3_freq_split_length <- split(cdr3_freq, cdr3_freq$length_seq)
cdr3_freq_split_length_wide <- lapply(cdr3_freq_split_length, function(x){
    l <- unique(x$length_seq)
    z <- split(x, x$IMGT)
    lapply(z, function(y){
        p <- unique(y$IMGT)
        df <- dcast(y, patient_id ~ AA, value.var = 'norm_freq_unique', fill = 0)
        as.data.table(df)[, pair := paste0(l, '_', p)]
    })
})
cdr3_lengths <- names(cdr3_freq_split_length_wide)

rm(cdr3_freq)
rm(cdr3_freq_split_length)

phenotypes <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/phenotypes.tsv')
hla_alleles_long <- fread(paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/reference_data/hla_msa/',hla_gene,'_long.tsv'))
if (phenotype == 'I'){
    hla_alleles_patients <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/ibd_hla_features.tsv')[gene == hla_gene]
} else if (phenotype == 'H'){
    hla_alleles_patients <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/healthy_hla_features.tsv')[gene == hla_gene]
} else {
    hla_alleles_patients <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/hla_features_healthy_and_ibd.tsv')[gene == hla_gene]
}
alleles_in_dataset <- unique(hla_alleles_patients$allele)

hla_var_sites <- sapply(grep(hla_gene, list.files('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/hla_matrices/as_in_Ishigaki/'), value = TRUE), function(x) {
    unlist(strsplit(x, '_'), use.names = FALSE)[[2]]
})
hla_var_sites <- unname(hla_var_sites)

if (cond_round == 0){
    for (l in cdr3_lengths){
        for (dt in cdr3_freq_split_length_wide[[l]]){

        pair <- dt$pair[1]

        cdr3_hla_matrix <- merge(dt, phenotypes, by = 'patient_id')
        amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
        aa_in_matrix <- colnames(cdr3_hla_matrix)[colnames(cdr3_hla_matrix) %in% amino_acids]
        Y_matrix <- cdr3_hla_matrix %>% dplyr::select(all_of(aa_in_matrix))
        Y_matrix <- as.data.frame(as.matrix(Y_matrix))
        #X_reduced <- "group"
        X_reduced <- sites_recategorate(hla_haplotypes_patients, site_matrix(significant_hits))
        formula_null <- as.formula(str_c('Y_matrix ~', X_reduced))
        permod0 <- adonis2(formula_null, data = cdr3_hla_matrix, parallel = (detectCores() - 2) , method="robust.aitchison", permutations = 20000)

        permanova_dt_all <- as.data.table(permod0[1,])
        permanova_dt_all$Site_hla <- c('null')

        aic_permod0 <- AICc_permanova2(permod0)
        permanova_dt_all$aic <- c(aic_permod0$AICc)
        permanova_dt_all$variance_explained <- permod0$R2[1]

        first <- TRUE
        for (hla_site in hla_var_sites){
            site_comb <- site_matrix(hla_site)
            conditional_matrix <- sites_recategorate(hla_alleles_patients, site_comb)
            #cond_covariates <- colnames(conditional_matrix)[-1]
            cond_covariates <- conditional_matrix %>% dplyr:: select(-1) %>% t() %>% as.data.frame() %>% arrange(desc(V1)) %>% slice(-1) %>% rownames()

            test <- merge(cdr3_hla_matrix, conditional_matrix, by = 'patient_id')
            Y_matrix <- test %>% dplyr::select(all_of(aa_in_matrix))
            Y_matrix <- as.data.frame(as.matrix(Y_matrix))

            X <- paste0(paste0(cond_covariates, collapse = '+'), '+group')
            #formula_full <- as.formula(str_c(paste0('Y_matrix ~', X), PCs, sep  = '+'))
            formula_full <- as.formula(str_c('Y_matrix ~', X))
            permod1 <- adonis2(formula_full, data = test, parallel = (detectCores() - 2), method="robust.aitchison", permutations = 20000)

            permanova_dt <- as.data.table(permod1[1,])
            aic_permod1 <- AICc_permanova2(permod1)
            permanova_dt$aic <- c(aic_permod1$AICc)
            permanova_dt$Site_hla <- hla_site
            permanova_dt$variance_explained <- permod1$R2[1]
            permanova_dt_all <- rbind(permanova_dt_all, permanova_dt)
            }
            
        permanova_dt_all$pair <- pair
        permanova_dt_all$HLA <- hla_gene
        permanova_dt_all <- permanova_dt_all %>% rename( Pvalue = grep("Pr", colnames(permanova_dt_all), value = TRUE)) %>% arrange(Pvalue)
        permanova_dt_all <- permanova_dt_all[, c('Length_cdr3', 'IMGT'):= tstrsplit(pair, '_', keep = c(1,2))]
        fwrite(permanova_dt_all, paste0(dir_results_cond,hla_gene,'_',pair,'_',phenotype,'.tsv'), sep = '\t')
        
    }
    significant_hits_with_length <- define_cond_hits(permanova_dt_all, 1)
    }
   
} else {
    manova_cond_all <- fread(paste0(dir_results_cond, phenotype, '/', hla_gene, '_conditional_round_',cond_round, '.tsv'))
    significant_hits_with_length <- define_cond_hits(manova_cond_all)
    min_P <- min(na.omit(manova_cond_all$Pvalue))
    rm(manova_cond_all)

    while (min_P < 0.05){
        first <- TRUE

        for (i in seq_along(significant_hits_with_length)){
            l <- names(significant_hits_with_length[i])
            significant_hits <- significant_hits_with_length[[l]]
            hla_reduced <- sites_recategorate(hla_alleles_patients, site_matrix(significant_hits))
            hla_reduced_with_pca <- merge(hla_reduced, pca_dt, by = 'patient_id')
            x_reduced_alleles <- colnames(hla_reduced)[-1]
            path_cond_manova <- paste0(dir_results_cond, phenotype, '/', hla_gene, '_conditional_round_',length(significant_hits), '.tsv')
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
                    if (using_groups == TRUE){
                        test <- merge(test, phenotypes, by = 'patient_id')
                    }   
                    tryCatch({
                        manova_cond <- conditional_fun(test, x_reduced_alleles, cond_covariates)
                        manova_cond$condition <- cond
                        setnames(manova_cond, old = grep('Pr', colnames(manova_cond), value = TRUE), new = 'Pvalue')
                        manova_cond <- manova_cond[, c('Length_cdr3', 'IMGT'):= tstrsplit(pair, '_', keep = c(1,2))]
                        if (first){
                            fwrite(manova_cond, path_cond_manova, sep = '\t', append = FALSE)
                            first <- FALSE

                        } else {
                            fwrite(manova_cond, path_cond_manova, sep = '\t', append = TRUE)
                        }
                        
                    }, error = function(e){
                        print(e)
                    })
                    
                }
            }
            }
        manova_cond_all <- fread(path_cond_manova, sep = '\t')
        significant_hits_with_length <- define_cond_hits(manova_cond_all)
        min_P <- min(na.omit(manova_cond_all$Pvalue))
        rm(manova_cond_all)
    }
}