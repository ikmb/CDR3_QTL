# inputs:
# phenotype 'both', 'I', 'H', 'UC', 'CD'
# hla_gene from hla_genes <- c('DQB1', 'DQA1', 'B', 'DPB1', 'DPA1', 'A', 'C')
# cond_round
# using_groups TRUE, FALSE
# path_manova
# downsampling TRUE, FALSE
# downsampling_size

source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')

define_cond_hits_wo_length <- function(manova_results, condit_round = 2, gene_names_in_cond = FALSE){
    #bonf_cor <- 0.05/nrow(manova_results)/2  # deviding by 2 because we have two models for each pair 
    if (condit_round == 1){

        if (!('Pvalue' %in% colnames(manova_results))){
            setnames(manova_results, old = grep('Pr', colnames(manova_results), value = TRUE), new = 'Pvalue')
        }
        sig_sites <- na.omit(manova_results) %>% 
            filter(variance_explained > 0) %>%
            filter(Pvalue == min(Pvalue)) %>% 
            filter(variance_explained == max(variance_explained)) %>% unique()
        sig_sites <- sig_sites$Site_hla
    } else {
        sig_sites <- na.omit(manova_results) %>% 
            filter(variance_explained > 0) %>%
            filter(Pvalue == min(Pvalue)) %>% 
            filter(variance_explained == max(variance_explained)) %>% unique()
        sig_sites <- sig_sites$Site_hla
        if (gene_names_in_cond | grepl(paste0(hla_genes, collapse = '|'),sig_sites$condition[1])){
            sig_sites_with_length <- lapply(sig_sites$condition, function(x) {
            condition <- strsplit(x, '_')[[1]]
            paired <- mapply(function(a, b) paste(a, b, sep = "_"),
                            condition[seq(1, length(condition), by = 2)],
                            condition[seq(2, length(condition), by = 2)])
            unname(paired)
        })
        } else {
            sig_sites_with_length <- lapply(sig_sites$condition, function(x) unlist(strsplit(x, '_')))
        } 

    }
    return(sig_sites)
}



dir_cond <- '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/only_expanded/'
if (using_groups == TRUE){
    dir_results_cond <- paste0(dir_cond, phenotype,'/with_groups/')
    phenotypes <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/phenotypes.tsv')
} else {
    dir_results_cond <- paste0(dir_cond, phenotype, '/')
    }

cdr3_freq <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_irt_only_expanded_excluded_germ_long.tsv')
if (phenotype == 'both'){
    cdr3_freq <- cdr3_freq
    freq_measure <- 'irt_freq_expand'
} else if (phenotype!= 'H' & phenotype != 'I'){
    cdr3_freq <- cdr3_freq[group %like% phenotype]
    freq_measure <- 'irt_freq_expand'
} else {
    cdr3_freq <- cdr3_freq[patient_id %like% phenotype]
    freq_measure <- 'irt_freq_expand_group'
}

if (downsampling){
    ids <- sample(unique(cdr3_freq$patient_id), downsampling_size)
    cdr3_freq <- cdr3_freq[patient_id %in% ids]
    dir_results_cond <- paste0(dir_results_cond, 'downsampled/')
} else {
    ids <- unique(cdr3_freq$patient_id)
}

if (!dir.exists(dir_results_cond)){
    dir.create(dir_results_cond, recursive = TRUE)
}

n_ind <- uniqueN(cdr3_freq$patient_id)

pca_dt <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')[patient_id %in% ids]

imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)]

cdr3_freq_split_imgt <- split(cdr3_freq, cdr3_freq$IMGT)
cdr3_freq_split_imgt_wide <- lapply(cdr3_freq_split_imgt, function(x){
    p <- unique(x$IMGT)
    df <- dcast(x, patient_id ~ AA, value.var = freq_measure, fill = 0)
        as.data.table(df)[, IMGT := p]
    })


rm(cdr3_freq)
rm(cdr3_freq_split_imgt)


hla_alleles_long <- fread(paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/reference_data/hla_msa/',hla_gene,'_long.tsv'))
hla_alleles_long$site <- paste0(hla_gene, '_', hla_alleles_long$site)

hla_alleles_patients <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/hla_features_healthy_and_ibd.tsv')
hla_alleles_patients <- hla_alleles_patients[gene == hla_gene][patient_id %in% ids]

alleles_in_dataset <- unique(hla_alleles_patients$allele)

hla_var_sites <- sapply(grep(paste0(hla_gene,'_'), list.files('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/hla_matrices/as_in_Ishigaki/'), value = TRUE), function(x) {
    gsub('_matrix.tsv', '', x)
})
hla_var_sites <- unname(hla_var_sites)

if (cond_round == 0){
    manova_cond_all <- na.omit(fread(path_manova))[HLA == hla_gene]

    if (!('Pvalue' %in% colnames(manova_cond_all))){
        setnames(manova_cond_all, old = grep('Pr', colnames(manova_cond_all), value = TRUE), new = 'Pvalue')
    }
    significant_hits <- define_cond_hits_wo_length(manova_cond_all, 1)
    significant_hits <- lapply(significant_hits, function(x) paste(hla_gene, x, sep = '_'))
} else {
    manova_cond_all <- fread(paste0(dir_results_cond, hla_gene, '_conditional_round_',cond_round, '.tsv'))
    significant_hits <- define_cond_hits_wo_length(manova_cond_all)
}
bonf <- 0.01 / uniqueN(manova_cond_all$pair)
min_P <- min(na.omit(manova_cond_all$Pvalue))
max_var_expl <- max(manova_cond_all$variance_explained)
rm(manova_cond_all)

while (min_P < bonf & max_var_expl > 0){
    first <- TRUE
    cond_round <- cond_round + 1
   
    hla_reduced <- sites_recategorate(hla_alleles_patients, site_matrix(significant_hits))
    hla_reduced_with_pca <- merge(hla_reduced, pca_dt, by = 'patient_id')
    x_reduced_alleles <- colnames(hla_reduced)[-1]
    path_cond_manova <- paste0(dir_results_cond, hla_gene, '_conditional_round_',cond_round, '.tsv')
    hla_var_sites_to_test <- hla_var_sites[!(hla_var_sites %in% significant_hits)]
    combinations_to_test <- lapply(hla_var_sites_to_test, function(x) c(significant_hits,x))

    for (j in seq_along(combinations_to_test)){
        sites <- combinations_to_test[[j]]
        cond <- paste(sites, collapse = '_')
        site_comb <- site_matrix(sites)
        conditional_matrix <- sites_recategorate(hla_alleles_patients, site_comb)
        cond_covariates <- colnames(conditional_matrix)[-1]

        for (dt in cdr3_freq_split_imgt_wide){
            cdr3_matrix <- merge(dt, hla_reduced_with_pca, by = 'patient_id')
            test <- merge(cdr3_matrix, conditional_matrix, by = 'patient_id')
            if (using_groups == TRUE){
                test <- merge(test, phenotypes, by = 'patient_id')
            }   
            tryCatch({
                manova_cond <- conditional_fun(test, x_reduced_alleles, cond_covariates)
                manova_cond$condition <- cond
                setnames(manova_cond, old = grep('Pr', colnames(manova_cond), value = TRUE), new = 'Pvalue')
                if (first & !file.exists(path_cond_manova)){
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
    manova_cond_all <- fread(path_cond_manova, sep = '\t')
    significant_hits <- define_cond_hits_wo_length(manova_cond_all)
    min_P <- min(na.omit(manova_cond_all$Pvalue))
    max_var_expl <- max(manova_cond_all$variance_explained)
    rm(manova_cond_all)
}