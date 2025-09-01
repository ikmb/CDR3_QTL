
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')

# inputs:
# phenotype 'both', 'I', 'H', 'UC', 'CD'
# downsampling TRUE, FALSE
# downsampling_size 1000
# using_groups TRUE, FALSE
# pcs <- 3
# hla_gene from hla_genes <- c('DQB1', 'DQA1', 'B', 'DPB1', 'DPA1', 'A', 'C')

dir_manova <- '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/only_expanded/'
if (using_groups == TRUE){
    dir_results_manova <- paste0(dir_manova, phenotype, '/with_groups/')
    phenotypes <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/phenotypes.tsv')
} else {
    dir_results_manova <- paste0(dir_manova, phenotype, '/')
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
    dir_results_manova <- paste0(dir_results_manova, 'downsampled/')
} else {
    ids <- unique(cdr3_freq$patient_id)
}

if (!dir.exists(dir_results_manova)){
    dir.create(dir_results_manova, recursive = TRUE)
}

n_ind <- uniqueN(cdr3_freq$patient_id)

pca_dt <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')[patient_id %in% ids]

imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)]

# ----- For fully expanded I don't split by length, but by IMGT only -----
#cdr3_freq_split_length <- split(cdr3_freq, cdr3_freq$length_seq)
#cdr3_freq_split_length_wide <- lapply(cdr3_freq_split_length, function(x){
#----------- this part is changed----
cdr3_freq_split_imgt <- split(cdr3_freq, cdr3_freq$IMGT)
cdr3_freq_split_imgt_wide <- lapply(cdr3_freq_split_imgt, function(x){
    p <- unique(x$IMGT)
    df <- dcast(x, patient_id ~ AA, value.var = freq_measure, fill = 0)
        as.data.table(df)[, IMGT := p]
    })


rm(cdr3_freq)
rm(cdr3_freq_split_imgt)

for (hla_gene in hla_genes){
    hla_paths <- grep(paste0(hla_gene, "_"), list.files('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/hla_matrices/with_pcs/',
        pattern = '.tsv', full.names = TRUE), value = TRUE)
    path_manova <- paste0(dir_results_manova, 'expanded_manova_results_', phenotype,'_',hla_gene, '.tsv')

    first <- TRUE

    for (dt in cdr3_freq_split_imgt_wide){
        for (hla in hla_paths){
            hla_with_pca <- fread(hla)[patient_id %in% ids]
            hla_site <- unique(hla_with_pca$G_S)
            cdr3_matrix <- merge(dt, hla_with_pca, by = 'patient_id')
            cdr3_matrix$pair <- paste0(cdr3_matrix$IMGT, "_", hla_site)
            if (using_groups == TRUE){
                cdr3_matrix <- merge(cdr3_matrix, phenotypes, by = 'patient_id') 
            }
            tryCatch({

            manova_df <- mlm_fun(cdr3_matrix, n_pcs = pcs)
            #setnames(manova_df, old = grep('Pr', colnames(manova_df), value = TRUE), new = 'Pvalue')
            manova_df <- manova_df[, c('Length_cdr3', 'IMGT', 'HLA', 'Site_hla'):= tstrsplit(pair, '_', keep = c(1,2,3))]

            if (first & !file.exists(path_manova)){
                fwrite(manova_df, path_manova, sep = '\t')
                first <- FALSE

            } else {
                fwrite(manova_df, path_manova, sep = '\t', append = TRUE)
            }

            }, error = function(e){
                print(e)
            })
            
    }
}
        
}

