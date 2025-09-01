
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

dir_manova <- '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/'
if (using_groups == TRUE){
    dir_results_manova <- paste0(dir_manova, phenotype, '/with_groups/')
    phenotypes <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/phenotypes.tsv')
} else {
    dir_results_manova <- paste0(dir_manova, phenotype, '/')
    }


cdr3_file_path <- '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv'
cdr3_freq <- cdr3_freq_by_phenotype(cdr3_file_path, phenotype)

ids <- unique(cdr3_freq$patient_id)
if (downsampling){
    ids <- sample(ids, downsampling_size)
    cdr3_freq <- cdr3_freq[patient_id %in% ids]
    dir_results_manova <- paste0(dir_results_manova, 'downsampled/')
} 

if (!dir.exists(dir_results_manova)){
    dir.create(dir_results_manova, recursive = TRUE)
}

n_ind <- uniqueN(cdr3_freq$patient_id)

pca_dt <- fread('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')[patient_id %in% ids]

imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)]

cdr3_freq_split_length <- split(cdr3_freq, cdr3_freq$Length_cdr3)
cdr3_freq_split_length_wide <- lapply(cdr3_freq_split_length, function(x){
    l <- unique(x$Length_cdr3)
    z <- split(x, x$IMGT)
    lapply(z, function(y){
        p <- unique(y$IMGT)
        df <- dcast(y, patient_id ~ AA, value.var = 'irt_freq_unique', fill = 0)
        as.data.table(df)[, pair := paste0(l, '_', p)]
    })
})
cdr3_lengths <- names(cdr3_freq_split_length_wide)

rm(cdr3_freq)
rm(cdr3_freq_split_length)

for (hla_gene in hla_genes){
    hla_paths <- grep(paste0(hla_gene, "_"), list.files('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/hla_matrices/with_pcs/',
        pattern = '.tsv', full.names = TRUE), value = TRUE)
    path_manova <- paste0(dir_results_manova, 'manova_results_', phenotype,'_',hla_gene, '.tsv')

    first <- TRUE

    for (l in cdr3_lengths){
            for (dt in cdr3_freq_split_length_wide[[l]]){
                for (hla in hla_paths){
                    hla_with_pca <- fread(hla)[patient_id %in% ids]
                    hla_site <- unique(hla_with_pca$G_S)
                    cdr3_matrix <- merge(dt, hla_with_pca, by = 'patient_id')
                    cdr3_matrix$pair <- paste0(cdr3_matrix$pair, "_", hla_site)
                    if (using_groups == TRUE){
                        cdr3_matrix <- merge(cdr3_matrix, phenotypes, by = 'patient_id') 
                    }
                tryCatch({

                manova_df <- mlm_fun(cdr3_matrix, n_pcs = pcs)
                #setnames(manova_df, old = grep('Pr', colnames(manova_df), value = TRUE), new = 'Pvalue')
                #manova_df <- manova_df[, c('Length_cdr3', 'IMGT', 'HLA', 'Site_hla'):= tstrsplit(pair, '_', keep = c(1,2,3,4))]

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
}

