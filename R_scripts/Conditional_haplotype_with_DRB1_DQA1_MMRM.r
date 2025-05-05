# inputs:
# phenotype 'both', 'I', 'H'
# hla_gene from hla_genes <- c('DQB1', 'DQA1', 'B', 'DPB1', 'DPA1', 'A', 'C')
# cond_round
# using_groups TRUE, FALSE
# path_manova


source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')


dir_results_cond <- paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/with_haplotypes/', phenotype, '/')
if (!file.exists(dir_results_cond)) {
    dir.create(dir_results_cond, recursive = TRUE)
    }
phenotypes <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/phenotypes.tsv')

dir_results_cond_permanova <- paste0(dir_results_cond, 'permanova/')
if (!file.exists(dir_results_cond_permanova)) {
    dir.create(dir_results_cond_permanova, recursive = TRUE)
    }

analysed_pairs <- strsplit(list.files(dir_results_cond_permanova), '_')
analysed_pairs <- unlist(lapply(analysed_pairs, function(x) paste0(x[2:3], collapse = '_')))

if (phenotype == 'both'){
    cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')
} else if (phenotype!= 'H' & phenotype != 'I'){
    cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')[group %like% phenotype]
} else {
    cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')[patient_id %like% phenotype]
}

ids <- unique(cdr3_freq$patient_id)
n_ind <- uniqueN(cdr3_freq$patient_id)

pca_dt <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/pca_hla_healthy_and_ibd.tsv')[patient_id %in% ids]

imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= (n_ind/2)][!(IMGT %in% imgt_to_discard)][length_seq %in% c('L13', 'L14', 'L15', 'L16', 'L17')]

#############
cdr3_freq$pair <- paste0(cdr3_freq$length_seq, '_', cdr3_freq$IMGT)
##############

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

rm(cdr3_freq)
rm(cdr3_freq_split_length)

all_hla <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/all_hla_variations_long.tsv')
hla_site_correlations <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/hla_site_correlations.tsv')

hla_alleles_long <- fread(paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/reference_data/hla_msa/DRB1_long.tsv'))
hla_alleles_long$site <- paste0('DRB1_', hla_alleles_long$site)
hla_alleles_patients <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/hla_features_healthy_and_ibd.tsv')
### !!!!!!!!
hla_alleles_patients <- hla_alleles_patients[gene == 'DRB1'][patient_id %in% ids] ### !!!! hla_gene == 'DRB1' because we condition on DRB1
### !!!!!!!!!
alleles_in_dataset <- unique(hla_alleles_patients$allele)

patients_allele <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/hla/hla_genotypes_matrix.csv')
hla_alleles_correlations_matrix <- cor(patients_allele[,-1], method = 'spearman')
hla_alleles_correlations <- melt(hla_alleles_correlations_matrix)
hla_alleles_correlations <- as.data.table(hla_alleles_correlations)

#all_hla_with_haplotypes <- merge(haplotypes, all_hla, by = 'allele')
#hla_haplotype_reference <- unique(all_hla_with_haplotypes[,.(haplotype,allele, Site_hla, AA)])
#hla_haplotypes_patients <- unique(unique(all_hla_with_haplotypes[, .(haplotype, patient_id, homo_hetero)])[,
    #homo_hetero := min(homo_hetero), by = c('haplotype', 'patient_id')])

hla_var_sites <- sapply(grep(paste0(hla_gene,'_'), list.files('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/hla_matrices/as_in_Ishigaki/'), value = TRUE), function(x) {
    gsub('_matrix.tsv','', x)
})
hla_var_sites <- unname(hla_var_sites)


#bonf <- 0.05 / nrow(na.omit(manova_cond_all))
#min_P <- min(na.omit(manova_cond_all$Pvalue))
#rm(manova_cond_all)

conditional_sites <- c('DRB1_13', 'DRB1_37', 'DRB1_71')
if (cond_round == 1){
    significant_hits <- conditional_sites
} else {
    manova_cond_all <- fread(path_manova)
    significant_hits_with_length <- define_cond_hits(manova_cond_all)
}

hla_reduced <- sites_recategorate(hla_alleles_patients, site_matrix(significant_hits))
hla_reduced_with_pca <- merge(hla_reduced, pca_dt, by = 'patient_id')
x_reduced_alleles <- colnames(hla_reduced)[-1]

path_cond_manova <- paste0(dir_results_cond, hla_gene, '_with_drb1_13_37_71_round_',cond_round, '.tsv')
hla_var_sites_to_test <- hla_var_sites[!(hla_var_sites %in% significant_hits)]
combinations_to_test <- lapply(hla_var_sites_to_test, function(x) c(significant_hits,x))

amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
first <- TRUE
for (l in cdr3_lengths){
    for (dx in cdr3_freq_split_length_wide[[l]]){
    
        pair <- dx$pair[1]

        cdr3_hla_matrix <- merge(dx, phenotypes, by = 'patient_id')
        cdr3_hla_matrix <- merge(cdr3_hla_matrix, hla_reduced, by = 'patient_id')


        for (j in seq_along(combinations_to_test)){
            sites <- combinations_to_test[[j]]
            dt <- all_hla[Site_hla %in% sites][
                , Site_hla_AA := paste0(Site_hla, '_', AA)][, .(patient_id, allele, homo_hetero, Site_hla_AA)]
            
            dt2 <- dt[grepl(paste(sites, collapse = "|"), Site_hla_AA)]
            
            dt_dqa1 <- dt2[grepl(paste0(hla_gene, '_'), Site_hla_AA)]
            dt_drb1 <- dt2[grepl("^DRB1", Site_hla_AA)]
            
            # Merge by patient (cross join)
            merged <- merge(
                dt_drb1, dt_dqa1,
                by = "patient_id", allow.cartesian = TRUE, suffixes = c("_drb1", "_dqa1")
            )
            merged <- merged[, .(
                Site_hla_AA_drb1 = paste(unique(Site_hla_AA_drb1), collapse = "_"),
                homo_hetero_drb1 = unique(homo_hetero_drb1),
                allele_dqa1 = unique(allele_dqa1),
                homo_hetero_dqa1 = unique(homo_hetero_dqa1),
                Site_hla_AA_dqa1 = unique(Site_hla_AA_dqa1)
                # Assuming col3 is always the same within group
            ), by = .(patient_id, allele_drb1)]
            # Create haplotype label
            merged[, haplo := paste0(Site_hla_AA_drb1, "_", Site_hla_AA_dqa1)]
            
            #dqa1_haplotypes <- hla_site_correlations[HLA1 != HLA2][Var1 %like% paste0(sites[[2]],'_')][Var2 %like% paste0(sites[[1]],'_')][value > -0.1]
            dqa1_haplotypes_alleles <- hla_alleles_correlations[value > 0.1][Var1!=Var2][Var1 %like% 'DRB1'][Var2 %like% hla_gene]
            
            valid <- merged[dqa1_haplotypes_alleles, 
                on = .(allele_dqa1 = Var2, allele_drb1 = Var1), 
                nomatch = 0
            ]
            valid_cleaned <- unique(valid[, .(patient_id, homo_hetero_drb1, homo_hetero_dqa1, haplo)])[, 
                homo_hetero := min(homo_hetero_drb1, homo_hetero_dqa1), by = patient_id][,
                homo_hetero := if (.N == 1) 2 else homo_hetero, by = patient_id] %>% arrange(patient_id)
            conditional_matrix <- dcast(valid_cleaned, patient_id ~ haplo, value.var = 'homo_hetero', fill = 0, fun.aggregate = sum)
            cond <- paste(sites, collapse = '_')
            cond_covariates <- colnames(conditional_matrix)[-1]

            test <- merge(cdr3_hla_matrix, conditional_matrix, by = 'patient_id')
            if (!using_groups){
                test <- test %>% dplyr::select(-c('group'))
            }
            
            tryCatch({
                manova_cond <- conditional_fun(test, x_reduced_alleles, cond_covariates, using_pcs = FALSE)
                manova_cond$condition <- cond
                setnames(manova_cond, old = grep('Pr', colnames(manova_cond), value = TRUE), new = 'Pvalue')
                manova_cond <- manova_cond[, c('Length_cdr3', 'IMGT'):= tstrsplit(pair, '_', keep = c(1,2))]
                if (first){
                    fwrite(manova_cond, path_cond_manova, sep = '\t')
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
sig_sites <- na.omit(manova_cond_all) %>% 
    group_by(Length_cdr3) %>% 
    filter(Pvalue == min(Pvalue)) %>% 
    filter(variance_explained == max(variance_explained)) %>% unique()
significant_hits_with_length <- lapply(sig_sites$condition, function(x) unlist(strsplit(x, ':')))

names(significant_hits_with_length) <- sig_sites$Length_cdr3
min_P <- min(na.omit(manova_cond_all$Pvalue))
#rm(manova_cond_all)
