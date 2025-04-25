source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')

amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
X_reduced <- "group"

cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')
n_ind <- uniqueN(cdr3_freq$patient_id)
imgt_to_discard <- c('P104', 'P105', 'P106', 'P117', 'P118')
cdr3_freq <- cdr3_freq[n_carriers >= 300][!(IMGT %in% imgt_to_discard)]
cdr3_freq$length_seq <- paste0('L', cdr3_freq$length_seq)
cdr3_freq$pair <- paste0(cdr3_freq$length_seq, '_', cdr3_freq$IMGT)


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
man_whit_results <- data.table(matrix(nrow = 0, ncol = 5) )
colnames(man_whit_results) <- c('datasets', 'pair', 'AA', 'p_val', 'cliff_delta')

for (phen in list(c('CD','HLH'))){
    first <- TRUE
    man_whit_results <- man_whit_results[0,]
    phen_dt <- phenotypes %>% filter(group %in% phen)
    dataset_combination <- paste0(phen[1], '_', phen[2])
    phen_dt$datasets <- dataset_combination
    for (l in cdr3_lengths){
        for (dt in cdr3_freq_split_length_wide[[l]]){
            pair <- dt$pair[1]
            cdr3_hla_matrix <- merge(dt, phen_dt, by = 'patient_id')
            aa_in_matrix <- colnames(cdr3_hla_matrix)[colnames(cdr3_hla_matrix) %in% amino_acids]
            Y_matrix <- cdr3_hla_matrix %>% dplyr::select(all_of(aa_in_matrix))
            Y_matrix <- as.data.frame(as.matrix(Y_matrix))
            
            formula_null <- as.formula(str_c('Y_matrix ~', X_reduced))
            permod0 <- adonis2(formula_null, data = cdr3_hla_matrix, parallel = (detectCores() - 2) , method="robust.aitchison", permutations = 9999)

            permanova_dt <- as.data.table(permod0[1,])
            permanova_dt$datasets <- dataset_combination

            aic_permod0 <- AICc_permanova2(permod0)
            permanova_dt$aic <- c(aic_permod0$AICc)
            permanova_dt$variance_explained <- permod0$R2[1]
            permanova_dt$pair <- pair

            if (first){
                permanova_dt_all <- permanova_dt
                first <- FALSE  
            } else {
                permanova_dt_all <- rbind(permanova_dt_all, permanova_dt)
            } 
            if (permanova_dt[[grep('Pr',  names(permanova_dt), value = TRUE)]] < 0.001) {
                            for (aa in aa_in_matrix) {
                                
                                man_whit <- wilcox.test(cdr3_hla_matrix[[aa]] ~ group, data = cdr3_hla_matrix)
                                if (man_whit$p.value < 0.05){
                                    cliff_delta <- cliff.delta(d = cdr3_hla_matrix[[aa]],
                                            f = cdr3_hla_matrix$group)
                                    man_whit_results_row <- c(dataset_combination, pair, aa, man_whit$p.value, cliff_delta$estimate )
                                    names(man_whit_results_row) <- c('datasets', 'pair', 'AA', 'p_val', 'cliff_delta')

                                    man_whit_results_row <- data.table(t(man_whit_results_row))
                                    man_whit_results <- rbind(man_whit_results, man_whit_results_row)
                                }
                            }}

    }
    }
    man_whit_results <- man_whit_results[,c("Length_cdr3", "IMGT") := tstrsplit(pair, "_", keep = c(1,2))]
    permanova_dt_all <- permanova_dt_all %>% rename( Pvalue = grep("Pr", colnames(permanova_dt_all), value = TRUE)) %>% arrange(Pvalue)
    permanova_dt_all <- permanova_dt_all[, c('Length_cdr3', 'IMGT'):= tstrsplit(pair, '_', keep = c(1,2))]
    fwrite(permanova_dt_all, paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/cdr3_compositions_permanova/',dataset_combination,'.tsv'), sep = '\t')
    fwrite(man_whit_results, paste0('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/cdr3_compositions_permanova/man_whit_results_for_',dataset_combination,'.tsv'), 
        sep = '\t')
}
