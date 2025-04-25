
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/hla_functions.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')

mmlm_phenotype <- function(cdr3_hla_matrix){
    name_pair <- unique(cdr3_hla_matrix$pair)

    amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
    Y_matrix <- as.matrix(cdr3_hla_matrix %>% dplyr::select(any_of(amino_acids)))
    Y <- paste0('cbind(', paste(colnames(Y_matrix), collapse = ','),") ~")
    X_full <- 'group'
    X_null <- '1'

    formula_full <- as.formula(str_c(Y, X_full))
    formula_null <- as.formula(str_c(Y, X_null))
    mod1 <- lm(formula_full, data = cdr3_hla_matrix)
    mod0 <- lm(formula_null, data = cdr3_hla_matrix)
    manova_results <- anova(mod0, mod1)

    mod1_mvlm <- mvlm(formula_full, cdr3_hla_matrix)
    var_exp_full <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
    p_val_full <- mod1_mvlm$pv["Omnibus Effect",1]

    mvlm_df <- data.table(P_val = p_val_full, variance_explained = var_exp_full)
    if (exists('mvlm_df')){
        manova_df <- cbind(data.table(manova_results), mvlm_df)
    }else{
        manova_df <- data.table(manova_results)
    }
    manova_df$pair <- name_pair
        return(manova_df)
    
}

dir_results <- '/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/'

phenotypes <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/phenotypes.tsv')
cdr3_freq <- fread('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3_all_freq_with_IRT.tsv')
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

rm(cdr3_freq)
rm(cdr3_freq_split_length)


first <- TRUE
for (i in seq_along(cdr3_freq_split_length_wide)){
    l <- names(cdr3_freq_split_length_wide)[i]
    for (dt in cdr3_freq_split_length_wide[[l]]){
        test <- merge(dt, phenotypes, by = 'patient_id')
        pair <- test$pair[1]
        mlmm_groups <- mmlm_phenotype(test)
        setnames(mlmm_groups, old = grep('Pr', colnames(mlmm_groups), value = TRUE), new = 'Pvalue')
        mlmm_groups <- mlmm_groups[, c('Length_cdr3', 'IMGT'):= tstrsplit(pair, '_', keep = c(1,2))]
        if (first){
            mlmm_groups_all <- mlmm_groups
            first <- FALSE

        } else {
            mlmm_groups_all <- rbind(mlmm_groups, mlmm_groups_all)
        }
    }     
}