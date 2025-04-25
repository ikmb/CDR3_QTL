dfs_together <- function(dir_path, path_out){
    manova_files <- list.files(path = dir_path, pattern = "*.tsv", full.names = TRUE)
    manova_df_all <- rbindlist(lapply(manova_files, fread))
    fwrite(manova_df_all, path_out, sep = '\t')    
}

define_variables <- function(pair_matrix){
    # method can be 'manova' or 'permanova' or 'conditional'
    amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
    Y_matrix <- as.matrix(pair_matrix %>% dplyr::select(any_of(amino_acids)))
    
    without_ref_allele <- pair_matrix %>% 
        summarise(across(contains('allele'), sum)) %>% 
        t() %>% 
        as.data.frame() %>% 
        arrange(desc(V1)) %>% 
        slice(-1) %>% 
        rownames()
    
    tryCatch({
        if (length(without_ref_allele)>=2){
            hla_alleles <- colnames(remove_correlated_columns(pair_matrix %>% 
                                                              dplyr::select(all_of(without_ref_allele)), 0.9))
        } 
        }, 
        error = function(e) {
        message("Error: ", e)
        },
        finally = {
            if (!exists('hla_alleles')){
                hla_alleles <- colnames(pair_matrix %>% dplyr::select(all_of(without_ref_allele)))}
            
        })
    

    return(list(hla_alleles, Y_matrix))
}


define_formulas <- function(hla_alleles, Y_matrix, n_pcs = 3, sites = NULL, method = 'manova'){
    
    Y <- paste0('cbind(', paste(colnames(Y_matrix), collapse = ','),") ~")
    hlas <- paste(hla_alleles, collapse = '+')
    pcs <- paste0('PC', seq(1,n_pcs), collapse = '+')

    if (method == 'manova'){
        X_full <- paste(hlas, pcs, 'group', sep = '+')
        X_null <- paste(pcs, 'group', sep = '+')

        formula_full <- as.formula(str_c(Y,X_full))
        formula_null <- as.formula(str_c(Y,X_null))
    } else if (method == 'permanova'){
        X_full <- paste(hlas, 'group', sep = '+')
        X_null <- paste('group')

        formula_full <- as.formula(str_c('Y_matrix ~', X_full))
        formula_null <- as.formula(str_c('Y_matrix ~', X_null))

    } else if (method == 'conditional'){
        comb_sites <- paste(sites, collapse = '+')
        X_full <- paste(comb_sites, pcs, 'group', sep = '+')
        X_null <- paste(hlas, pcs, 'group', sep = '+')
  
        formula_full <- as.formula(str_c(Y, X_full))
        formula_null <- as.formula(str_c(Y, X_null))

    }
        
    return(list(formula_full, formula_null))
    
}



mlm_fun <- function(cdr3_hla_matrix, writing_mode = 'no', dir_results = "work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/irt_freq_unique/", n_pcs = 3){
    
    bonf <- 0.0001/nrow(cdr3_hla_matrix)
    variables <- define_variables(cdr3_hla_matrix)
    hla_alleles <- variables[[1]]
    Y_matrix <- variables[[2]]

    formulas <- define_formulas(hla_alleles, Y_matrix, n_pcs = 3, method = 'manova')
        
    mod1 <- lm(formulas[[1]], data = cdr3_hla_matrix)
    mod0 <- lm(formulas[[2]], data = cdr3_hla_matrix)
    name_pair <- unique(cdr3_hla_matrix$pair)
    manova_results <- anova(mod1, mod0)
    manova_df <- as.data.table(manova_results)
    manova_df[, pair := name_pair]
    
    tryCatch({
    
    if (sum(is.na(rowSums(coef(mod1)))) > 0){
        
        allele_to_exclude <- rownames(coef(mod1))[is.na(rowSums(coef(mod1)))]
        new_hla_alleles <- setdiff(variables$hla_alleles, allele_to_exclude)
        formulas <- define_formulas(new_hla_alleles, Y_matrix, n_pcs = 3, method = 'manova')
    
    } 
    
    }, 
    error = function(e) {
        # Handle the error
        message("Error: ", e)
    },
    finally = {
        mod1_mvlm <- mvlm(formulas[[1]], cdr3_hla_matrix)
        mod0_mvlm <- mvlm(formulas[[2]], cdr3_hla_matrix)
        
        var_exp_full <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
        p_val_full <- mod1_mvlm$pv["Omnibus Effect",1]
        
        var_exp_null <- mod0_mvlm$pseudo.rsq["Omnibus Effect",1]
        p_val_null <- mod0_mvlm$pv["Omnibus Effect",1]
        
        variance_explained <- var_exp_full - var_exp_null
        names(variance_explained) <- 'var_expl'
        P_val <-  mapply(c, p_val_null, p_val_full)
        colnames(P_val) <- 'P_val'
        Omnibus = mapply(c, var_exp_null, var_exp_full)
        colnames(Omnibus) <- 'Omnibus'
        
    
        mvlm_df <- data.table(Omnibus, model = c('null', 'full'),P_val, variance_explained)
        
    })
        
       
    if (writing_mode == 'yes'){

        mod1_summary <- summary(mod1)
        mod0_summary <- summary(mod0)

        dir_results_summary <- paste0(dir_results,'/mmlm_summary/')
        if (!file.exists(dir_results_summary)) {
            dir.create(dir_results_summary, recursive = TRUE)
        }

        mod1_mvlm_summary <- summary(mod1_mvlm)
        mod0_mvlm_summary <- summary(mod0_mvlm)
        
        dir_results_summary_mvlm <- paste0(dir_results,'/mvlm_summary/')
        if (!file.exists(dir_results_summary_mvlm)) {
            dir.create(dir_results_summary_mvlm, recursive = TRUE)
            }

        if (na.omit(manova_df)$`Pr(>F)` < bonf){
        
    
        lapply(summary(mod1), function(x){
            capture.output(x, file=paste0(dir_results_summary,name_pair,'_mod1.txt'), append = TRUE)})
        lapply(summary(mod0), function(x){
            capture.output(x, file=paste0(dir_results_summary,name_pair,'_mod0.txt'), append = TRUE)})
        lapply(summary(mod1_mvlm), function(x){
            capture.output(x, file=paste0(dir_results_summary_mvlm,name_pair,'_mod1.txt'), append = TRUE)})
        lapply(summary(mod0_mvlm), function(x){
            capture.output(x, file=paste0(dir_results_summary_mvlm,name_pair,'_mod0.txt'), append = TRUE)}) 
    } 
}
    manova_df <- cbind(manova_df, mvlm_df)
    return(manova_df) 

    }                    


permanova_fun <- function(cdr3_hla_matrix, cond_covariates, x_reduced_alleles, n_permut = 999, n_cor = 23, writing_mode = 'no'){
    
    #variables <- define_variables(cdr3_hla_matrix)
    #hla_alleles <- variables[[1]]
    #Y_matrix <- variables[[2]]

    #formulas <- define_formulas(hla_alleles, Y_matrix,  method = 'permanova')
    n_pcs <- 3
    amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
    aa_in_matrix <- colnames(cdr3_hla_matrix)[colnames(cdr3_hla_matrix) %in% amino_acids]
    Y_matrix <- cdr3_hla_matrix %>% dplyr::select(all_of(aa_in_matrix))
    Y_matrix <- as.data.frame(as.matrix(Y_matrix))
    X <- paste0(cond_covariates,  collapse = '+')
    X_reduced <- paste0(x_reduced_alleles,  collapse = '+')
    PCs <- paste0('PC', seq(1,n_pcs), collapse = '+')

    formula_full <- as.formula(str_c(paste0('Y_matrix ~', X), sep   = '+'))
    formula_null <- as.formula(str_c('Y_matrix ~', X_reduced))

    permod1 <- adonis2(formula_full, data = cdr3_hla_matrix, parallel = n_cor, method="robust.aitchison", permutations = n_permut)
    permod0 <- adonis2(formula_null, data = cdr3_hla_matrix, parallel = n_cor, method="robust.aitchison", permutations = n_permut)
    
    pair_name <- cdr3_hla_matrix$pair[1]
    permanova_dt <- rbind(as.data.table(permod1[1,]), as.data.table(permod0[1,]))
    permanova_dt$model <- c('full', 'null')
    permanova_dt$pair <- pair_name
    permanova_dt$var_exp <- permod1$R2[1] - permod0$R2[1]
    aic_permod1 <- AICc_permanova2(permod1)
    aic_permod0 <- AICc_permanova2(permod0)
    permanova_dt$aic <- c(aic_permod1$AICc, aic_permod0$AICc)
    
    if (writing_mode == 'yes'){
        path_bd <- '/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/permanova_results/PERMANOVA_9999perm.db'
        conn <- dbConnect(RSQLite::SQLite(), dbname = path_bd)
        dbWriteTable(conn, path_bd, permanova_dt, append = TRUE)
        dbDisconnect(conn)}
        
    return(permanova_dt)
}


### CONDITIONAL PART

define_cond_hits <- function(manova_results, condit_round = 2){
    #bonf_cor <- 0.05/nrow(manova_results)/2  # deviding by 2 because we have two models for each pair 
    if (condit_round == 1){
        sig_sites <- na.omit(manova_results) %>% 
            group_by(Length_cdr3) %>% 
            filter(Pvalue == min(Pvalue)) %>% 
            filter(variance_explained == max(variance_explained)) %>% unique()
        sig_sites_with_length <- sig_sites$Site_hla
    } else {
        sig_sites <- na.omit(manova_results) %>% 
            group_by(Length_cdr3) %>% 
            filter(Pvalue == min(Pvalue)) %>% 
            filter(variance_explained == max(variance_explained)) %>% unique()
        sig_sites_with_length <- lapply(sig_sites$condition, function(x) unlist(strsplit(x, '_')))
    }
    names(sig_sites_with_length) <- sig_sites$Length_cdr3
    return(sig_sites_with_length)
}

site_matrix <- function(sites){
    ref_sites <- data.table::dcast(hla_alleles_long[site %in% sites], allele ~ site, value.var = 'AA')[allele %in% unique(hla_alleles_patients$allele)]
    combin <- apply(ref_sites[,-1], 1, function(x) paste(sites, x, sep = '_', collapse = '_'))
    ref_sites$comb <- combin
    ref_sites <- ref_sites[, .(allele, comb)]
    return(ref_sites)
}

sites_recategorate <- function(hla_table, site_combinations){
    hla_patient_comb <- merge(hla_table, site_combinations, by = 'allele')
    hla_patient_comb_matrix <- data.table::dcast(hla_patient_comb, patient_id ~ comb, value.var = 'homo_hetero', fun.aggregate = sum)
    return(hla_patient_comb_matrix)
}

site_matrix_haplotype <- function(sites){
    ref_sites <- data.table::dcast(hla_haplotype_reference[Site_hla %in% sites], haplotype ~ Site_hla, value.var = 'AA')
    setcolorder(ref_sites, c('haplotype', sites))
    combin <- apply(ref_sites[,-1], 1, function(x) paste(colnames(ref_sites)[-1], x, sep= '_', collapse = '_'))
    ref_sites$comb <- combin
    ref_sites <- ref_sites[, .(haplotype, comb)]
    return(ref_sites)
}

sites_recategorate_haplotype <- function(hla_haplotype_table, site_combinations){
    hla_patient_haplotype_comb <- merge(hla_haplotype_table, site_combinations, by = 'haplotype')
    hla_patient_haplotype_comb_matrix <- data.table::dcast(hla_patient_haplotype_comb, patient_id ~ comb, value.var = 'homo_hetero', fun.aggregate = sum)
    return(hla_patient_haplotype_comb_matrix)
}


conditional_fun <- function(cdr3_hla_matrix, x_reduced_sites, conditional_sites){
    try({
        
        name_pair <- unique(cdr3_hla_matrix$pair)
        n_pcs <- 3 #  change when needed
        amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
        Y_matrix <- as.matrix(cdr3_hla_matrix %>% dplyr::select(any_of(amino_acids)))
        Y <- paste0('cbind(', paste(colnames(Y_matrix), collapse = ','),") ~")
        x_reduced <- paste(x_reduced_sites, collapse = '+')
        pcs <- paste0('PC', seq(1,n_pcs), collapse = '+')

        comb_sites <- paste(conditional_sites, collapse = '+')

        if ('group' %in% colnames(cdr3_hla_matrix) ){
            X_full <- paste(comb_sites, pcs, 'group', sep = '+')
            X_null <- paste(x_reduced, pcs, 'group', sep = '+')
        } else {
            X_full <- paste(comb_sites, pcs, sep = '+')
            X_null <- paste(x_reduced, pcs, sep = '+')
        }

        formula_full <- as.formula(str_c(Y, X_full))
        formula_null <- as.formula(str_c(Y, X_null))
        mod1 <- lm(formula_full, data = cdr3_hla_matrix)
        mod0 <- lm(formula_null, data = cdr3_hla_matrix)

        manova_results <- anova(mod0, mod1)

        if (sum(is.na(rowSums(coef(mod1)))) > 0){
        
            allele_to_exclude <- rownames(coef(mod1))[is.na(rowSums(coef(mod1)))]
            new_comb_sites <- paste(setdiff(conditional_sites, allele_to_exclude), collapse = '+')

        } 
        if (sum(is.na(rowSums(coef(mod0)))) > 0){
            allele_to_exclude <- rownames(coef(mod0))[is.na(rowSums(coef(mod0)))]
            new_x_reduced <- paste(setdiff(x_reduced_sites, allele_to_exclude), collapse = '+')
        }
        
        try({
            if ('group' %in% colnames(cdr3_hla_matrix) ){
                X_full <- paste(new_comb_sites, pcs, 'group', sep = '+')
                X_null <- paste(new_x_reduced, pcs, 'group', sep = '+')
            } else {
                X_full <- paste(new_comb_sites, pcs, sep = '+')
                X_null <- paste(new_x_reduced, pcs, sep = '+')
            }

            formula_full <- as.formula(str_c(Y, X_full))
            formula_null <- as.formula(str_c(Y, X_null))
            mod1_mvlm <- mvlm(formula_full, cdr3_hla_matrix)
            mod0_mvlm <- mvlm(formula_null, cdr3_hla_matrix)
            
            var_exp_full <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
            p_val_full <- mod1_mvlm$pv["Omnibus Effect",1]
            
            var_exp_null <- mod0_mvlm$pseudo.rsq["Omnibus Effect",1]
            p_val_null <- mod0_mvlm$pv["Omnibus Effect",1]
            
            variance_explained <- var_exp_full - var_exp_null
            names(variance_explained) <- 'var_expl'
            P_val <-  mapply(c, p_val_null, p_val_full)
            colnames(P_val) <- 'P_val'
            Omnibus <- mapply(c, var_exp_null, var_exp_full)
            colnames(Omnibus) <- 'Omnibus'
            
            mvlm_df <- data.table(Omnibus, model = c('null', 'full'),P_val, variance_explained)

        

        }, silent = TRUE)

        if (exists('mvlm_df')){
            manova_df <- cbind(data.table(manova_results), mvlm_df)
        }else{
            manova_df <- data.table(manova_results)
        }
    
        manova_df$pair <- name_pair
        return(manova_df)

    }, silent = TRUE)
    
}


mlm_fun_conditional <- function(cdr3_hla_matrix, dir_results, sites, n_pcs = 3, writing_mode = 'no'){
    
    very_significant <- 0.01/100000
    name_pair <- unique(cdr3_hla_matrix$pair)

    variables <- define_variables(cdr3_hla_matrix)
    hla_alleles <- variables[[1]]
    Y_matrix <- variables[[2]]

    formulas <- define_formulas(hla_alleles, Y_matrix, n_pcs = n_pcs, sites = sites, method = 'conditional')

    mod1 <- lm(formulas[[1]], data = cdr3_hla_matrix)
    mod0 <- lm(formulas[[2]], data = cdr3_hla_matrix)

    manova_results <- anova(mod0, mod1)
    
    tryCatch({
        
        if (sum(is.na(rowSums(coef(mod1)))) > 0){

            allele_to_exclude <- rownames(coef(mod1))[is.na(rowSums(coef(mod1)))]
            new_hla_alleles <- setdiff(hla_alleles, allele_to_exclude)
            new_sites <- setdiff(sites, allele_to_exclude)
            formulas <- define_formulas(new_hla_alleles, Y_matrix, n_pcs = 3, sites = new_sites, method = 'conditional')

        } 
        
        }, 
        error = function(e) {
            # Handle the error
            message("Error: ", e)
        },
        finally = {
            mod1_mvlm <- mvlm(formulas[[1]], cdr3_hla_matrix)
            mod0_mvlm <- mvlm(formulas[[2]], cdr3_hla_matrix)
            
            var_exp_full <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
            p_val_full <- mod1_mvlm$pv["Omnibus Effect",1]
            
            var_exp_null <- mod0_mvlm$pseudo.rsq["Omnibus Effect",1]
            p_val_null <- mod0_mvlm$pv["Omnibus Effect",1]
            
            variance_explained <- var_exp_full - var_exp_null
            names(variance_explained) <- 'var_expl'
            P_val <-  mapply(c, p_val_null, p_val_full)
            colnames(P_val) <- 'P_val'
            Omnibus = mapply(c, var_exp_null, var_exp_full)
            colnames(Omnibus) <- 'Omnibus'
            
            mvlm_df <- data.table(Omnibus, model = c('null', 'full'),P_val, variance_explained)
            manova_df <- cbind(data.table(manova_results), mvlm_df)
            manova_df$pair <- name_pair

        })
        
    if (writing_mode == 'yes'){
        mod1_summary <- summary(mod1)
        mod0_summary <- summary(mod0)
    
        dir_results_summary <- paste0(dir_results,'/mmlm_summary/')
        if (!file.exists(dir_results_summary)) {
            dir.create(dir_results_summary, recursive = TRUE)
            }

        mod1_mvlm_summary <- summary(mod1_mvlm)
        mod0_mvlm_summary <- summary(mod0_mvlm)
        
        dir_results_summary_mvlm <- paste0(dir_results,'/mvlm_summary/')
        if (!file.exists(dir_results_summary_mvlm)) {
            dir.create(dir_results_summary_mvlm, recursive = TRUE)
            }

        try({
            if (na.omit(manova_df)$Pr..F. < very_significant){
            
                lapply(summary(mod1), function(x){
                    capture.output(x, file=paste0(dir_results_summary,name_pair,'_mod1.txt'), append = TRUE)})
                lapply(summary(mod0), function(x){
                    capture.output(x, file=paste0(dir_results_summary,name_pair,'_mod0.txt'), append = TRUE)})
                lapply(summary(mod1_mvlm), function(x){
                    capture.output(x, file=paste0(dir_results_summary_mvlm,name_pair,'_mod1.txt'), append = TRUE)})
                lapply(summary(mod0_mvlm), function(x){
                    capture.output(x, file=paste0(dir_results_summary_mvlm,name_pair,'_mod0.txt'), append = TRUE)}) 
            } 
          
            }, silent = TRUE)
    }


    return(manova_df)                
         
}


mlm_fun_permutation <- function(cdr3_hla_matrix_pca){
    
    amino_acids <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') 
    name_pair <- unique(cdr3_hla_matrix_pca$pair)

    # for v3 I use another approach, I already removed correlated columns from hla_matrices
    #without_ref_allele <- cdr3_hla_matrix_pca %>% 
    #        summarise(across(contains('allele'), sum)) %>% t() %>% as.data.frame() %>% arrange(desc(V1)) %>% slice(-1) %>% rownames()
    without_ref_allele <- colnames(cdr3_hla_matrix_pca %>% select(contains('allele')))
    
    aa_in_matrix <- paste0(colnames(cdr3_hla_matrix)[colnames(cdr3_hla_matrix) %in% amino_acids], collapse = ',')
    
    
    formula_full <- as.formula(paste0(paste0('cbind(',aa_in_matrix,') ~'), paste0(without_ref_allele, collapse = '+'), '+PC1+PC2+PC3'))
    formula_null <- as.formula(paste0('cbind(',aa_in_matrix,') ~', 'PC1+PC2+PC3'))
    
    mod1 <- lm(formula_full, data = cdr3_hla_matrix_pca)
    mod0 <- lm(formula_null, data = cdr3_hla_matrix_pca)
    
    manova <- anova(mod1, mod0)
    
    mod1_mvlm <- mvlm(formula_full, cdr3_hla_matrix_pca)
    mod0_mvlm <- mvlm(formula_null, cdr3_hla_matrix_pca)
    
    var_exp_full <- mod1_mvlm$pseudo.rsq["Omnibus Effect",1]
    p_val_full <- mod1_mvlm$pv["Omnibus Effect",1]
    
    var_exp_null <- mod0_mvlm$pseudo.rsq["Omnibus Effect",1]
    p_val_null <- mod0_mvlm$pv["Omnibus Effect",1]
    
    variance_explained <- var_exp_full - var_exp_null
    names(variance_explained) <- 'var_expl'
    P_val <-  mapply(c, p_val_null, p_val_full)
    colnames(P_val) <- 'P_val'
    Omnibus = mapply(c, var_exp_null, var_exp_full)
    colnames(Omnibus) <- 'Omnibus'
    
    manova_df <- data.frame(manova) %>% 
        mutate(pair = name_pair)
    mvlm_df <- data.frame(Omnibus, model = c('null', 'full'),P_val, variance_explained)
    
    manova_df <- bind_cols(manova_df, mvlm_df)
    
    return(manova_df)
 
            
}



setting_allele_reference <- function(cdr3_hla_matrix){
    allele_columns <- grep("allele", names(cdr3_hla_matrix), value = TRUE)

    for (col in allele_columns) {
        cdr3_hla_matrix[[col]] <- as.factor(cdr3_hla_matrix[[col]])
        # Find the most frequent level
        most_frequent <- names(sort(table(cdr3_hla_matrix[[col]]), decreasing = TRUE))[1]
        
        # Set the most frequent level as the reference
        cdr3_hla_matrix[[col]] <- relevel(cdr3_hla_matrix[[col]], ref = most_frequent)
    }
    return(cdr3_hla_matrix)
}

remove_correlated_columns <- function(df_matrix, cutoff_value){

    cor_matrix <- cor(df_matrix)
    columns_to_remove <- findCorrelation(cor_matrix, cutoff = cutoff_value, exact = FALSE, names = TRUE)
    df_wo_correlations <- df_matrix[, c(columns_to_remove) := NULL]
    
    return(df_wo_correlations)
}

