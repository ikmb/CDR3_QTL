source('libraries.R')
source('cdr3_functions.R')

#----------function with QC and preprocessing of TCR-seq data------
cdr3_preprocessing <- function(path_to_parquet_file, sample_id, path_out){
    
    cdr3_df <- read_parquet(path_to_parquet_file)  %>% 
        filter(sequenceStatus == 'In') %>%
        select(c(aminoAcid, count, vFamilyName, vGeneName, jFamilyName, jGeneName)) %>% 
        mutate(patient_id = sample_id, count = as.integer(count)) %>%
        rename(amino_acids = aminoAcid, v_gene = vGeneName, j_gene = jGeneName) %>%
        mutate(length = nchar(amino_acid)) %>%
        separate(v_gene, c('v_gene', 'v_gene_2'), sep = '/') %>%  
        mutate(v_gene = ifelse(is.na(v_gene), vFamilyName, v_gene)) %>%
        mutate(j_gene = ifelse(is.na(j_gene), jFamilyName, j_gene)) %>%
        select(patient_id, amino_acid, count, v_gene, j_gene, length) %>%
        group_by(patient_id, amino_acid, v_gene, j_gene, length) %>%
        summarise(count = sum(count), .groups = 'keep') %>% 
        ungroup()

    write_parquet(cdr3_df, paste0(path_out,'preprocessed_',id,'.parquet'))    
    return(cdr3_df)
}

#---------------------------------------------------------------------

path_healthy <- '../../TRB_samples//raw_data/'
path_ibd <- '../../TCellData//IBD_CDR3_QTL/IBD_TCR_REP/'
path_out <- '../data/preprocessed_v2/'

for (f in list.files(path_healthy)){
    path <- paste0(path_healthy,f)
    patient_id <- paste0('H', unlist(strsplit(unlist(strsplit(f, '_'))[[2]], '\\.'))[[1]])
    cdr3_dt <- data.table(read_parquet(path)) 
    cdr3_preprocessing(path, patient_id, path_out)
    summary_cdr3_dt <- rbind(summary_cdr3_dt, 
                             data.table(patient_id = patient_id, 
                                        total_number_of_seq = cdr3_dt[,.N], 
                                        productive_CDR3 = cdr3_dt[sequenceStatus == 'In',.N], 
                                        n_unique = uniqueN(cdr3_dt[sequenceStatus == 'In',aminoAcid])))
}

for (f in list.files(path_ibd)){
    path <- paste0(path_ibd,f)
    patient_id <- paste0('I', unlist(strsplit(unlist(strsplit(f, '_'))[[2]], '\\.'))[[1]])
    cdr3_dt <- data.table(read_parquet(path)) 
    cdr3_preprocessing(path, patient_id, path_out)
    summary_cdr3_dt <- rbind(summary_cdr3_dt, 
                             data.table(patient_id = patient_id, 
                                        total_number_of_seq = cdr3_dt[,.N], 
                                        productive_CDR3 = cdr3_dt[sequenceStatus == 'In',.N], 
                                        n_unique = uniqueN(cdr3_dt[sequenceStatus == 'In',aminoAcid])))
}

for (indiv in list.files(path_out)){
    cdr3_df <- read_parquet(paste0(path_out,indiv)) 
    write.table(cdr3_df, '../data/cdr3_preprocessed_healthy_and_ibd.tsv', sep = '\t', quote = FALSE, row.names = FALSE, append = TRUE, col.names = FALSE)
    }

write_tsv(summary_cdr3_dt, '../data/summary_cdr3_dt.tsv')

ids_low_repertoire <- summary_cdr3_dt  %>% filter(total_number_of_seq < 30000) %>% pull(patient_id)

write.table(ids_low_repertoire, '../data/ids_low_depths.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

#---Read preprocessed CDR3, removing duplicated samples, and filtering lengths 12 - 18----

cdr3 <- data.table(read_tsv('../data/cdr3_preprocessed_healthy_and_ibd.tsv', col_names = FALSE))

colnames(cdr3) <- c('patient_id', 'amino_acids', 'v_gene', 'j_gene', 'length', 'count')

write_tsv(cdr3, '../data/cdr3_preprocessed_healthy_and_ibd.tsv')
ids_in_cdr3 <- unique(pull(cdr3_qc, patient_id))
write.table(ids_in_cdr3, '../data/ids_in_cdr3.txt', col.names = FALSE, row.names = FALSE, quote = FALSE)


ids_duplicates <- read.table('../data/ids_to_remove_duplicates.txt')$V1

cdr3_qc <- cdr3 %>% 
    filter(!(patient_id %in% ids_low_repertoire)) %>%
    filter(!patient_id %in% ids_duplicates) %>%
    filter(length >= 12 & length<=18)

write_tsv(cdr3_qc, '../data/cdr3_qc_healthy_and_ibd.tsv')

cdr3_grouped <- cdr3_qc %>% group_by(v_gene, j_gene, length) %>% group_split()
for (i in seq(1,length(cdr3_grouped))){
    cdr3_group <- cdr3_grouped[[i]]
    v_j_exclude(cdr3_group)
    }






