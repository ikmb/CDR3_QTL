source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries.R')
source('/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3-QTL_functions.R')
library(DBI)

#hla_features <- read_tsv('../data/hla_features_wo_duplicates_healthy_and_ibd_those_with_cdr3.tsv')
#pairs <- '../cdr3_hla_pairs_v3/'
#dir_pair_path <- paste0('../cdr3_hla_pairs_9PCs/')
#dir_results <- paste0('../manova_results_9PCs/')

pcs <- 3
dir_downsampling <- '/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/CDR3_downsampling_experiment/unique/'
dir_content <- list.files(dir_downsampling, full.names = TRUE)
directories <- dir_content[file.info(dir_content)$isdir]

conn <- dbConnect(RSQLite::SQLite(), "/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/CDR3_downsampling_experiment/unique/cdr3-downsampling-unique.db")


for (path_in in directories){
    sample_size <- tail(unlist(strsplit(path_in, '/')), 1)
    pairs <- paste0(path_in, '/manova_results/cdr3_hla_pairs/irt_freq_unique/')
    pairs_paths <- paste0(pairs, list.files(pairs))
    dir_results <- paste0(path_in, '/manova_results/')

    for (f in pairs_paths){
        
        pair_matrix <- fread(f)
        manova_df <- mlm_fun(pair_matrix, dir_results, n_pcs = pcs)
        manova_df$sample_size <- sample_size
        dbWriteTable(conn, 
            "/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/CDR3_downsampling_experiment/unique/cdr3-downsampling-unique.db", 
            manova_df, append = TRUE)
        #if (first){
        #    manova_df_all <- manova_df[0,]
        #    first <- FALSE
        #}
       # if (is.null(manova_df)){
         #   next
        #} else {
        #    manova_df_all <- rbind(manova_df_all, manova_df) 
        #    }  

    }
    }

# Export table to CSV
db_path <- "/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/CDR3_downsampling_experiment/unique/cdr3-downsampling-unique.db"
conn <- dbConnect(RSQLite::SQLite(), db_path)
tables <- dbListTables(conn)
df <- dbReadTable(conn,  tables)

output_path <- "/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/CDR3_downsampling_experiment/unique/Downsampling_manova_all.csv"
write.csv(df, output_path, row.names = FALSE)

dbDisconnect(conn)