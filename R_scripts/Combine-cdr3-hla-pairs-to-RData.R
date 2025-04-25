
#path_in <- "/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/cdr3_hla_pairs/"
#directories <- list.files("/work_beegfs/sukmb667/projects/cdr3-qtl/healthy_and_ibd/cdr3_hla_pairs/", full.names = TRUE)
#pairs_dir <- directories[file.info(directories)$isdir]

for (tsv_dir in pairs_dir) {
    prep_mode <- strsplit(tsv_dir, "/")[[1]][length(strsplit(tsv_dir, "/")[[1]])]
    # List all .tsv files in the directory
    tsv_files <- list.files(tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

    # Read all .tsv files into a list of data.tables
    data_table_list <- lapply(tsv_files, fread)

    save(data_table_list, 
        file = paste0(path_in,prep_mode,".RData"))
}