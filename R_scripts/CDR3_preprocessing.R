source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/libraries_analysis.R')
source('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/scripts/R_functions/cdr3_functions.R')

#------preprocessing .parquet files, treating each seqience as unique-------

path_healthy <- '/work_ikmb/sukmb667/projects/cdr3-qtl/TRB_samples//raw_data/'
path_ibd <- '/work_ikmb/sukmb667/projects/cdr3-qtl/TCellData//IBD_CDR3_QTL/IBD_TCR_REP/'
path_out <- '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/data/cdr3/minimal_preprocessed/'
seq_depth <- 30000 # this threshold is used to remove low depth repertoires,
# which are not suitable for further analysis

if (!file.exists(path_out)) {
        dir.create(path_out, recursive = TRUE)
        }

for (f in list.files(path_healthy)){
    path <- paste0(path_healthy,f)
    sample_id <- paste0('H', unlist(strsplit(unlist(strsplit(f, '_'))[[2]], '\\.'))[[1]])
    cdr3_dt <- as.data.table(read_parquet(path)) 
    if (nrow(cdr3_dt) > seq_depth){   # to remove low depth reprtoires
        cdr3_dt_preproccesing(cdr3_dt, sample_id, path_out)
    }

}

for (f in list.files(path_ibd)){
    path <- paste0(path_ibd,f)
    sample_id <- paste0('I', unlist(strsplit(unlist(strsplit(f, '_'))[[2]], '\\.'))[[1]])
    cdr3_dt <- as.data.table(read_parquet(path)) 
    if (nrow(cdr3_dt) > seq_depth){ # to remove low depth reprtoires
        cdr3_dt_preproccesing(cdr3_dt, sample_id, path_out)
    }

}

first_file <- TRUE
for (indiv in list.files(path_out)){
    cdr3_dt_preproc <- as.data.table(read_parquet(paste0(path_out,indiv))) 
    fwrite(cdr3_dt_preproc, '../data/cdr3_preprocessed_healthy_and_ibd.tsv', 
           sep = '\t', quote = FALSE, row.names = FALSE, append = TRUE, 
           col.names = first_file)
    first_file <- FALSE
    }


