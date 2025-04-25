
hla_genes = c('A','C','B','DRB1','DQA1', 'DQB1', 'DPA1','DPB1')

for (gene_name in hla_genes){
    hla <- read.table(paste0('../../reference_data/hla_msa/',gene_name,'_prot.txt'), sep = ' ', header = FALSE, col.names = c('allele', 'aa'), skip = 1) 
    hla_reference <- hla %>% 
        separate(aa, into = paste0(seq(0, lengths(strsplit(hla$aa[1], '')))), sep = "") 
    hla_reference <- hla_reference[vapply(hla_reference, function(x) x[1] != '.', logical(1L))]
    colnames(hla_reference) <- c('allele', seq(0, ncol(hla_reference)-2))                                          
    hla_reference_variable_sites <- hla_reference[vapply(hla_reference, function(x) length(unique(x)) > 1, logical(1L))]
                                                         
    hla_reference_long <- hla_reference_variable_sites %>% 
        pivot_longer(names_to = 'site', values_to = 'AA', cols = -1) %>%
        filter(AA != '*') 
    write_tsv(hla_reference_long, paste0("../../reference_data/hla_msa/",gene_name, "_long.tsv"))
} 

                                                         
                                                         