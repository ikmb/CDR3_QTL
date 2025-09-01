#hla_gene <- 'B'

files_to_rbind <- data.table(paths = c('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/with_groups/both/',
    '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/H/',
    '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/UC/', 
    '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/I/',
    '/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/conditional_analysis/CD/'), 
    group = c('ALL', 'HBD', 'UC', 'IBD', 'CD'))

resulting_cdr3qtls <- rbindlist(Map(function(x,y) {
    rbindlist(lapply(grep(paste0(hla_gene, '_'), 
    list.files(x, pattern = '.tsv', full.names = TRUE), value = TRUE), function(z) determine_pval_from_each_round(z, hla_gene, y)))
}, files_to_rbind$paths, files_to_rbind$group))

files_to_manovas <- data.table(
    paths = c('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/main_manova_as_in_Ishigaki.tsv', 
        paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/H/manova_results_H_',hla_gene,'.tsv'),
        paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/UC/manova_results_UC_',hla_gene,'.tsv'), 
        paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/I/manova_results_I_',hla_gene,'.tsv'),
        paste0('/work_ikmb/sukmb667/projects/cdr3-qtl/healthy_and_ibd/manova_results/CD/manova_results_CD_',hla_gene,'.tsv')), 
    group = c('ALL', 'HBD', 'UC', 'IBD', 'CD'))

resulting_first_cdr3qtls <- rbindlist(Map(function(x,y) { 
    determine_pval_from_main_manova(x, hla_gene, y)}, 
    files_to_manovas$paths, files_to_manovas$group))

resulting_cdr3qtls_for_plot <- rbind(resulting_first_cdr3qtls, resulting_cdr3qtls) %>%
    filter(Pvalue <= 0.05/20000) %>%
    mutate(Pvalue = ifelse(Pvalue == 0 | Pvalue == Inf , smallest_number, Pvalue))