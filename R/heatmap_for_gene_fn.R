
# Aine Fairbrother-Browne
# 2020

heatmap_of_gene_corrs = function(gene, summary_brain, summary_controls){
  
  # fn. to gen heatmap to show mt-nuc correlations of a gene
  
  if(grepl('ENS', gene)){
    gene_sym = convert_sym_ens(gene)
  }
  if(!grepl('ENS', gene)){
    gene_sym = gene
    gene = convert_sym_ens(gene, input_ENS=F)
  }
  
  if(length(gene)>1){
    for(g in gene){
      if(g %in% summary_brain$nuc_gene){
        gene=g
      }
    }
  }
  
  print(gene)
  print(gene_sym)
  
  ### preparing brain region data -----------------------------------------------------

  # get correlations 
  summary_gene_corrs = summary_brain %>%
    dplyr::filter(nuc_gene == gene) %>%
    dplyr::select(matches("corr|gene")) %>%
    tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs") %>%
    dplyr::arrange(region) %>% 
    dplyr::mutate(pair=paste0(nuc_gene, '_', mt_gene)) %>% 
    dplyr::mutate(region=gsub('_corrs', '', region))
  
  # adding pvalues
  summary_gene = summary_brain %>%
    dplyr::filter(nuc_gene == gene) %>%
    dplyr::select(matches("pval|gene")) %>% 
    tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval") %>%
    dplyr::arrange(region) %>% 
    dplyr::mutate(pair=paste0(nuc_gene, '_', mt_gene)) %>% 
    dplyr::mutate(region=gsub('_pvals', '', region)) %>% 
    dplyr::left_join(summary_gene_corrs, by=c('pair', 'nuc_gene', 'mt_gene', 'region', 'gene_name'))
  
  ### preparing control region data -----------------------------------------------------
  
  ctrl_gene_corrs = summary_controls %>%
    dplyr::filter(nuc_gene == gene) %>%
    dplyr::select(matches("corr|gene")) %>%
    dplyr::mutate(pair=paste0(nuc_gene, '_', mt_gene)) %>% 
    tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs") %>% 
    dplyr::mutate(region=gsub('_spearman_corrs','',gsub('heart', 'Heart', gsub('muscle', 'Muscle', region)))) 
    
  master = summary_controls %>%
    dplyr::filter(nuc_gene == gene) %>%
    dplyr::select(matches("pvals|gene")) %>%
    dplyr::mutate(pair=paste0(nuc_gene, '_', mt_gene)) %>% 
    tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval") %>% 
    dplyr::mutate(region=gsub('_spearman_pvals','',gsub('heart', 'Heart', gsub('muscle', 'Muscle', region)))) %>% 
    dplyr::left_join(ctrl_gene_corrs, by=c('pair', 'nuc_gene', 'mt_gene', 'region')) %>% 
    dplyr::bind_rows(summary_gene) %>% 
    dplyr::mutate(sig_indicator = ifelse(pval < 0.05 & pval > 0.01, '*',
                                         ifelse(pval < 0.01 & pval > 0.001, '**',
                                                ifelse(pval < 0.001, '***', '')))) %>% 
    dplyr::mutate(region=gsub('basalganglia', 'BG', region)) %>% 
    mutate(mt_gene_sym=recode(mt_gene,
                              "ENSG00000198888"="MT-ND1",
                              "ENSG00000198763"="MT-ND2",
                              "ENSG00000198840"="MT-CO1",
                              "ENSG00000198886"="MT-CO2",
                              "ENSG00000212907"="MT-ATP8",
                              "ENSG00000198786"="MT-ATP6",
                              "ENSG00000198695"="MT-CO3",
                              "ENSG00000198899"="MT-ND3",
                              "ENSG00000228253"="MT-ND4L",
                              "ENSG00000198804"="MT-ND4",
                              "ENSG00000198712"="MT-ND5",
                              "ENSG00000198938"="MT-ND6",
                              "ENSG00000198727"="MT-CYB"))

  # #plotting heatmap
  return(ggplot(master, aes(x=reorder(region, corrs), y=reorder(mt_gene_sym, corrs), fill=corrs)) +
           theme_minimal(base_size=11) +
           theme(plot.title = element_text(hjust = 0.5)) +
           geom_tile() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # rotates x axis
           ggtitle(paste0(gene_sym,' | ', gene)) +
           theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0, size=9)) +
           scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits=c(-1,1)) +
           geom_text(aes(label=sig_indicator), size=3) +
           labs(y="mtDNA-encoded gene", x="", fill=expression(rho), subtitle = "* 0.01<p<0.05, ** 0.001<p<0.01, *** p<0.001"))
}

# # #test
# gene = "ENSG00000197548"
# start.time = Sys.time()
# heatmap_of_gene_corrs(gene=gene, summary_brain=summary_brain, summary_controls=summary_controls) #ENSG00000005194
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)


