
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
  
  ### preparing brain region data -----------------------------------------------------
  
  # filtering for gene R value and pivoting data
  summary_gene_corrs = summary_brain %>% 
    dplyr::filter(nuc_gene == gene) %>% 
    dplyr::select(matches("corr|gene")) %>% 
    tidyr::gather(key="region", value="corrs", -c("mt_gene", "nuc_gene", "gene_name")) %>% 
    dplyr::arrange(region)
  
  # getting pvalues 
  summary_gene_pval = summary_brain %>% 
    dplyr::filter(nuc_gene == gene) %>% 
    dplyr::select(matches("pval")) %>% 
    tidyr::gather(key="region", value="pval") %>% 
    dplyr::arrange(region)
  
  summary_gene = cbind(summary_gene_corrs %>% dplyr::select(-gene_name), pval=summary_gene_pval$pval)
  
  ### preparing control region data -----------------------------------------------------
  
  # filtering for gene R value and pivoting data
  ctrl_gene_corrs = summary_controls %>% 
    dplyr::filter(nuc_gene == gene) %>% 
    dplyr::select(matches("corr|gene")) %>% 
    tidyr::gather(key="region", value="corrs", -c("mt_gene", "nuc_gene"))
  
  # getting pvalues 
  ctrl_gene_pval = summary_controls %>% 
    dplyr::filter(nuc_gene == gene) %>% 
    dplyr::select(matches("pvals|gene")) %>% 
    tidyr::gather(key="region", value="pval", -c("mt_gene", "nuc_gene"))
  
  summary_gene_ctrl = cbind(ctrl_gene_corrs, ctrl_gene_pval["pval"])
  
  summary_gene = rbind(summary_gene, summary_gene_ctrl)
  
  # add sig indicator
  summary_gene = summary_gene %>% mutate(sig_indicator = ifelse(summary_gene$pval < 0.05 & summary_gene$pval > 0.01, '*',
                                                                ifelse(summary_gene$pval < 0.01 & summary_gene$pval > 0.001, '**',
                                                                       ifelse(summary_gene$pval < 0.001, '***', ''))))
  
  # tidying region names
  summary_gene$region = gsub("Brain", "", 
                             gsub("_corrs", "", 
                                  gsub("basalganglia", "BG", 
                                       gsub("_spearman", "", 
                                            gsub("muscle", "Muscle", 
                                                 gsub("heart", "Heart", summary_gene$region))))))
  
  # # converting mt_gene genes
  summary_gene$mt_gene = convert_sym_ens(summary_gene$mt_gene, same_length=TRUE)
  
  #plotting heatmap
  p = ggplot(summary_gene, aes(x=reorder(region, corrs), y=reorder(mt_gene, corrs), fill=corrs)) +
    theme_classic(base_size=11) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # rotates x axis
    ggtitle(gene_sym) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0, size=9)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits=c(-1,1)) +
    geom_text(aes(label=sig_indicator), size=3) +
    labs(y="mtDNA-encoded gene", x="", fill=expression(rho), subtitle = "* 0.01<p<0.05, ** 0.001<p<0.01, *** p<0.001")
  
  return(p)
  
}

# #test
# start.time = Sys.time()
# print(heatmap_of_gene_corrs(gene="SOD2", summary_brain=summary_brain, summary_controls=summary_controls)) #ENSG00000005194
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)

# this fn. alone takes just 4s



