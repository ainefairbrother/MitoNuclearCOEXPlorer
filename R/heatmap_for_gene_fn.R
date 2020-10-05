# fn. to gen heatmap to show mt-nuc correlations of a gene:
# x-axis = brain region, y-axis = mt-genes
# gene - ens id
# summary_brain - needs to be pre-imported
# include control tissues
# supply a significance indicator 
# convert mt genes

heatmap_of_gene_corrs = function(gene, summary_brain, summary_controls){
  
  if(grepl('ENS', gene)){
    gene_sym = convert_sym_ens(gene)
  }
  if(!grepl('ENS', gene)){
    gene_sym = gene
    gene = convert_sym_ens(gene, input_ENS=F)
  }
  
  if(gene %in% summary_brain$nuc_gene){
    
    ### preparing brain region data -----------------------------------------------------
    
    # filtering for gene R value and pivoting data
    summary_gene_corrs = summary_brain %>% 
      dplyr::filter(nuc_gene == gene) %>% 
      dplyr::select(matches("corr|gene")) %>% 
      tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs") %>% 
      dplyr::arrange(region)
    
    # getting pvalues 
    summary_gene_pval = summary_brain %>% 
      dplyr::filter(nuc_gene == gene) %>% 
      dplyr::select(matches("pval|gene")) %>% 
      tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval") %>% 
      dplyr::arrange(region)
    
    summary_gene = dplyr::bind_cols(list(summary_gene_corrs, summary_gene_pval[,'pval']))
    
    # df structure test
    if(dim(summary_gene)[1] != length(unique(summary_gene$region))*13){
      print(
        paste("Incorrect row number for summary_gene. Table does not contain a per-mt gene pair for", gene, "for every tissue")
      )
      #break()
    }
    
    ### preparing control region data -----------------------------------------------------
    if(gene %in% summary_controls$nuc_gene){
      
      # filtering for gene R value and pivoting data
      ctrl_gene_corrs = summary_controls %>%
        dplyr::filter(nuc_gene == gene) %>%
        dplyr::select(matches("corr|gene")) %>%
        tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs")
      
      # getting pvalues
      ctrl_gene_pval = summary_controls %>%
        dplyr::filter(nuc_gene == gene) %>%
        dplyr::select(matches("pval|gene")) %>%
        tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval")
      
      summary_gene_ctrl = dplyr::bind_cols(list(ctrl_gene_corrs, ctrl_gene_pval[,'pval']))
      
      summary_gene = summary_gene %>% dplyr::bind_rows(summary_gene, summary_gene_ctrl)
    }
    
    
    # add sig indicator
    summary_gene = summary_gene %>% mutate(sig_indicator = ifelse(summary_gene$pval < 0.05 & summary_gene$pval > 0.01, '*',
                                                      ifelse(summary_gene$pval < 0.01 & summary_gene$pval > 0.001, '**',
                                                             ifelse(summary_gene$pval < 0.001, '***', ''))))
    
    # tidying region names
    summary_gene$region = gsub("Brain", "", 
                               gsub("_corrs", "", 
                                    gsub("basalganglia", "BG", 
                                         gsub("_spearman", "", summary_gene$region))))
    
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
    
  } else {
    return("")
  }
}

#test
#heatmap_of_gene_corrs(gene="PINK1", summary_brain=summary_brain, summary_controls=summary_controls) #ENSG00000005194
