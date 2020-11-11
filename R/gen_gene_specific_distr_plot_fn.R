
# Aine Fairbrother-Browne
# 2020

genDistributionPlotWithGene = function(gene, summary_brain){
  
  # fn that generates a distribution plot of all gene-mt values (per region) with a red line for gene X to indicate where it 
  # lies on the distribution
  
  if(grepl('ENS', gene)){
    gene_sym = convert_sym_ens(gene, load_genespace=F)
  }
  if(!grepl('ENS', gene)){
    gene_sym = gene
    gene = convert_sym_ens(gene, input_ENS=F, load_genespace=F)
  }
  
  # filtering for gene R value and pivoting data
  summary_brain_clean = summary_brain %>% 
    dplyr::select(matches("corr|nuc_gene")) %>% 
    tidyr::gather(key="region", value="R_value", -c("nuc_gene")) %>% 
    dplyr::mutate(region = gsub('_corrs', '', gsub('Brain', '', gsub('basalganglia', 'BG', region))))
    
  gene_mean_across_mt_df = summary_brain_clean %>% 
    dplyr::filter(nuc_gene == gene) %>% 
    dplyr::group_by(region) %>%
    dplyr::mutate(gene_mean_across_mt = round(mean(R_value), 4))
  
  gene_mean_label_df = data.frame(
    region = unique(gene_mean_across_mt_df$region),
    label = unique(gene_mean_across_mt_df$gene_mean_across_mt))
  
  distPlot = ggplot(data=summary_brain_clean, aes(x=R_value)) +
    theme_classic(base_size=11) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="top") +
    geom_density(fill='#ECDEFF', colour='white') +
    facet_wrap(facets=vars(region), nrow=4, ncol=3) +
    xlim(-1,1) +
    xlab(expression(rho)) +
    ylab('Density') +
    
    # add line to indicate mean rho of gene across 13 mt genes
    geom_vline(
      data = gene_mean_label_df,
      mapping = aes(
      xintercept = label,
      colour='mean correlation value across 13 mtDNA genes'),
      linetype='dashed',
      size=0.5) +
    
    # add annotation to distplot: mean rho
    geom_label(
      size=3,
      data = gene_mean_label_df,
      mapping = aes(label=paste0("\U03BC", '=', label),
                    x = as.numeric(label),
                    y = 0.3),
      fill='darkgrey',
      colour = 'white',
      fontface = 'bold'
    )  +
    
    # add legend 
    scale_color_manual(name = '', values = c('mean correlation value across 13 mtDNA genes'= '#3016FF'),
                       labels = paste(gene_sym, 'mean correlation value across 13 mtDNA genes'))

  return(distPlot)
  
}

# start.time = Sys.time()
# genDistributionPlotWithGene(gene='ATG7', summary_brain=summary_brain)
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)

