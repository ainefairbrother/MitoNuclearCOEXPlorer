# fn that generates a distribution plot of all gene-mt values (per region) with a red line for gene X to indicate where it lies on the distribution

genDistributionPlotWithGene = function(gene, summary_brain){
  
  #setwd('/home/abrowne/shiny/test_shiny/')
  
  if(grepl('ENS', gene)){
    gene_sym = convert_sym_ens(gene, load_genespace=F)
  }
  if(!grepl('ENS', gene)){
    gene_sym = gene
    gene = convert_sym_ens(gene, input_ENS=F, load_genespace=F)
  }
  
  # cleaning data
  # melt summary_brain
  summary_brain = summary_brain %>% dplyr::select(matches('corr|nuc_gene')) %>% pivot_longer(2:13, names_to='region', values_to='R_value')
  
  # tidy region names
  summary_brain$region = gsub('_corrs', '', gsub('Brain', '', gsub('basalganglia', 'BG', summary_brain$region)))
  
  region_names = as.vector(unique(summary_brain$region))
  
  prod_region_plot = function(val_){
    
    region_df = summary_brain %>% dplyr::filter(region == val_)
    
    gene_df = summary_brain %>% dplyr::filter(nuc_gene == gene) %>% filter(region == val_)
    
    mean_mt_R = mean(gene_df$R_value)
    
    #vertical_lines = filter_df$R_value

    distPlot = ggplot(data=region_df, aes(x=R_value)) +
      theme_classic(base_size=11) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_density(fill='aquamarine3') +
      ggtitle(val_) +
      xlim(-1,1) +
      xlab(expression(rho)) +
      ylab('Density') +
      geom_vline(aes_q(
        xintercept = mean_mt_R,
        colour='mean correlation value across 13 mtDNA genes'),
        linetype='dashed', size=0.5) +
      theme(legend.position = c(0.9,0.9)) +
      scale_color_manual(name = '', values = c('mean correlation value across 13 mtDNA genes'= 'aquamarine4'),
                         labels = paste(gene_sym, 'mean correlation value across 13 mtDNA genes'))
    
    return(distPlot)
  }
  
  prod_region_mean = function(val_){
    
    region_df = summary_brain %>% dplyr::filter(region == val_)
    
    gene_df = summary_brain %>% dplyr::filter(nuc_gene == gene) %>% filter(region == val_)
    
    mean_mt_R = mean(gene_df$R_value)
    
    return(mean_mt_R)
  }
  
  # if gene is present, run both functions
  if(gene %in% summary_brain$nuc_gene){

    plots = mclapply(region_names, prod_region_plot, mc.cores=3)
    names(plots) = region_names
    means = mclapply(region_names, prod_region_mean, mc.cores=3)
    names(means) = region_names
    
    # sorting plots by mean(mt)
    means_sorted = means[order(unlist(means), decreasing = T)]
    plots_sorted = plots[names(means_sorted)]
    
    # plotting panel of region plots
    return(ggpubr::ggarrange(plotlist=plots_sorted, common.legend=TRUE))

  }
  
}

# start.time = Sys.time()
# genDistributionPlotWithGene(gene='ENSG00000197548', summary_brain=summary_brain)
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)

# 3 cores 49 secs

