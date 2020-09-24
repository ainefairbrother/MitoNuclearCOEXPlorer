# fn that generates a distribution plot of all gene-mt values (per region) with a red line for gene X to indicate where it lies on the distribution

genDistributionPlotWithGene = function(gene, summary_brain){
  
  library(tidyr)
  library(ggpubr)
  library(ggplot2)
  library(dplyr)
  
  setwd('/home/abrowne/shiny/test_shiny/')
  
  # importing convert fn. 
  source("./R/convert_sym_ens.R")
  
  if(grepl('ENS', gene)){
    gene_sym = convert_sym_ens(gene, load_genespace=F)
  }
  if(!grepl('ENS', gene)){
    gene_sym = gene
    gene = convert_sym_ens(gene, input_ENS=F, load_genespace=F)
  }
  
  if(gene %in% summary_brain$nuc_gene){
  
    # melt summary_brain
    summary_brain = summary_brain %>% dplyr::select(matches('corr|nuc_gene')) %>% pivot_longer(2:13, names_to='region', values_to='R_value')
  
    # tidy region names
    summary_brain$region = gsub('_corrs', '', gsub('Brain', '', gsub('basalganglia', 'BG', summary_brain$region)))
  
    # gene subset
    means = vector('list', length(unique(summary_brain$region)))
    plots = vector('list', length(unique(summary_brain$region)))
    for(i in 1:length(unique(summary_brain$region))){
      
      # importing convert fn. 
      source('./R/convert_sym_ens.R')
      
      # converting ens --> sym 
      gene_sym = convert_sym_ens(gene)
  
      region_label = unique(summary_brain$region)[i]
      
      region_df = summary_brain %>% dplyr::filter(region == region_label)
  
      gene_df = summary_brain %>% dplyr::filter(nuc_gene == gene) %>% filter(region == region_label)
      
      mean_mt_R = mean(gene_df$R_value)
      means[[i]] = mean_mt_R
      names(means)[[i]] = region_label
  
      #vertical_lines = filter_df$R_value
      gene_legend = paste(gene_sym, 'mean correlation value across 13 mtDNA genes')
  
      p = ggplot(data=region_df, aes(x=R_value)) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_density(fill='aquamarine3') +
        ggtitle(region_label) +
        xlim(-1,1) +
        xlab(expression(rho)) + 
        ylab('Density') +
        geom_vline(aes_q(
          xintercept = mean_mt_R, 
          colour=gene_legend), 
          linetype='dashed', size=0.5) +
        theme(legend.position = c(0.9,0.9)) +
        scale_color_manual(name = '', values = c(gene_legend= 'aquamarine4'))
  
      plots[[i]] = p
      names(plots)[[i]] = region_label
    }
  
    # sorting plots by mean(mt)
    means_sorted = means[order(unlist(means), decreasing = T)] 
    plots_sorted = plots[names(means_sorted)]
  
    return(ggpubr::ggarrange(plotlist=plots_sorted, common.legend=TRUE))
    
  } else{
    
    return()
  }

}

start.time = Sys.time()
genDistributionPlotWithGene(gene='ENSG00000197548', summary_brain=summary_brain)
end.time = Sys.time()
time.taken = end.time - start.time
print(time.taken)
