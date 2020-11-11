
# Aine Fairbrother-Browne
# 2020

densMode = function(x){
  td = density(x)
  maxDens = which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}

test_list_for_enrichment = function(gene_list, iters, filt_bkg_for_protein_coding=F, return_which='plot'){
  
  # function to take in a list of user-inputted genes
  # determines enrichment of this gene set in high -ve and +ve mito relationships
  
  # split input string into list of genes
  gene_list_spl = strsplit(gene_list, ', ')[[1]]
  iters = as.double(iters)
  
  stopifnot(length(gene_list_spl) > 1)
  
  if(grepl('ENS', gene_list_spl[1])){
    gene_sym = convert_sym_ens(gene_list_spl, input_ENS=T)
  }
  if(!grepl('ENS', gene_list_spl[1])){
    gene_sym = gene_list_spl
    gene_list_spl = convert_sym_ens(gene_list_spl, input_ENS=F)
  }
  
  # test for gene list being protein coding
  grch_list = grch %>%
    dplyr::filter(gene_id %in% gene_list_spl) %>% 
    dplyr::select(c(gene_id, gene_biotype))
  
  # get region names and corresponding col labels 
  region_names = colnames(summary_brain %>%
                            dplyr::select(matches('corr')))
  
  colnames(summary_brain) = gsub('_corrs', '', 
                                 gsub('Brain', '', 
                                      gsub('basalganglia', 'BG', colnames(summary_brain))))
  region_names = gsub('_corrs', '', 
                      gsub('Brain', '', 
                           gsub('basalganglia', 'BG', region_names)))
  
  # define df to collect output 
  rn = c('gene_set_med', 'random_set_med', 'num_above_rand_med', 'num_below_rand_med', 'p', 'sig_indicator')
  emptydf = as.data.frame(matrix(0L, nrow=length(rn), ncol=length(region_names)))
  colnames(emptydf) = region_names
  rownames(emptydf) = rn
  
  # define p val cutoff
  pcutoff = 0.05/(length(region_names))
  
  # parallelise across regions, calculating the deviance of the list from random for each region
  prod_region_plot = function(region, return_which_sub){
    
    # subset the summary table for region in current loop
    summary_region = summary_brain %>% 
      dplyr::select(region, 'nuc_gene')
    
    # filtering background list for biotype - depending on checkbox selection
    # preparing background list
    if(filt_bkg_for_protein_coding == T){
      grch_filt = grch %>%
        dplyr::filter(gene_id %in% unique(summary_region$nuc_gene)) %>%
        dplyr::filter(gene_biotype == 'protein_coding')
      nuc_gene_list = as.character(unique(grch_filt$gene_id))
      
    } else {
      nuc_gene_list = as.character(unique(summary_region$nuc_gene))
    }
    
    # subset the region subset for genes in the gene panel
    summary_region_gene_panel_subset = summary_region %>% 
      dplyr::filter(nuc_gene %in% gene_list_spl)
    
    # get gene panel median
    pop_med = median(summary_region_gene_panel_subset[,1])
    
    # define empty array to hold medians for all random iterations
    med_arr = replicate(iters, 0)
    
    # define empty matrix to hold random vectors of length=gene list for each iteration 
    rand_arr = matrix(0L, nrow=iters, ncol=length(gene_list_spl)*13)
    
    # repeat random picking and median computation for n iters
    for(j in 0:iters){
      
      # select random nuclear genes of size of target gene list
      rand_selection = sample(x=nuc_gene_list, size=length(gene_list_spl), replace=F)
      
      # then get all pairs associated with the random selection
      get_all_pairs = summary_region %>% 
        filter(nuc_gene %in% rand_selection)
      
      # add check and kill fn if not true
      stopifnot(dim(get_all_pairs)[1] == length(gene_list_spl)*13)
      
      # collect rand vecs
      rand_vals = as.vector(get_all_pairs[,c(region)])
      rand_arr[j,] = rand_vals
      
      # add medians to array
      med_arr[j] = median(rand_vals)
    }
    
    # median value of random medians
    med_of_med_arr = median(med_arr)
    
    # calculate where list median lies in pop. of random list medians
    above_pop_med = sum(med_arr > pop_med)
    below_pop_med = sum(med_arr < pop_med)
    
    # assign values to df
    emptydf['gene_set_med', region] = pop_med
    emptydf['random_set_med', region] = med_of_med_arr
    emptydf['num_above_rand_med', region] = above_pop_med
    emptydf['num_below_rand_med', region] = below_pop_med
    
    # calculate p
    if(pop_med < med_of_med_arr){
      p = as.numeric((iters-above_pop_med)/iters)
      if(p == 0){ # instead of a pvalue of 0, it gives p<1/number of iters
        emptydf['p', region] = 1/iters
      } else{
        emptydf['p', region] = p
      }
    }
    
    if(pop_med > med_of_med_arr){
      p = as.numeric((iters-below_pop_med)/iters)
      if(p == 0){
        emptydf['p', region] = 1/iters
      } else{
        emptydf['p', region] = p
      }
    }
    
    p = as.numeric(emptydf['p', region])
    
    if(return_which_sub == 'plot'){
      
      if((emptydf['p', region] < 0.05) & (emptydf['p', region] > pcutoff)){
        emptydf['sig_indicator', region] = '*'
      }
      if(emptydf['p', region] < pcutoff){
        emptydf['sig_indicator', region] = '**'
      }
      
      # downsampling rand_arr to 10% to allow graphics to plot more efficiently 
      if((iters <= 1000) & (iters > 100)){
        rand_arr = rand_arr[sample(x=nrow(rand_arr), size=100, replace=FALSE), ]
      }
      if(iters == 10000){
        rand_arr = rand_arr[sample(x=nrow(rand_arr), size=500, replace=FALSE), ]
      }
      
      # make into longform df for plotting and transpose
      rand_arr = rand_arr %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        gather()
      
      pl = ggplot() +
        theme_minimal(base_size=10) +
        xlim(-1,1) + 
        ggtitle(region) +
        xlab(expression(rho)) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.key = element_rect(fill = "transparent", colour = "transparent"),
              legend.background = element_blank(),
              legend.box.background = element_blank(),
              legend.position = c(0.15,0.9)) +
        
        # plot the random distributions of gene sets for the current region
        geom_density(data=rand_arr, 
                     aes(x=value, group=key, colour='Random gene sets'), alpha=0.001, fill='lightgrey') +
        
        # plot a line to indicate the median r value of the median of all the iterations
        geom_vline(xintercept=as.numeric(med_of_med_arr, 
                                         colour='Random gene set median'), colour='aquamarine4', linetype='dashed', size=0.5) + 
        
        # plot the gene panel distribution for the current region 
        geom_density(data=summary_region_gene_panel_subset, 
                     aes_q(x=summary_region_gene_panel_subset[,region], colour='Target gene set'), alpha=0.5, fill='lightgrey') +
        
        # plot a line to indicate the median r value of the gene set
        geom_vline(xintercept=as.numeric(pop_med), 
                   colour='blue3', linetype='dashed', size=0.5) + 
        
        scale_colour_manual(name='', values=c('Target gene set'='blue','Random gene sets'='aquamarine3'))
      
      if(p == 1/iters){
        pl = pl + labs(subtitle = paste0(paste('P < ', p), 
                                         "\n", 
                                         paste('Median of all random sets = ',
                                               round(med_of_med_arr, 4)),
                                         "\n", 
                                         paste('Median of target set = ', 
                                               round(pop_med, 4))))
        
      }
      if(p != 1/iters){
        pl = pl + labs(subtitle = paste0(paste('P = ', p), 
                                         "\n", 
                                         paste('Median of all random sets = ',
                                               round(med_of_med_arr, 4)),
                                         "\n", 
                                         paste('Median of target set = ', 
                                               round(pop_med, 4))))
      }
      
      return(pl)
    }else{
      return(p)
    }
  }
  
  # RUN ANALYSIS AND RENDER PLOT
  print("Running analysis...")
  start.time = Sys.time()

  if(return_which == 'plot'){
    plot_vect = lapply(region_names, prod_region_plot, return_which_sub='plot')
    fig = ggpubr::ggarrange(plotlist=plot_vect, common.legend=F, ncol=3, nrow=4)
    end.time = Sys.time()
    time.taken = end.time - start.time
    print(paste("Runtime =", time.taken))
    return(ggpubr::annotate_figure(fig,
                                   bottom = text_grob("Analysis and visualisation provided by the MitoNuclearCOEXPlorer tool [https://snca.atica.um.es/MitoNuclearCOEXPlorer/]", color = "black",
                                                      hjust = 1, x = 1, face = "italic", size = 10)))
    
  }else{
    pval = lapply(region_names, prod_region_plot, return_which_sub='pval')
    names(pval) = region_names
    end.time = Sys.time()
    time.taken = end.time - start.time
    print(paste("Runtime =", time.taken))
    return(pval)
  }
  
}

# start.time = Sys.time()
# ls = "ALKBH1, C1QBP, CDK5RAP1"
# test_list_for_enrichment(ls, iters=100, T)
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)


