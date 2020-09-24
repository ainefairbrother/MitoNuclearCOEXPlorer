# function to take in a list of user-inputted genes
# determines enrichment of this gene set in high -ve and +ve mito relationships

densMode = function(x){
  td = density(x)
  maxDens = which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}

test_list_for_enrichment = function(gene_list, iters=10000){
  
  setwd("/home/abrowne/shiny/test_shiny/")
  
  # import libs
  library(sjPlot)
  library(ggplot2)
  library(gridExtra)
  
  # split input string into list of genes
  #gene_list = strsplit(gene_list, ", ")[[1]]
  
  # get region names and corresponding col labels 
  region_names = gsub("_corr", "", gsub("Brain", "", colnames(summary_brain)))[3:14]
  colnames(summary_brain)[3:14] = region_names
  
  # define df to collect output 
  rn = c('gene_set_med', 'random_set_med', 'num_above_rand_med', 'num_below_rand_med', 'p', 'sig_indicator')
  emptydf = as.data.frame(matrix(0L, nrow=length(rn), ncol=length(region_names)))
  colnames(emptydf) = region_names
  rownames(emptydf) = rn
  
  # define p val cutoff
  pcutoff = 0.05/(length(region_names))
  
  # define vect to hold plots 
  plot_vect = list()
  
  # loop through regions, calculating the deviance of the list from random for each region 
  for(i in 1:length(region_names)){
    
    # subset the summary table for region in current loop
    summary_region = summary_brain[,c(region_names[i], "nuc_gene")]
    
    # get list of nuclear genes
    nuc_gene_list = unique(summary_region$nuc_gene)
    
    # subset the region subset for genes in the gene panel
    summary_region_gene_panel_subset = summary_region[summary_region$nuc_gene %in% gene_list, ]
    
    # get gene panel median
    pop_med = median(summary_region_gene_panel_subset[,1])
    
    # define empty array to hold medians for all random iterations
    med_arr = replicate(iters, 0)
    
    # define empty matrix to hold random vectors of length=gene list for each iteration 
    rand_arr = matrix(0L, nrow=iters, ncol=length(gene_list)*13)
    
    # repeat random picking and median computation for n iters
    for(j in 0:iters){
      
      rand_selection = sample(nuc_gene_list, length(gene_list))
      get_all_pairs = summary_region %>% 
        filter(nuc_gene %in% rand_selection)
      
      # add check and kill fn if not true
      stopifnot(dim(get_all_pairs)[1] == length(gene_list)*13)
      
      # collect rand vecs
      rand_vals = as.vector(get_all_pairs[,c(region_names[i])])
      rand_arr[j,] = rand_vals
      
      # add medians to array
      med_arr[j] = median(rand_vals)
    }
    
    rand_arr = as.data.frame(rand_arr)
    rand_arr = rand_arr %>% gather()
    
    # median value of random medians
    med_of_med_arr = median(med_arr)
    
    # calculate where list median lies in pop. of random list medians
    above_pop_med = sum(med_arr > pop_med)
    below_pop_med = sum(med_arr < pop_med)
    
    # assign values to df
    emptydf["gene_set_med", region_names[i]] = pop_med
    emptydf["random_set_med", region_names[i]] = med_of_med_arr
    emptydf["num_above_rand_med", region_names[i]] = above_pop_med
    emptydf["num_below_rand_med", region_names[i]] = below_pop_med
    
    # calculate p
    if(pop_med < med_of_med_arr){
      p = as.numeric((iters-above_pop_med)/iters)
      if(p == 0){
        emptydf["p", region_names[i]] = 1/iters
      } else{
        emptydf["p", region_names[i]] = p
      }
    }
    
    
    if(pop_med > med_of_med_arr){
      p = (iters-below_pop_med)/iters
      print(p)
      if(p == 0){
        emptydf["p", region_names[i]] = 1/iters
      } else{
        emptydf["p", region_names[i]] = p
      }
    }
    
    if((emptydf["p", region_names[i]] < 0.05) & (emptydf["p", region_names[i]] > pcutoff)){
      emptydf["sig_indicator", region_names[i]] = "*"
    }
    if(emptydf["p", region_names[i]] < pcutoff){
      emptydf["sig_indicator", region_names[i]] = "**"
    }
    
    p = as.numeric(emptydf["p", region_names[i]])
    
    # Density plots with semi-transparent fill
    pl = ggplot() +
      theme_classic(base_size=8) +
      xlim(-1,1) + 
      labs(title=region_names[i]) +
      
      # plot the random distributions of gene sets for the current region
      geom_density(data=rand_arr, aes(x=value, group=key, colour='Random gene sets'), alpha=0.001, fill='lightgrey') +
      # plot a line to indicate the median r value of the median of all the iterations
      geom_vline(xintercept=as.numeric(med_of_med_arr), colour='aquamarine4', linetype='dashed', size=0.5) + 
      
      # plot the gene panel distribution for the current region 
      geom_density(data=summary_region_gene_panel_subset, 
                   aes_q(x=summary_region_gene_panel_subset[,region_names[i]], colour='Target gene set'), alpha=0.5, fill='lightgrey') +
      
      # plot a line to indicate the median r value of the gene set
      geom_vline(xintercept=as.numeric(pop_med), colour='blue3', linetype='dashed', size=0.5) + 
      
      # add stats to plot
      annotate("text", x=0.9, y=(densMode(rand_arr[, "value"])$y)*2, hjust=1, label=paste("Gene set median = ",   round(pop_med, 4)),       size=2.5) +
      annotate("text", x=0.9, y=(densMode(rand_arr[, "value"])$y)*2-0.2, hjust=1, label=paste("Random set median = ", round(med_of_med_arr, 4)), size=2.5) +
      theme(legend.position = c(0.15,0.9)) +
      scale_colour_manual(name="", values=c("Target gene set"="blue","Random gene sets"="aquamarine3"))
    
    if(p == 1/iters){
      pl = pl + annotate("text", x=0.9, y=(densMode(rand_arr[, "value"])$y)*2+0.2, hjust=1, label=paste("P < ", p), size=2.5)
    }
    if(p != 1/iters){
      pl = pl + annotate("text", x=0.9, y=(densMode(rand_arr[, "value"])$y)*2+0.2, hjust=1, label=paste("P = ", p), size=2.5)
    }
    
    plot_vect[[i]] = pl
  }
  
  #return(plot_vect[[1]])
  
  return(sjPlot::plot_grid(x=plot_vect, margin=c(0.1, 0.1, 0.1, 0.1), tags=replicate(length(region_names), "")))
  
}

# test
ls = read.csv("/home/abrowne/Gene_lists/panel_app_lists/Parkinson Disease and Complex Parkinsonism.tsv", sep="\t")
ls = as.vector(ls$EnsemblId.GRch38.)
test_list_for_enrichment(ls, iters=1000)







