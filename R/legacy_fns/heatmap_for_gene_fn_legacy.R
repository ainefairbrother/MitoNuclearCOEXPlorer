# fn. to gen heatmap to show mt-nuc correlations of a gene:
# x-axis = brain region, y-axis = mt-genes
# gene - ens id
# summary_brain - needs to be pre-imported
# include control tissues
# supply a significance indicator 
# convert mt genes

heatmap_of_gene_corrs = function(gene, summary_brain, summary_controls){
  
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  
  setwd("/home/abrowne/shiny/test_shiny/")
  
  # importing correlation summary datasets 
  if(!exists("summary_brain")){
    summary_brain = read.fst("./data/GTEx_brain_summary_table.fst")
  }
  if(!exists("summary_controls")){
    summary_controls = read.fst("./data/GTEx_control_tissues_summary_table.fst")
  }
  if(!exists("Homo_sapiens.GRCh38.97")){
    # importing gtf for name conversion
    Homo_sapiens.GRCh38.97 = rtracklayer::import("./data/Homo_sapiens.GRCh38.97.gtf")
  }
  
  if(gene %in% summary_brain$nuc_gene){
    
    # importing convert fn. 
    source("./R/convert_sym_ens.R")
    
    # converting ens --> sym 
    gene_sym = convert_sym_ens(gene)
    
    ### preparing brain region data -----------------------------------------------------
    
    # filtering for gene R value and pivoting data
    summary_gene_corrs = summary_brain %>% 
      dplyr::filter(nuc_gene == gene) %>% 
      dplyr::select(matches("corr|gene")) %>% 
      tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs")
    
    # getting pvalues 
    summary_gene_pval = summary_brain %>% 
      dplyr::filter(nuc_gene == gene) %>% 
      dplyr::select(matches("pval|gene")) %>% 
      tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval")
    
    ### preparing control region data -----------------------------------------------------
    if(gene %in% summary_controls$nuc_gene){
      
      # filtering for gene R value and pivoting data
      ctrl_gene_corrs = summary_controls %>% 
        dplyr::filter(nuc_gene == gene) %>% 
        dplyr::select(matches("corr|gene")) %>% 
        tidyr::pivot_longer(cols=contains("corr"), names_to="region", values_to="corrs")
      
      print(ctrl_gene_corrs[1:5,])
      
      # getting pvalues 
      ctrl_gene_pval = summary_controls %>% 
        dplyr::filter(nuc_gene == gene) %>% 
        dplyr::select(matches("pval|mt_gene")) %>% 
        tidyr::pivot_longer(cols=contains("pval"), names_to="region", values_to="pval")
      
      summary_gene = bind_cols(summary_gene_corrs, summary_gene_pval["pval"])
      ctrl_gene = bind_cols(ctrl_gene_corrs, ctrl_gene_pval["pval"])
      
      master = bind_rows(summary_gene, ctrl_gene)
      
      #master$pval = round(master$pval, 2)
      
      # add sig indicator 
      master = master %>% mutate(sig_indicator = ifelse(master$pval < 0.05 & master$pval > 0.01, '*', 
                                                        ifelse(master$pval < 0.01 & master$pval > 0.001, '**', 
                                                               ifelse(master$pval < 0.001, '***', ''))))
      
      # tidying region names
      master$region = gsub("Brain", "", gsub("_corr", "", gsub("basalganglia", "BG", master$region)))
      
      # converting mt_gene genes
      master$mt_gene_gene = convert_sym_ens(master$mt_gene, load_genespace=F, same_length=T)
    } else{
      
      master = bind_cols(summary_gene_corrs, summary_gene_pval["pval"])
      
      
      master = bind_rows(summary_gene, ctrl_gene)
      master = master %>% mutate(sig_indicator = ifelse(master$pval < 0.05 & master$pval > 0.01, '*', 
                                                        ifelse(master$pval < 0.01 & master$pval > 0.001, '**', 
                                                               ifelse(master$pval < 0.001, '***', ''))))
      
      # tidying region names
      master$region = gsub("Brain", "", gsub("_corr", "", gsub("basalganglia", "BG", master$region)))
      
      # converting mt_gene genes
      master$mt_gene_gene = convert_sym_ens(master$mt_gene_gene, load_genespace=F, same_length=T)
    }
    
    #plotting heatmap
    return(ggplot(master, aes(x=reorder(region, corrs), y=reorder(mt_gene_gene, corrs), fill=corrs)) +
             theme_classic() + 
             geom_tile() +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # rotates x axis
             ggtitle(gene_sym) +
             theme(plot.title = element_text(hjust = 0.5)) +
             scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits=c(-1,1)) +
             geom_text(aes(label=sig_indicator), size=3) +
             labs(y="mtDNA-encoded gene", x="", fill="R value"))
    
  }
}

heatmap_of_gene_corrs(gene="ENSG00000198695", summary_brain=summary_brain, summary_controls=summary_controls)
