
# Aine Fairbrother-Browne
# 2020

GenFigSingleGene = function(arg1){
  
  # start.time = Sys.time()
  
  a = genDistributionPlotWithGene(arg1, summary_brain=summary_brain)
  
  b = heatmap_of_gene_corrs(arg1, summary_brain, summary_controls=summary_controls)
  
  p = ggpubr::ggarrange(b,a,nrow=2,ncol=1,labels=c("A", "B"))
  
  return(ggpubr::annotate_figure(p,
                                 bottom = text_grob("Plot produced using the MitoNuclearCOEXPlorer tool [https://snca.atica.um.es/MitoNuclearCOEXPlorer/]", color = "black",
                                                    hjust = 1, x = 1, face = "italic", size = 10)))
}

# start.time = Sys.time()
# GenFigSingleGene('ATG7')
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)
# 19s