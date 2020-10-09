
# Aine Fairbrother-Browne
# 2020

GenFigSingleGene = function(arg1){
  
  # start.time = Sys.time()
  
  a = genDistributionPlotWithGene(arg1, summary_brain=summary_brain)
  
  b = heatmap_of_gene_corrs(arg1, summary_brain, summary_controls=summary_controls)

  p = ggpubr::ggarrange(b,a,nrow=2,ncol=1,labels=c("A", "B"))
  
  # end.time = Sys.time()
  # time.taken = end.time - start.time
  # print(time.taken)

  return(p)
  
}

# start.time = Sys.time()
# GenFigSingleGene('ENSG00000197548')
# end.time = Sys.time()
# time.taken = end.time - start.time
# print(time.taken)
# 1.14 mins