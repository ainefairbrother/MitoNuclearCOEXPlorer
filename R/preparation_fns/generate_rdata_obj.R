
######## GENERATE R IMG ######## 

rm(list=ls())

wd = "/home/abrowne/shiny/test_shiny/"
setwd(wd)

library(rtracklayer)
library(fst)

summary_brain = read.fst("./data/GTEx_brain_summary_table.fst")
summary_contols = read.fst("./data/GTEx_control_tissues_summary_table.fst")

Homo_sapiens.GRCh38.97 = rtracklayer::import("/home/abrowne/gtf/Homo_sapiens.GRCh38.97.gtf")
grch = as.data.frame(Homo_sapiens.GRCh38.97)

# # importing fns
# source("./R/convert_sym_ens.R")
# source("./R/gen_gene_specific_distr_plot_fn.R")
# source("./R/single_gene_gen_fig.R")
# source("./R/test_list_for_enrichment_fn.R")

save.image("./app_data.Rdata")

rm(list=ls())

