
######## GENERATE R IMG ######## 

rm(list=ls())

wd = "/home/abrowne/shiny/mito_nuc_shiny/"
setwd(wd)

library(rtracklayer)
library(fst)

summary_brain = read.fst("./data/GTEx_brain_summary_table.fst")
summary_controls = read.fst("./data/GTEx_control_tissues_summary_table.fst")
grch = read.fst("./data/grch3897.fst")

#grch = as.data.frame(rtracklayer::import("/home/abrowne/gtf/Homo_sapiens.GRCh38.97.gtf"))

save.image("./data/app_data.Rdata")

rm(list=ls())

