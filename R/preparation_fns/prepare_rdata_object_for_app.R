# script set up the R image which is uploaded and hosted w/ the app

######## PREP R IMG CONTENTS ######## 
rm(list=ls())

wd = "/home/abrowne/shiny/mito_nuc_shiny/"
gtf_loc = "/home/abrowne/gtf/Homo_sapiens.GRCh38.97.gtf"
summary_brain_loc = "/home/abrowne/projects/GTEx_6p/Alan_mapped_GTEx_6p/processed_files/AlanMap_GTExBrain_corr_SPEARMANsummary_table.csv"
summary_controls_loc = "/home/abrowne/projects/GTEx_6p/control_tissue_data/controlssummary_table.csv"

setwd(wd)

library(rtracklayer)
library(fst)

pval_correct_summary_table = function(summary_table, pval_col_name="pval"){
  
  # get pvalue columns
  pval_cols = colnames(summary_table)[grepl(pval_col_name, colnames(summary_table))]
  
  for(col in pval_cols){
    #print(summary_table[1:4,col])
    summary_table[col] = p.adjust(summary_table[,col], method = "fdr", n = length(summary_table[,col]))
    #print(summary_table[1:4,col])
  }
  
  return(summary_table)
}

# read in data
summary_brain = read.csv(summary_brain_loc)
summary_controls = read.csv(summary_controls_loc)

# pval correct
summary_brain = pval_correct_summary_table(summary_brain)
summary_controls = pval_correct_summary_table(summary_controls)

# fix col names
if(length(unique(summary_controls$nuc_gene)) < length(unique(summary_controls$mt_gene))){
  colnames(summary_controls)[colnames(summary_controls) == 'nuc_gene'] = 'temp1'
  colnames(summary_controls)[colnames(summary_controls) == 'mt_gene'] = 'temp2'
  colnames(summary_controls)[colnames(summary_controls) == 'temp1'] = 'mt_gene'
  colnames(summary_controls)[colnames(summary_controls) == 'temp2'] = 'nuc_gene'
}

grch = as.data.frame(rtracklayer::import("/home/abrowne/gtf/Homo_sapiens.GRCh38.97.gtf"))
write.fst(grch, "./data/grch3897.fst")

write_csv(summary_brain, "./data/GTEx_brain_summary_table.csv")
write_csv(summary_controls, "./data/GTEx_control_tissues_summary_table.csv")

write.fst(summary_brain, "./data/GTEx_brain_summary_table.fst")
write.fst(summary_controls, "./data/GTEx_control_tissues_summary_table.fst")

rm(list=ls())
