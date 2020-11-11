
# Aine Fairbrother-Browne
# 2020

convert_sym_ens = function(id_list, input_ENS=T, same_length=F, load_genespace=F){
  
  library(parallel)
  cores = parallel::detectCores(all.tests = FALSE, logical = TRUE)
  
  # define fn 
  get_conversion = function(gene_id){
    if(input_ENS == T){
      conversion = grch[(grch$gene_id == gene_id) & (grch$type == "gene"),]$gene_name
      if(length(conversion) > 0 ){
        return(conversion)
      }
      else{
        return(gene_id)
      }
    }
    if(input_ENS == F){
      conversion = grch[(grch$gene_name %in% id_list) & (grch$type == "gene"),]$gene_id #Homo_sapiens.GRCh38.97
      if(length(conversion) > 0 ){
        return(conversion)
      }
      else{
        return(gene_id)
      }
    }
  }
  
  # conditionally load genespace - could be loaded priorly if many lists are being passed
  if(load_genespace == T){
    grch = read.fst("./data/grch3897.fst")
  }
  
  if(same_length==T){
    converted_out_list = mclapply(id_list, get_conversion, mc.cores=cores/4)
    return(unlist(converted_out_list))
  }
  if(same_length==F){
    if(input_ENS == T){
      return(grch[(grch$gene_id %in% id_list) & (grch$type == "gene"),]$gene_name)
    }
    else{
      return(grch[(grch$gene_name %in% id_list) & (grch$type == "gene"),]$gene_id)
    }
  }
}