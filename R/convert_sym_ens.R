
convert_sym_ens = function(id_list, input_ENS=T, same_length=F, load_genespace=F){
  
  # library(BiocManager)
  # options(repos = BiocManager::repositories())
  
  # import libs
  library(plyr)
  library(parallel)
  
  setwd("/home/abrowne/shiny/test_shiny/")
  
  # define fn 
  get_conversion = function(gene_id){
    if(input_ENS == T){
      conversion = Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$gene_id == gene_id & Homo_sapiens.GRCh38.97$type == "gene"]$gene_name
      if(length(conversion) > 0 ){
        return(conversion)
      }
      else{
        return(gene_id)
      }
    }
    if(input_ENS == F){
      conversion = Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$gene_name %in% id_list & Homo_sapiens.GRCh38.97$type == "gene"]$gene_id
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
    Homo_sapiens.GRCh38.97 = rtracklayer::import("./data/Homo_sapiens.GRCh38.97.gtf")
  }
  
  if(same_length==T){
    converted_out_list = mclapply(id_list, get_conversion, mc.cores=50)
    return(unlist(converted_out_list))
  }
  if(same_length==F){
    if(input_ENS == T){
      return(Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$gene_id %in% id_list & Homo_sapiens.GRCh38.97$type == "gene"]$gene_name)
    }
    else{
      return(Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$gene_name %in% id_list & Homo_sapiens.GRCh38.97$type == "gene"]$gene_id)
    }
  }
}

