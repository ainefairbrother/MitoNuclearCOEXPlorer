# shiny interface to allow querying of correlation data

# ---------------------- imports ---------------------------------------

# allows bioconductor pkgs to be found by shiny.io server

#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
options(warn=-1)

setwd("/home/abrowne/shiny/test_shiny/")

library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(rtracklayer)
library(fst)

#load("./shiny_app_env_image.Rdata")

# importing fns
source("/home/abrowne/shiny/test_shiny/R/heatmap_for_gene_fn.R")
source("/home/abrowne/shiny/test_shiny/R/gen_gene_specific_distr_plot_fn.R")
source("/home/abrowne/shiny/test_shiny/R/test_list_for_enrichment_fn.R")

# ---------------------- shiny setup ---------------------------------
# user interface
ui = fluidPage(
  
  # anything put as an argument into fluidPage will be  visible in the app 
  titlePanel('Mito-nuc brain browser'), 
  #headerPanel('Query mitocondrial-nuclear relationships by nuclear gene'), 
  
  # sidebarPanel layout aesthetic
  sidebarPanel(
    
    # Text at the top of the side panel
    helpText("Enter a Ensembl gene code to query \n 
             the mitochondrial-nuclear correlation data"), 
    
    # input box to insert single gene
    textInput(inputId = "gene_text", h3("Gene ID"), value = ""), #ENSG00000197548
    
    # input box to insert gene list
    textInput(inputId = "gene_list_text", h3("List of gene IDs"), value = c(""))
    #ENSG00000101986,ENSG00000141385, ENSG00000003393, ENSG00000214274, ENSG00000122359, ENSG00000142192, ENSG00000100299, ENSG00000142192, ENSG00000105409, ENSG00000123191, ENSG00000148090, ENSG00000131943, ENSG00000006283, ENSG00000162063, ENSG00000250479, ENSG00000106153, ENSG00000083937, ENSG00000114859, ENSG00000128973, ENSG00000068120, ENSG00000047457, ENSG00000182578, ENSG00000174080, ENSG00000135929, ENSG00000172817, ENSG00000117593, ENSG00000204843, ENSG00000101152, ENSG00000116675, ENSG00000130816, ENSG00000111361, ENSG00000119718, ENSG00000070785, ENSG00000115211, ENSG00000145191, ENSG00000118402, ENSG00000112425, ENSG00000100225, ENSG00000112367, ENSG00000087086, ENSG00000089280, ENSG00000131979, ENSG00000131095, ENSG00000030582, ENSG00000213614, ENSG00000049860, ENSG00000135486, ENSG00000166033, ENSG00000136156, ENSG00000131398, ENSG00000171385, ENSG00000164976, ENSG00000155980, ENSG00000188906, ENSG00000143669, ENSG00000186868, ENSG00000187566, ENSG00000074181, ENSG00000141458, ENSG00000119655, ENSG00000123240, ENSG00000125779, ENSG00000116288, ENSG00000100311, ENSG00000113721, ENSG00000108518, ENSG00000158828, ENSG00000184381, ENSG00000185345, ENSG00000171867, ENSG00000080815, ENSG00000143801, ENSG00000011275, ENSG00000107290, ENSG00000168575, ENSG00000145335, ENSG00000142168, ENSG00000021574, ENSG00000104133, ENSG00000161011, ENSG00000159082, ENSG00000120948, ENSG00000183735, ENSG00000205090, ENSG00000095970, ENSG00000011295, ENSG00000011600, ENSG00000188021, ENSG00000124164, ENSG00000165280, ENSG00000197969, ENSG00000069329, ENSG00000196998, ENSG00000143324, ENSG00000169083, ENSG00000111676, ENSG00000130638, ENSG00000124788, ENSG00000204842, ENSG00000066427, ENSG00000163635, ENSG00000147894, ENSG00000141837, ENSG00000160213, ENSG00000165060, ENSG00000197386, ENSG00000154118, ENSG00000101361, ENSG00000156475, ENSG00000112592
    
    ),
  
  mainPanel(
    tabsetPanel(type="tabs", 
                tabPanel("Correlation heatmap", plotOutput("heatmap")), 
                tabPanel("Region-wise R distribution", plotOutput("distributions")),
                tabPanel("Enrichment of gene list in mt-nuc pairs", plotOutput("enrichments"))
    ),
    
    width=5
  )
)

# backend
server = function(input, output){
  
  output$heatmap = renderPlot({
    
    heatmap_of_gene_corrs(gene=input$gene_text, summary_brain=summary_brain, summary_controls=summary_controls)
    
  },
  height = 800, 
  width = 1200
  )
  
  output$distributions = renderPlot({
    
    genDistributionPlotWithGene(gene=input$gene_text, summary_brain=summary_brain, summary_controls=summary_controls)
    
  },
  height = 800, 
  width = 1200
  )
  
  output$enrichments = renderPlot({
    
    test_list_for_enrichment(gene_list=input$gene_list_text)
    
  },
  height = 800, 
  width = 1200
  )
  
}

# initiate app
shinyApp(ui=ui, server=server)



