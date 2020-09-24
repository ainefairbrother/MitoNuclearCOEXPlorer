# shiny interface to allow querying of correlation data

# ---------------------- imports ---------------------------------------

# allows bioconductor pkgs to be found by shiny.io server

setwd("/home/abrowne/shiny/test_shiny/")
#load("./app_data.Rdata")

options(warn=-1)
library(BiocManager)
library(tidyverse)
library(rtracklayer)
library(fst)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(plyr)
library(dplyr)
library(parallel)
library(ggpubr)
options(gsubfn.engine = "R")

# importing fns
source("./R/convert_sym_ens.R")
source("./R/gen_gene_specific_distr_plot_fn.R")
source("./R/single_gene_gen_fig.R")
source("./R/test_list_for_enrichment_fn.R")


# ---------------------- shiny initialise ------------------------------
myapp = shinyApp(
  
  # user interface
  ui = navbarPage("Mito-nuclear brain browser",
                  
                  #theme = shinytheme("yeti"),
                  
                  # tags$style(HTML(".navbar-default .navbar-nav > li > a[data-value='Explore data with gene'] {color: #7a89b8;}")),
                  # tags$style(HTML(".navbar-default .navbar-nav > li > a[data-value='Explore data with gene set'] {color: #7a89b8;}")),
                  
                  tabPanel("Welcome",
                           
                           fluidPage(
                             
                             br(),
                             
                             titlePanel(
                               h2(strong("Mito-nuclear brain browser"), align = "center")
                             ),
                             
                             #div(img(src="mt_nuc_schemat.png", width="10%"), style="text-align: center;"),
                             
                             br(),
                             
                             fluidRow(
                               column(8,
                                      offset=2,
                                      wellPanel(
                                        h4(strong("Welcome to the mito-nuclear brain browser. This tool allows exploration of relationships between target gene(s) of
                                                  nuclear origin and the mitochondrial genome")),
                                        br(),
                                        p("Levaraging transcriptomic data from the GTEx project, nuclear-mitochondrial
                                          co-ordination was quantified by correlating expression levels of nDNA- and mtDNA-encoded genes."),
                                        br(),
                                        div(img(src="corr_mat_schematic.png", width="60%"), style="text-align: center;"),
                                        br(),
                                        p("The rationale here being that co-expressed genes may be directly regulating each other, controlled by the same regulatory regime,
                                          functionally related, members of the same pathway or part of the same protein complex, and as such,
                                          the correlation of their expression is meaningful."),
                                        p(strong("See below for analyses options offered by the mito-nuclear brain browser tool and example outputs."))
                                      )
                               ),
                               
                               column(5,
                                      offset=1,
                                      wellPanel(
                                        h4(strong('Explore data with gene')),
                                        hr(),
                                        p("At the resolution of a single nuclear gene,",em("x,")," the mito-nuclear brain browser can run analyses to determine:"),
                                        p(strong("1."), "The Spearman correlation value of,",em("x,")," with each mtDNA-encoded protein coding gene in 12
                                          CNS tissues and two other body tissues. This is output in the form of a heatmap. This example uses the ATG7 gene:"),
                                        br(),
                                        div(img(src="example_heatmap.png", width="100%"), style="text-align: center;"),
                                        br(),
                                        br(),
                                        p(strong("2."), "A visualisation of the position of gene, ",em("x"), ", the overall nuclear-mitochondrial correlation distribution for 12 GTEx CNS regions.
                                          This is output in the form of a 12-panel density plot. This example uses the ATG7 gene:"),
                                        br(),
                                        div(img(src="example_dist.png", width="100%"), style="text-align: center;"),
                                        br()
                                      )
                               ),
                               
                               column(5,
                                      wellPanel(
                                        h4(strong('Explore data with gene set')),
                                        hr(),
                                        p("The gene set exploration function takes a target gene set,",em("t,")," (n>1) and tests whether it has a non-random
                                          association with the mitochondrial genome."),
                                        p("This is acheived by calculating whether the median correlation value of,",em("t,")," is more extreme than the medians
                                          of randomly selected sets. The tool reports the following statistics on the top right of the density plot output:"),
                                        p(strong("i."),"The median mito-nuclear correlation value of all the random set medians (green dashed line),"),
                                        p(strong("ii."),"The median mito-nuclear correlation value of the,",em("t,")," (blue dashed line),"),
                                        p(strong("iii."),"The P-value for the shift of ",em("t "),"from random gene sets."),
                                        p("This example uses a Parkinson's Disease implicated gene set, with the background filtered for protein coding genes and
                                          10,000 bootstrap lists:"),
                                        br(),
                                        div(img(src="example_target_rand_list.png", width="100%"), style="text-align: center;"),
                                        br(),
                                        h4('Resources'),
                                        p("Gene sets can be downloaded from:"),
                                        tags$ul(
                                          tags$li(tagList(a("MSigDB", href="https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp")), " for gene sets relating to molecular pathways"),
                                          tags$li(tagList(a("Genomics England PanelApp", href="https://panelapp.genomicsengland.co.uk/panels/")), " for clinically curated gene sets relating to human disorders")
                                        )
                                      ))))),
                  
                  tabPanel(strong("Explore data with gene"),
                           
                           sidebarPanel(
                             offset=3,
                             h4(strong('Explore data with gene')),
                             hr(),
                             h4('Gene ID'),
                             helpText("Enter an Ensembl gene code or gene symbol", style = "font-size:11px;"),
                             textInput(inputId = "gene_text", label=NULL, value = ""),
                             submitButton(text = "Submit", icon = NULL, width = NULL),
                             hr(),
                             h4('Download'),
                             downloadButton("downloadSINGLEGENE", "Download plot"),
                             width=3
                           ),
                           
                           mainPanel(
                             width=4,
                             fluidRow(
                               plotOutput(outputId = "single_gene_plots") %>% withSpinner(color="#0dc5c1")
                             )
                           )
                  ),
                  
                  tabPanel(strong("Explore data with gene set"),
                           
                           sidebarPanel(
                             
                             h4(strong("Explore data with gene set")),
                             hr(),
                             h4("Target gene set"),
                             helpText("Enter a list of comma and space separated ENS IDs or gene symbols (reccomended min. n=20)", style = "font-size:11px;"),
                             helpText("e.g. ENSG00000101986, ENSG00000141385, ENSG00000003393, ENSG00000214274,...", style = "font-size:11px;"),
                             # input box to insert gene list
                             textInput(inputId = "gene_list_text", label=NULL, value="", width="100%"),
                             h4("Background gene set"),
                             helpText("If your gene list is protein coding, it is preferable to match its biotype to that of the background list", style = "font-size:11px;"),
                             # checkbox to filter background list
                             checkboxInput(inputId="bkg", label="Filter background for protein coding genes only", value=FALSE), 
                             h4("Bootstraps"),
                             helpText("Enter the number of random gene sets to compare the target set against.
                                      1000 recommended for accuracy of P-values, but is slower to render
                                      ", style = "font-size:11px;"),
                             radioButtons(inputId="iters", label=NULL, choices = c("10"=10, "100"=100, "1000"=1000, "10000 (~15mins)"=10000), selected = "10", inline = FALSE, width = NULL),
                             br(),
                             submitButton(text = "Submit", icon = NULL, width = NULL),
                             hr(),
                             h4('Download'),
                             downloadButton("downloadENRICHMENT", "Download plot"),
                             width=3),
                           
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Analyse",
                                        plotOutput("enrichments") %>% withSpinner(color="#0dc5c1")
                               ),
                               tabPanel("Help", 
                                        fluidRow(
                                          br(),
                                          wellPanel(
                                            h4("Target set details"),
                                            htmlOutput("helpmsg") %>% withSpinner(color="#0dc5c1"),
                                            br(),
                                            htmlOutput("helpmsg1"),
                                            tags$head(tags$style("#helpmsg1{color: red;}"))
                                          )
                                        ),
                                        fluidRow(
                                          wellPanel(
                                            h4("Gene set analysis: method details"),
                                            p("For each GTEx CNS region, r, and the input gene set, l, the median nuclear-mitochondrial correlation value of l for r was calculated. 
                                              The distribution of nuclear-mitochondrial pairs is inclusive of all mitochondrial correlations for each gene, making the size of the distribution 
                                              (length of l)*13. To generate empirical distributions, a random sample of nuclear genes of length l is selected from the list of all genes expressed 
                                              in all GTEx CNS regions (15,001) and all correlations with mitochondrial genes are included.
                                              A two-tailed test is then carried out to determine whether l has a more extreme median nuclear-mitochondrial correlation value than could be expected by chance. 
                                              To do this, the median of l is compared against k (the number of bootstrap sets used) medians of randomly selected gene sets. 
                                              P-values are calculated as follows, where n is the number of random medians less extreme than the median of l:"),
                                            p("P = (k-n)/k"),
                                            p("Calculation of a 0 P-value is a limitation imposed by the number of bootstraps, and so in these cases a P-value of P<1/(number of bootstraps) is reported.")
                                          )))),
                             width=6)),
                  
                  tabPanel('Download',
                           
                           fluidPage(
                             fluidRow(
                               column(8,
                                      offset=2,
                                      wellPanel(
                                        h4('Download the per-region Spearman correlation values and corresponding p-values for all mito-nuclear gene pairs'),
                                        br(),
                                        br(),
                                        p('Download correlations and p-values for 12 CNS regions (97.7 MB)'),
                                        downloadButton("downloadCNS", "Download"),
                                        br(),
                                        br(),
                                        p('Download correlations and p-values for heart and muscle (23 MB)'),
                                        downloadButton("downloadCTRL", "Download")
                                      ))))),
                  
                  tabPanel('Cite us',
                           
                           fluidPage(
                             fluidRow(
                               column(8,
                                      offset=2,
                                      wellPanel(
                                        p('[Link to pre-print here]')
                                      ))))),
                  
                  tabPanel("Contact",
                           
                           fluidPage(
                             fluidRow(
                               column(8,
                                      offset=2,
                                      wellPanel(
                                        h3("Mito-nuclear brain browser was developed by Aine Fairbrother-Browne"),
                                        p('A collaboration between the Hodgkinson and Ryten laboratories'),
                                        br(),
                                        h4('For any questions related to this resource or publication please contact:'),
                                        p('Aine Fairbrother-Browne for technical issues and general questions about the project - aine.fairbrother-browne.18@ucl.ac.uk'),
                                        p('Sonia García-Ruiz for technical issues relating to the user interface (UI) - s.ruiz@ucl.ac.uk'),
                                        br(),
                                        h4("Ryten Lab (University College London):"),
                                        p("UCL Great Ormond Street Institute of Child Health"),
                                        p("30 Guilford Street"),
                                        p("London WC1N 1EH"),
                                        p(tagList(a("Visit the Ryten Lab", href="https://snca.atica.um.es/"))),
                                        br(),
                                        h4("Hodgkinson Lab (King's College London):"),
                                        p("Guy's Hospital"),
                                        p("Great Maze Pond"),
                                        p("London SE1 9RT"),
                                        p(tagList(a("Visit the Hodgkinson Lab", href="https://www.hodgkinsonlab.org/")))
                                        
                                      )))))
  ),
  
  # backend
  server = function(input, output){
    
    output$single_gene_plots = renderPlot({
      
      # check to see if input gene is present in input data
      gene_sym_check = summary_brain %>%
        dplyr::filter(input$gene_text %in% nuc_gene)
      gene_ens_check = summary_brain %>%
        dplyr::filter(input$gene_text %in% gene_name)
      
      validate(need(nrow(gene_sym_check) > 0 | nrow(gene_ens_check) > 0, "Gene not found in data"))
      validate(need(input$gene_text != "", "Please provide a gene ID (ENS/symbol)"))
      
      GenFigSingleGene(input$gene_text)
      
    },
    height = 1000,
    width = 800
    )
    
    output$helpmsg <- renderText({
      
      gene_list = as.character(unique(strsplit(input$gene_list_text, ", ")[[1]]))
      genes_in_data = summary_brain %>% dplyr::filter(nuc_gene %in% gene_list) %>% dplyr::select(nuc_gene)
      genes_in_data = genes_in_data$nuc_gene
      genes_in_data = as.character(unique(genes_in_data))
      genes_not_found = as.character(setdiff(gene_list, genes_in_data))
      
      grch_list = grch %>%
        dplyr::filter(gene_id %in% genes_in_data) %>% 
        dplyr::select(c(gene_id, gene_biotype))
      
      msg = HTML(paste(
        paste0(strong("Size of target set: "), length(gene_list)),
        paste0(strong("Number of target set genes found in database: "), paste0(length(unique(genes_in_data)), "/", length(unique(gene_list)))),
        paste0(strong("Target set genes not found in database: "), genes_not_found),
        paste0(strong("Biotypes present in target set: "), unique(grch_list$gene_biotype)),
        sep="<br/>"))
      
      return(msg)
      
    })
    
    output$helpmsg1 <- renderText({
      
      gene_list = as.character(unique(strsplit(input$gene_list_text, ", ")[[1]]))
      genes_in_data = summary_brain %>% dplyr::filter(nuc_gene %in% gene_list) %>% dplyr::select(nuc_gene)
      genes_in_data = genes_in_data$nuc_gene
      genes_in_data = as.character(unique(genes_in_data))
      
      if((length(genes_in_data) < 20) & (length(gene_list) > 0)){
        msg1 = HTML("Warning: gene set less than reccomended set size")
        return(msg1)
      }else(
        return("")
      )
    })
    
    output$enrichments = renderPlot({
      
      gene_list = as.character(unique(strsplit(input$gene_list_text, ", ")[[1]]))
      genes_in_data = summary_brain %>% dplyr::filter(nuc_gene %in% gene_list) %>% dplyr::select(nuc_gene)
      genes_in_data = genes_in_data$nuc_gene
      genes_in_data = as.character(unique(genes_in_data))
      
      validate(need(input$gene_list_text != "", "
                    Please provide a list of gene IDs (ENS/symbol) in the form: gene1, gene2, gene3"))
      validate(need(length(genes_in_data) != 0, "
                    Please provide a list of gene IDs (ENS/symbol) in the form: gene1, gene2, gene3"))
      validate(need(length(gene_list) >= 2, "
                    Less than 2 genes entered"))
      validate(need(length(genes_in_data) >= 2, "
                    Less than 2 genes found in data"))
      
      test_list_for_enrichment(gene_list=input$gene_list_text, iters=input$iters, filt_bkg_for_protein_coding=input$bkg)
      
    },
    height = 950, 
    width = 1100
    )
    
    # downloads
    output$downloadCNS <- downloadHandler(
      filename = function(){
        paste("correlations_12_CNS_regions","csv",sep=".")
      },
      content = function(con){
        file.copy('./data/GTEx_brain_summary_table.csv', con)
      })
    
    output$downloadCTRL <- downloadHandler(
      filename = function(){
        paste("correlations_heart_and_muscle","csv",sep=".")
      },
      content = function(con){
        file.copy('./data/GTEx_control_tissues_summary_table.csv', con)
      })
    
    output$downloadSINGLEGENE <- downloadHandler(
      filename = function(){
        paste0(input$gene_text, "_", Sys.time(), "_mito_nuclear_bb.png", sep="")
      },
      content=function(con){
        ggsave(con, device = "png", width=30, height=35, units="cm")
      }
    )
    
    output$downloadENRICHMENT <- downloadHandler(
      filename = function(){
        paste0("mito_nuc_BB_plot_n=", length(strsplit(input$gene_list_text, ', ')[[1]]), "_genes_bootstraps=", input$iters, "_", Sys.time(), ".png", sep="")
      },
      content=function(file){
        ggsave(file, device = "png", width=40, height=40, units="cm")
      }
    )
  }
  
)

# initiate app
#runApp(myapp)



