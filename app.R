# Phylostratigraphy App
VERSION = "0.2.2"

# Requisite libraries:
library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(queryup)
library(readxl)
library(org.Hs.eg.db)
library(phylostratr)
library(tidyverse)
library(topGO)
library(GO.db)
library(fgsea)
library(AnnotationDbi)
library(pbapply)
library(STRINGdb)
library(quarto)
library(DT)
library(tidyverse)
library(shinythemes)

#setwd("~/Documents/JainAnalytics/Kolabtree/Proposal1/")
#MasterGeneLists <- read_excel(path = "MasterGeneLists.xlsx")
#saveRDS(MasterGeneLists, file = "MasterGeneLists.rds")
MasterGeneLists <- readRDS("MasterGeneLists_Updated.rds")
# strata <- readRDS("9606_strata.rds")
# Merge results into a single hittable
#results <- merge_besthits(strata)
#ph <- stratify(results, classify_by_adjusted_pvalue(0.001))
#saveRDS(ph, file = "ph_0.001.rds")
ph <- readRDS("ph_0.001.rds")
levels(ph$mrca_name)[match("cellular organisms", levels(ph$mrca_name))] <-
  "Cellular Organisms"
Final_Taxa_Emergence_Timeline_with_Cellular_Organisms <- read.csv("Final_Taxa_Emergence_Timeline_with_Cellular_Organisms.csv")

levels(ph$mrca_name) <- paste0(levels(ph$mrca_name), ": ", paste0(Final_Taxa_Emergence_Timeline_with_Cellular_Organisms$Emergence..MYA., " MYA"))

#up_res <- UniProt.ws::queryUniProt(query = ph$qseqid, fields = c("id", "gene_names"))
#up <- UniProt.ws::UniProt.ws()

#pathways <- fgsea::gmtPathways("c5.go.v2023.2.Hs.symbols.gmt")
#saveRDS(pathways, file = "PhylostratigraphyApp/c5.go.v2023.2.Hs.symbols.rds")
pathways <- readRDS("c5.go.v2023.2.Hs.symbols.rds")
#gene_uniprot_conv_df <- UniPrgot.ws::select(up, keys = ph$qseqid, to = "Gene_Name")
#saveRDS(gene_uniprot_conv_df, "PhylostratigraphyApp/gene_uniprot_conv_df.rds")
#gene_uniprot_conv_df <- readRDS("PhylostratigraphyApp/uniprot_gene_conv_df.rds")
#gene_uniprot_conv_df <- read_delim("uniprotkb_organism_id_9606_AND_model_or_2024_05_12.tsv.gz", 
#                                                                     delim = "\t", escape_double = FALSE, 
#                                                                     trim_ws = TRUE)
#gene_uniprot_conv_df$Gene <- sapply(gene_uniprot_conv_df$`Gene Names`, function(x) strsplit(x, split = " ")[[1]][1])
#saveRDS(gene_uniprot_conv_df, file = "PhylostratigraphyApp/uniprot_gene_conv_df.rds")
gene_uniprot_conv_df <- readRDS("uniprot_gene_conv_df.rds")

ph$gene <- gene_uniprot_conv_df$Gene[match(ph$qseqid, gene_uniprot_conv_df$Entry)]
all_genes <- ph$gene[!is.na(ph$gene)]

ph_filt <- ph |>
  dplyr::filter(!duplicated(gene))


# chunks <- round(seq(from = 1, to = length(ph$qseqid), length.out = 100))

# for(i in 1:length(chunks)){
#   gene_uniprot_conv_df <- UniProt.ws::select(up, keys = ph$qseqid[(chunks[i]):(chunks[i+1]-1)], to = "Gene_Name")
#   levels(ph$mrca_name)[match("cellular organisms", levels(ph$mrca_name))] <- "Cellular Organisms"
#
# }
# doshiny()
# pboptions(
#   type = "shiny",
#   title = "Shiny progress",
#   label = "Computing functional enrichment for genes...")

##### UI #####
ui <- fluidPage(
  theme = shinytheme('flatly'),
  # Application title
  titlePanel(paste0("PhyloFastStrat v", VERSION)),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      shinyWidgets::switchInput(inputId = "runmode", label = "Run Mode", offLabel = "Preset Gene Lists", onLabel = "Custom Gene Lists", value = FALSE),
      conditionalPanel(
        condition = "!input.runmode",
        selectInput(
          inputId = "disease_of_interest",
          label = "Disease-Gene List:",
          choices = c(colnames(MasterGeneLists)),
          multiple = F,
          selected = NULL,
          selectize = T
        )
      ),
      conditionalPanel(
        condition = "input.runmode",
        textInput(inputId = "other_disease_of_interest", "Gene set name:"),
        selectizeInput(inputId = "selectized_genes", label = "Genes to query:",
                       choices = NULL,
                       multiple = TRUE,
                       selected = NULL,
                       options = list(
                         splitOn = I("(function() { return/[,; ]/; })()"),
                         plugins = list('remove_button'),
                         delimiter = " ",
                         persist = TRUE
                       )),
        fileInput(
          inputId = "custom_file",
          label = "Custom Gene List:",
          accept = ".csv",
          buttonLabel = "Upload"
        ),
        actionButton(inputId = "saverds", label = "Cache custom data"),
        textOutput("cache_message")
      ),
      actionButton("run", "Run analysis on selected genes"),
      htmlOutput("numGenesMappedUP"),
      headerPanel(""),
      shinyWidgets::switchInput(inputId = "no_isoforms", label = "Exclude alternative isoforms", onLabel = "Yes", offLabel = "No", value = FALSE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      type = "tabs",
      tabPanel("Linear Plots",
               plotOutput("ph_plot"),
               plotOutput("rate_plot")),
      tabPanel("Results table",
               DT::DTOutput("phylo_table")
               ),
      tabPanel("Gene map",
               selectInput("ps", label = "Select a phylostrata to visualize:", choices = unique(ph$mrca_name)),
               div(
                 style = "position: relative;",
                 withSpinner(plotOutput("string_plot", height = "800px"), type = 8)
               )
               ),
      tabPanel(
        "Functional enrichment plots",
        plotly::plotlyOutput("GOplot2", width = "100%", height = "700px") %>% withSpinner(),
        numericInput(inputId = "labeltop", label = "Number of top GO-BPs to label:", value = 5, min = 1, max = 10),
        numericInput(inputId = "cappval", label = "Max -log P-value to visualize", value = 100, min = 1, max = 100),
        plotOutput("GOplot_static", height = "700px")
        #plotly::plotlyOutput("duplicated_go_plot", width = "100%", height = "500px") %>% withSpinner()
      )
    )
    )
  )
)

##### SERVER #####
server <- function(input, output, session) {
  # Dynamic selectize on server-side
  updateSelectizeInput(session, inputId = "selectized_genes", choices = all_genes, server = TRUE, selected = NULL)

  # If run mode is switched, change disease of interest selection
  # observeEvent(input$runmode, {
  #   updateSelectInput(session, "disease_of_interest")
  #   #print("Triggered reset.")
  #   updateTextInput(session, "mytext", value = "test")
  # })
  
  # Set genes to query based on input, and refresh with "run" button
  ##### genes_of_interest ####
  genes_of_interest_reactive <- eventReactive(input$run, {
    if(!is.null(input$selectized_genes) & input$runmode){
      genes_of_interest <- input$selectized_genes
    }else if(!is.null(input$disease_of_interest) & !input$runmode){
      genes_of_interest <- MasterGeneLists[[input$disease_of_interest]]
    }else{
      custom_gene_list <- readr::read_csv(file = input$custom_file$datapath)
      genes_of_interest <- custom_gene_list[,1]
    }
    # Save a copy of genes_of_interest to compare input to what is processed
    genes_of_interest_input <- genes_of_interest[!is.na(genes_of_interest)]
    
    # Clean genes of interest
    genes_of_interest <-
      genes_of_interest[!is.na(genes_of_interest)]
    genes_of_interest <-
      genes_of_interest[!grepl("RNU", genes_of_interest)]
    genes_of_interest <- unique(genes_of_interest)
    genes_of_interest <-
      sapply(genes_of_interest, function(x)
        strsplit(x, split = " ", fixed = T)[[1]][1])
    
    # How many genes were found in the phylostrata object?
    genes_of_interest_found <- ph$gene[ph$gene %in% genes_of_interest]
    
    output$numGenesMappedUP <- renderUI({
      # Message to show gene mapping to phylostrat results
      HTML(paste0("Of ", length(genes_of_interest_input), " input genes, ", length(genes_of_interest), " gene IDs were cleaned and used.", "<br/>",
        "From ", length(genes_of_interest), " cleaned genes, ", length(genes_of_interest_found), " UniProt IDs were mapped (including isoforms)."))
    })
    
    genes_of_interest
  })
  
  ##### disease_name #####
  disease_name <- reactive({
    if(!is.null(input$other_disease_of_interest) & input$runmode){
      input$other_disease_of_interest
    }else{
      input$disease_of_interest
    }
  })
  
  observeEvent(input$saverds, {
    if(is.na(input$other_disease_of_interest)){
      cache_name <- input$custom_file$datapath
    }else{
      cache_name <- input$other_disease_of_interest
    }
    
    MasterGeneLists[[cache_name]] = c(genes_of_interest_reactive(), rep(NA, nrow(MasterGeneLists) - length(genes_of_interest_reactive())))
    
    saveRDS(MasterGeneLists, file = "MasterGeneLists.rds")
    output$cache_message <- renderText({"Successfully cached data!"})
  })
  
  ##### phylotable #####
  output$phylo_table <- DT::renderDT({
    #input$run
    if(input$no_isoforms){
      ph <- ph_filt
    }
    ph <- ph %>%
      filter(gene %in% genes_of_interest_reactive()) %>%
      #filter(mrca_name %in% input$ps) %>%
      group_by(Gene = gene, Phylostrata = mrca_name) %>%
      summarize(`Total Isoforms` = n())
    ph
  }, 
  extensions = 'Buttons', 
  server = F,
  options = list(
    paging = TRUE,
    searching = TRUE,
    #fixedColumns = TRUE,
    # scrollY = 400,
    autoWidth = TRUE,
    ordering = TRUE,
    lengthMenu = c(10, 25, 50, 100),
    dom = 'tB',
    buttons = c('copy', 'csv', 'excel', "pageLength")
  )
  )
  
  ##### plot of PS emergence ##### 
  output$ph_plot <- renderPlot({
    # Filter alt isoforms if specified
    if(input$no_isoforms){
      ph <- ph_filt
    }
    # Set up ph object
    genes_of_interest <- genes_of_interest_reactive()
    ph$is_of_interest <- ph$gene %in% genes_of_interest
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")
    ph$gene <-
      gene_uniprot_conv_df$Gene[match(ph$qseqid, gene_uniprot_conv_df$Entry)]

    # Reframe data
    ph_by_ps <- ph %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name = as.factor(mrca_name),
        annot_of_interest = as.factor(annot_of_interest)
      ) %>%
      dplyr::count(mrca_name = as.factor(mrca_name),
                   annot_of_interest,
                   .drop = F) %>%
      dplyr::mutate(gene_relevance = annot_of_interest)
    
    # Generate plot
    ggplot(
      data = ph_by_ps,
      aes(
        x = mrca_name,
        y = log10(n),
        group = gene_relevance,
        fill = gene_relevance,
        label = round(n, digits = 2)
      )
    ) +
      geom_line(aes(color = gene_relevance)) +
      geom_label(
        label.size = 0.05,
        alpha = 0.7,
        size = 3,
        vjust = -0.25,
        show.legend = F
      ) +
      labs(
        y = bquote(Number ~ of ~ novel ~ genes ~ (log[10] ~ scaled)),
        x = "Phylostrata",
        title = paste0(
          "Comparison of evolution rates of ",
          disease_name(),
          " genes to all other genes"
        ),
        fill = "",
        color = ""
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      ylim(c(0, 5)) +
      ggpubr::theme_pubr(x.text.angle = 90,
                         margin = T,
                         base_size = 11)
    
  })
  
  ##### string plot #####
  output$string_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    # Add error handling
    tryCatch({
      string_db <- STRINGdb::STRINGdb$new(species = 9606, version = "12", score_threshold = 200, input_directory = ".")
      ph$is_of_interest <- ph$gene %in% genes_of_interest_reactive()
      ph$annot_of_interest <-
        ifelse(ph$is_of_interest,
               disease_name(),
               "All Human Genes")
      
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        paste0(hcl(h = hues, l = 65, c = 100)[1:n], "FF")
      }
      
      withProgress(message = "Pulling data from STRINGdb...", value = 0, {
        ph_mapped <- string_db$map(ph %>% filter(is_of_interest & mrca_name %in% input$ps) %>% as.data.frame(), "gene", removeUnmappedRows = T)
        incProgress(amount = 0.5, message = "Generating plot...")
        string_db$plot_network(ph_mapped$STRING_id, required_score = 800)
        incProgress(amount = 0.5, message = "Done!")
      })
    }, error = function(e) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, paste("Error generating STRING plot:", e$message), cex = 1.2)
    })
  })
  
  ##### normalized evo rate  #####
  output$rate_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    ph$is_of_interest <- ph$gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")

    ph_by_ps <- ph %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name = as.factor(mrca_name),
        annot_of_interest = as.factor(annot_of_interest)
      ) %>%
      dplyr::count(mrca_name = as.factor(mrca_name),
                   annot_of_interest,
                   .drop = F) %>%
      dplyr::mutate(gene_relevance = annot_of_interest)
    
    
    fixed_ratio <- ph_by_ps %>%
      group_by(gene_relevance) %>%
      summarize(total_in_group = sum(n)) %>%
      pull(total_in_group)
    fixed_ratio <- min(fixed_ratio) / sum(fixed_ratio)
    
    ggplot(
      data = ph_by_ps %>%
        pivot_wider(
          id_cols = c(mrca_name),
          names_from = gene_relevance,
          values_from = n
        ) %>%
        mutate(Ratio = (.[[3]] / .[[2]]) / fixed_ratio) %>%
        filter(!is.na(Ratio)),
      aes(x = mrca_name, y = Ratio)
    ) +
      geom_hline(yintercept = 1, color = "grey80") +
      geom_point() +
      ggpubr::theme_pubr(x.text.angle = 90) +
      geom_line(group = "a") +
      labs(
        x = "Phylostrata",
        y = paste0(
          "Ratio of Evolutionary Rate \n(Disease Genes : Total Human Genome)"
        ),
        title = paste0(
          "Normalized ratio of novel gene emergence for ",
          disease_name(),
          " to overall human genes"
        )
      )
  })
  
  ##### GO reactive data ####
  GO_df_reactive <- reactive({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    ph$is_of_interest <- ph$gene %in% genes_of_interest_reactive()
    ph$annot_of_interest <-
      ifelse(ph$is_of_interest,
             disease_name(),
             "All Human Genes")

    ph_of_interest <- ph %>%
      filter(is_of_interest == TRUE)
    
    gene_list_by_strata <-
      lapply(unique(ph_of_interest$ps), function(x) {
        ph_of_interest$gene[ph_of_interest$ps == x]
      })
    names(gene_list_by_strata) <- unique(ph_of_interest$ps)
    
    # This takes a while, but these are the best results...
    annFUN.org.symbol <-
      topGO::annFUN.org(whichOnto = "BP",
                        mapping = "org.Hs.eg.db",
                        ID = "symbol")
    
    GO_df_by_strata <-
      lapply(names(gene_list_by_strata), function(i) {
        gene_set_to_enrich <- unique(gene_list_by_strata[[i]])
        gene_set_to_enrich <-
          gene_set_to_enrich[!is.na(gene_set_to_enrich)]
        if (length(gene_set_to_enrich) < 1) {
          return(NA)
        }
        
        allGenes <- unique(unlist(annFUN.org.symbol))
        geneList <-
          as.factor(as.integer(allGenes %in% gene_set_to_enrich))
        names(geneList) <- allGenes
        geneList_numeric <- as.numeric(geneList[geneList == 1])
        names(geneList_numeric) <- names(geneList)[geneList == 1]
        
        pathways <- annFUN.org.symbol
        
        allres <- fgsea(pathways = pathways,
                        stats = geneList_numeric,
                        minSize = 1)
        
        # OLD WAY USING TOPGO - 
        # GOdata <- new(
        #   "topGOdata",
        #   ontology = "BP",
        #   allGenes = geneList,
        #   annot = annFUN.org,
        #   mapping = "org.Hs.eg.db",
        #   ID = "symbol"
        # )
        # GOgraph <- graph(GOdata)
        # gs <- geneScore(GOdata, whichGenes = gene_set_to_enrich)
        
        # resultFisher <-
        #   runTest(GOdata, algorithm = "classic", statistic = "fisher")
        # resultKS <- runTest(GOdata, algorithm = "elim", statistic = "ks")
        
        #allres <- GenTable(GOdata, classic = resultFisher, numChar = 1000)
        ## Add progress:
        #incProgress(1/length(gene_list_by_strata), detail = paste0("Computing enrichment for ", i))
        
        return(allres)
        
        # return(
        #   res %>%
        #     dplyr::mutate(
        #       Pathway_Name = gsub(
        #         pattern = "GO(BP|MF|CC)",
        #         x = gsub("_", x = pathway, replacement = " "),
        #         replacement = ""
        #       )
        #     ) %>%
        #     dplyr::mutate(Pathway_Name = str_to_title(Pathway_Name)) %>%
        #     dplyr::slice_max(abs(NES), n = 30)
        # )
      }) 
    
    
    names(GO_df_by_strata) <-
      ph$mrca_name[match(names(gene_list_by_strata), ph$ps)]
    
    GO_res_df <-
      bind_rows(lapply(GO_df_by_strata, as.data.frame), .id = "PS") %>%
      mutate(PS = factor(PS, levels = levels(ph$mrca_name))) %>%
      mutate(Term = AnnotationDbi::Term(pathway)) %>%
      filter(!is.na(pval))
    
    GO_res_df
  })
  
  ##### Dup GO term plot #####
  #output$duplicated_go_plot <- plotly::renderPlotly({
    # ph$is_of_interest <- ph$qseqid %in% up_of_interest_reactive()
    # ph$annot_of_interest <-
    #   ifelse(ph$is_of_interest,
    #          input$disease_of_interest,
    #          "All Human Genes")
    # ph$gene <-
    #   gene_uniprot_conv_df$To[match(ph$qseqid, gene_uniprot_conv_df$From)]
    # 
    # ph_of_interest <- ph %>%
    #   filter(is_of_interest == TRUE)
    # 
    # gene_list_by_strata <-
    #   lapply(unique(ph_of_interest$ps), function(x) {
    #     ph_of_interest$gene[ph_of_interest$ps == x]
    #   })
    # names(gene_list_by_strata) <- unique(ph_of_interest$ps)
    # 
    # annFUN.org.symbol <-
    #   topGO::annFUN.org(whichOnto = "BP",
    #                     mapping = "org.Hs.eg.db",
    #                     ID = "symbol")
    # 
    # GO_df_by_strata <-
    #   lapply(names(gene_list_by_strata), function(i) {
    #     gene_set_to_enrich <- unique(gene_list_by_strata[[i]])
    #     gene_set_to_enrich <-
    #       gene_set_to_enrich[!is.na(gene_set_to_enrich)]
    #     if (length(gene_set_to_enrich) < 2) {
    #       return(NA)
    #     }
    #     
    #     allGenes <- unique(unlist(annFUN.org.symbol))
    #     geneList <-
    #       as.factor(as.integer(allGenes %in% gene_set_to_enrich))
    #     names(geneList) <- allGenes
    #     geneList_numeric <- as.numeric(geneList[geneList == 1])
    #     names(geneList_numeric) <- names(geneList)[geneList == 1]
    #     
    #     res <- fgsea(pathways = pathways,
    #                  stats = geneList_numeric,
    #                  minSize = 0)
    #     
    #     # GOdata <- new(
    #     #   "topGOdata",
    #     #   ontology = "BP",
    #     #   allGenes = geneList,
    #     #   annot = annFUN.org,
    #     #   mapping = "org.Hs.eg.db",
    #     #   ID = "symbol"
    #     # )
    #     # GOgraph <- graph(GOdata)
    #     # gs <- geneScore(GOdata, whichGenes = gene_set_to_enrich)
    #     
    #     # resultFisher <-
    #     #   runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #     # resultKS <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    #     
    #     # allres <- GenTable(GOdata, classic = resultFisher, numChar = 1000)
    #     
    #     # return(allres)
    #     
    #     return(
    #       res %>%
    #         dplyr::mutate(
    #           Pathway_Name = gsub(
    #             pattern = "GO(BP|MF|CC)",
    #             x = gsub("_", x = pathway, replacement = " "),
    #             replacement = ""
    #           )
    #         ) %>%
    #         dplyr::mutate(Pathway_Name = str_to_title(Pathway_Name)) %>%
    #         dplyr::slice_max(abs(NES), n = 30)
    #     )
    #   })
    
  #   GO_res_df <- GO_df_reactive()
  #   
  #   if (nrow(GO_res_df) == 0) {
  #     plot_out <- ggplot()
  #   } else{
  #     GO_res_df_dup <- GO_res_df %>%
  #       dplyr::filter(base::duplicated(Term) |
  #                       base::duplicated(Term, fromLast = T)) %>%
  #       dplyr::filter(!is.na(Term)) %>%
  #       dplyr::filter(pval < 0.1)
  #     
  #     # IF too many terms
  #     if(nrow(GO_res_df_dup) > 20){
  #       # Get the terms that appear the most number of times
  #       term_counts <- table(GO_res_df_dup$Term)
  #       top20_terms <- names(term_counts)[order(term_counts, decreasing = T)[1:20]]
  #       GO_res_df_dup <- GO_res_df_dup %>% filter(Term %in% top20_terms)
  #     }
  #     
  #     plot_out <- ggplot(GO_res_df_dup,
  #                        aes(x = PS, y = Term, group = Term)) +
  #       geom_line() +
  #       geom_point(aes(color = -log10(as.numeric(pval)))) +
  #       #scale_x_continuous(n.breaks = 27) +
  #       ggpubr::theme_pubr(x.text.angle = 90, base_size = 8, legend = "none") +
  #       labs(
  #         x = "Phylostrata",
  #         y = "GO Term",
  #         title = paste0(
  #           "BPs with multiple evolutionary origins related to ",
  #           disease_name()
  #         ),
  #         color = "Hypergeometric Test -log10 P-value"
  #       )
  #     
  #     if (nrow(GO_res_df_dup) == 0) {
  #       plot_out <- plot_out +
  #         annotate(
  #           geom = "text",
  #           x = 1,
  #           y = 1,
  #           label = "No repeated GO terms over phylostrata."
  #         )
  #     }
  #   }
  #   
  #   plotly::ggplotly(p = plot_out)
  # })
  
  ##### Plotly GO #####
  output$GOplot2 <- plotly::renderPlotly({
    GO_res_df <- GO_df_reactive()
    
    GO_res_df_top <- GO_res_df %>%
      group_by(PS) %>%
      arrange(desc(padj)) %>%
      #slice_min(n = input$maxGoTerms, order_by = classic) %>%
      ungroup()
    
    plot2 <- ggplot(GO_res_df_top, aes(
      x = PS,
      y = -log10(as.numeric(pval)),
      label = `Term`,
      color = PS
    )) +
      # geom_hline(
      #   yintercept = -log10(0.05),
      #   linetype = 2,
      #   color = "grey70"
      # ) +
      geom_point(show.legend = F, position = position_jitter(0.1)) +
      geom_label(
        data = GO_res_df_top,
        aes(label = `Term`, y = -log(as.numeric(pval))),
        show.legend = F,
        size = 3,
        hjust = 0.5,
        vjust = 1
      ) +
      scale_x_discrete(limits = levels(GO_res_df$PS)) +
      ggpubr::theme_pubr(x.text.angle = 90, legend = "none") +
      labs(
        x = "Phylostrata (PS)",
        y = "-log P-value (Fisher test)",
        title = "Top enriched GO terms by phylostrata",
        subtitle = "Top 3 in each PS are labeled"
      )

    plotly::ggplotly(plot2, tooltip = "label")
  })
  
  ##### Static GO ####
  output$GOplot_static <- renderPlot({
    GO_res_df <- GO_df_reactive()
    
    if(input$cappval < 100){
      GO_res_df$pval[-log10(GO_res_df$pval) > input$cappval] <- 10^(-1*input$cappval)
    }
    
    GO_res_df_top <- GO_res_df %>%
      group_by(PS) %>%
      arrange(desc(pval)) %>%
      slice_min(n = input$labeltop, order_by = pval, with_ties = F) %>%
      ungroup()
    
    ggplot(GO_res_df, aes(
      x = PS,
      y = -log10(as.numeric(pval)),
      label = `Term`,
      color = PS
    )) +
      # geom_hline(
      #   yintercept = -log10(0.05),
      #   linetype = 2,
      #   color = "grey70"
      # ) +
      geom_point(show.legend = F, position = position_jitter(0.1)) +
      ggrepel::geom_label_repel(
        data = GO_res_df_top,
        aes(label = `Term`, y = -log10(as.numeric(pval))),
        show.legend = F,
        size = 3,
        hjust = 0.5,
        vjust = 1, max.overlaps = 30
      ) +
      scale_x_discrete(limits = levels(GO_res_df$PS)) +
      ggpubr::theme_pubr(x.text.angle = 90, legend = "none") +
      labs(
        x = "Phylostrata (PS)",
        y = "-log P-value (Fisher test)",
        title = "Top enriched GO terms by phylostrata",
        subtitle = "Top 5 in each PS are labeled"
      )
    #plotly::ggplotly(plot2, tooltip = "label")
  })
  
}

##### RUN #####
shinyApp(ui = ui, server = server)
