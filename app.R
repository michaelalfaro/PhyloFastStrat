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
library(digest)  # For creating cache keys
library(GGally)  # For correlation matrix plot
library(gridExtra)  # For arranging plots
library(grDevices)  # For PDF creation

# Create cache directory if it doesn't exist
if (!dir.exists("cache")) {
  dir.create("cache")
}

# Function to get cached STRINGdb results
get_cached_string_results <- function(genes, ps) {
  # Create a unique cache key based on genes and phylostrata
  cache_key <- digest::digest(paste(sort(genes), ps, collapse = "_"))
  cache_file <- file.path("cache", paste0(cache_key, ".rds"))
  
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  return(NULL)
}

# Function to save STRINGdb results to cache
save_to_cache <- function(genes, ps, results) {
  cache_key <- digest::digest(paste(sort(genes), ps, collapse = "_"))
  cache_file <- file.path("cache", paste0(cache_key, ".rds"))
  saveRDS(results, cache_file)
}

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
          label = "Disease-Gene Lists:",
          choices = c(colnames(MasterGeneLists)),
          multiple = TRUE,
          selected = NULL,
          selectize = T
        )
      ),
      conditionalPanel(
        condition = "input.runmode",
        textInput(inputId = "other_disease_of_interest", "Gene set name:", placeholder = "Enter a name for your gene set"),
        div(
          style = "margin-bottom: 15px;",
          textAreaInput(
            inputId = "pasted_genes",
            label = "Paste genes (one per line or comma-separated):",
            placeholder = "Paste your genes here...",
            height = "150px",
            resize = "vertical"
          )
        ),
        selectizeInput(inputId = "selectized_genes", label = "Or select genes from list:",
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
          label = "Or upload a CSV file:",
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
               div(
                 style = "margin-bottom: 10px;",
                 downloadButton("download_plots", "Save All Plots as PDF", 
                              style = "margin-bottom: 10px;")
               ),
               plotOutput("ph_plot"),
               plotOutput("rate_plot"),
               plotOutput("correlation_matrix")),
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
      # Handle multiple disease sets
      genes_of_interest <- lapply(input$disease_of_interest, function(x) {
        genes <- MasterGeneLists[[x]]
        genes[!is.na(genes)]
      })
      names(genes_of_interest) <- input$disease_of_interest
    }else if(!is.null(input$pasted_genes) && input$pasted_genes != ""){
      # Handle pasted genes
      genes_of_interest <- unlist(strsplit(input$pasted_genes, split = "[\n,]+"))
      genes_of_interest <- trimws(genes_of_interest)
      genes_of_interest <- genes_of_interest[genes_of_interest != ""]
      genes_of_interest <- list(genes_of_interest)
      names(genes_of_interest) <- ifelse(
        input$other_disease_of_interest != "",
        input$other_disease_of_interest,
        paste0("Custom Gene Set ", format(Sys.time(), "%Y%m%d_%H%M%S"))
      )
    }else{
      custom_gene_list <- readr::read_csv(file = input$custom_file$datapath)
      genes_of_interest <- list(custom_gene_list[,1])
      names(genes_of_interest) <- ifelse(
        input$other_disease_of_interest != "",
        input$other_disease_of_interest,
        paste0("Custom Gene Set ", format(Sys.time(), "%Y%m%d_%H%M%S"))
      )
    }
    
    # Set default name if none provided
    if(input$runmode && (is.null(input$other_disease_of_interest) || input$other_disease_of_interest == "")) {
      updateTextInput(session, "other_disease_of_interest", value = paste0("Custom Gene Set ", format(Sys.time(), "%Y%m%d_%H%M%S")))
    }
    
    # Clean genes of interest
    genes_of_interest <- lapply(genes_of_interest, function(genes) {
      genes <- genes[!is.na(genes)]
      genes <- genes[!grepl("RNU", genes)]
      genes <- unique(genes)
      genes <- sapply(genes, function(x) strsplit(x, split = " ", fixed = T)[[1]][1])
      genes
    })
    
    # How many genes were found in the phylostrata object?
    genes_of_interest_found <- lapply(genes_of_interest, function(genes) {
      ph$gene[ph$gene %in% genes]
    })
    
    # Calculate total genes for display
    total_genes <- sum(sapply(genes_of_interest, length))
    total_mapped <- sum(sapply(genes_of_interest_found, length))
    
    output$numGenesMappedUP <- renderUI({
      HTML(paste0("Of ", total_genes, " input genes, ", total_genes, " gene IDs were cleaned and used.", "<br/>",
        "From ", total_genes, " cleaned genes, ", total_mapped, " UniProt IDs were mapped (including isoforms)."))
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
  
  # Add reactive values to store plots
  plot_values <- reactiveValues(
    ph_plot = NULL,
    rate_plot = NULL,
    correlation_matrix = NULL
  )
  
  ##### plot of PS emergence ##### 
  output$ph_plot <- renderPlot({
    # Filter alt isoforms if specified
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    # Set up ph object for multiple gene sets
    genes_of_interest <- genes_of_interest_reactive()
    ph$is_of_interest <- FALSE
    ph$annot_of_interest <- "All Human Genes"
    
    # Add each gene set
    for(set_name in names(genes_of_interest)) {
      ph$is_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- TRUE
      ph$annot_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- set_name
    }

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
    p <- ggplot(
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
        title = "Comparison of evolution rates across gene sets",
        fill = "",
        color = ""
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      ylim(c(0, 5)) +
      ggpubr::theme_pubr(x.text.angle = 90,
                         margin = T,
                         base_size = 11)
    
    # Store the plot
    plot_values$ph_plot <- p
    p
  })
  
  ##### normalized evo rate  #####
  output$rate_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    # Set up ph object for multiple gene sets
    genes_of_interest <- genes_of_interest_reactive()
    ph$is_of_interest <- FALSE
    ph$annot_of_interest <- "All Human Genes"
    
    # Add each gene set
    for(set_name in names(genes_of_interest)) {
      ph$is_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- TRUE
      ph$annot_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- set_name
    }

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
    
    # Calculate ratios for each gene set independently
    ratio_data <- ph_by_ps %>%
      group_by(gene_relevance) %>%
      mutate(total_in_group = sum(n)) %>%
      ungroup() %>%
      group_by(mrca_name) %>%
      mutate(
        total_genes = sum(n),
        ratio = n / total_genes
      ) %>%
      mutate(
        normalized_ratio = ifelse(gene_relevance == "All Human Genes", 
                                1,  # Baseline is always 1
                                ratio / (total_in_group / sum(total_in_group)))
      )
    
    # Create separate data frames for baseline and gene sets
    baseline_data <- ratio_data %>% filter(gene_relevance == "All Human Genes")
    gene_sets_data <- ratio_data %>% filter(gene_relevance != "All Human Genes")
    
    p <- ggplot() +
      # Plot baseline
      geom_line(data = baseline_data, 
                aes(x = mrca_name, y = normalized_ratio, group = gene_relevance),
                color = "grey50", linetype = "dashed", size = 1) +
      # Plot gene sets
      geom_line(data = gene_sets_data,
                aes(x = mrca_name, y = normalized_ratio, 
                    color = gene_relevance, group = gene_relevance),
                size = 1) +
      geom_point(data = gene_sets_data,
                 aes(x = mrca_name, y = normalized_ratio, 
                     color = gene_relevance, group = gene_relevance)) +
      geom_hline(yintercept = 1, color = "grey80") +
      ggpubr::theme_pubr(x.text.angle = 90) +
      labs(
        x = "Phylostrata",
        y = "Normalized Ratio of Gene Emergence",
        title = "Normalized ratio of novel gene emergence across gene sets",
        color = "Gene Set"
      ) +
      scale_color_brewer(palette = "Set1") +
      annotate("text", x = Inf, y = Inf, 
               label = "Baseline (All Human Genes)", 
               hjust = 1.1, vjust = 1.1, 
               color = "grey50", size = 3)
    
    # Store the plot
    plot_values$rate_plot <- p
    p
  })
  
  ##### correlation matrix #####
  output$correlation_matrix <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    # Set up ph object for multiple gene sets
    genes_of_interest <- genes_of_interest_reactive()
    ph$is_of_interest <- FALSE
    ph$annot_of_interest <- "All Human Genes"
    
    # Add each gene set
    for(set_name in names(genes_of_interest)) {
      ph$is_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- TRUE
      ph$annot_of_interest[ph$gene %in% genes_of_interest[[set_name]]] <- set_name
    }

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
    
    # Calculate ratios for each gene set independently
    ratio_data <- ph_by_ps %>%
      group_by(gene_relevance) %>%
      mutate(total_in_group = sum(n)) %>%
      ungroup() %>%
      group_by(mrca_name) %>%
      mutate(
        total_genes = sum(n),
        ratio = n / total_genes
      ) %>%
      mutate(
        normalized_ratio = ifelse(gene_relevance == "All Human Genes", 
                                1,  # Baseline is always 1
                                ratio / (total_in_group / sum(total_in_group)))
      )
    
    # Filter out the baseline and reshape data for correlation matrix
    correlation_data <- ratio_data %>%
      filter(gene_relevance != "All Human Genes") %>%
      select(mrca_name, gene_relevance, normalized_ratio) %>%
      pivot_wider(names_from = gene_relevance, 
                 values_from = normalized_ratio) %>%
      select(-mrca_name)  # Remove mrca_name column immediately after reshaping
    
    # Only show correlation matrix if we have more than one gene set
    if(ncol(correlation_data) > 1) {
      # Create scatter plot matrix
      plot_data <- correlation_data %>%
        as.data.frame() %>%
        select_if(is.numeric)  # Ensure we only have numeric columns
      
      # Create pairs plot
      p <- pairs(plot_data,
            labels = colnames(plot_data),
            main = "Scatter Plot Matrix of Normalized Gene Emergence Rates",
            pch = 16,
            col = rgb(0, 0, 0, 0.5),
            cex.labels = 1.5,  # Increased label size
            cex.axis = 1.2,    # Increased axis text size
            upper.panel = function(x, y, ...) {
              # Check if correlation is significant
              cor_test <- cor.test(x, y)
              p_val <- cor_test$p.value
              if(p_val < 0.05) {
                # Add light blue background for significant correlations
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
                     col = rgb(0.9, 0.9, 1, 0.3))
              }
              points(x, y, ...)
              abline(lm(y ~ x), col = "red", lty = 2)
              # Add correlation coefficient and p-value
              cor_val <- cor_test$estimate
              # Format p-value with scientific notation if needed
              p_text <- ifelse(p_val < 0.001, 
                             sprintf("p = %.1e", p_val),
                             sprintf("p = %.3f", p_val))
              text(mean(range(x)), max(y), 
                   sprintf("r = %.2f\n%s", cor_val, p_text),
                   cex = 1.2, col = "red", font = 2)  # Increased text size and made bold
            },
            lower.panel = function(x, y, ...) {
              # Check if correlation is significant
              cor_test <- cor.test(x, y)
              p_val <- cor_test$p.value
              if(p_val < 0.05) {
                # Add light blue background for significant correlations
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
                     col = rgb(0.9, 0.9, 1, 0.3))
              }
              points(x, y, ...)
              abline(lm(y ~ x), col = "red", lty = 2)
              # Add correlation coefficient and p-value
              cor_val <- cor_test$estimate
              # Format p-value with scientific notation if needed
              p_text <- ifelse(p_val < 0.001, 
                             sprintf("p = %.1e", p_val),
                             sprintf("p = %.3f", p_val))
              text(mean(range(x)), max(y), 
                   sprintf("r = %.2f\n%s", cor_val, p_text),
                   cex = 1.2, col = "red", font = 2)  # Increased text size and made bold
            },
            diag.panel = function(x, ...) {
              usr <- par("usr")
              on.exit(par(usr))
              par(usr = c(usr[1:2], 0, 1.5))
              h <- hist(x, plot = FALSE)
              breaks <- h$breaks
              nB <- length(breaks)
              y <- h$counts
              y <- y/max(y)
              rect(breaks[-nB], 0, breaks[-1], y, col = "grey80")
            })
      
      # Store the plot data and parameters
      plot_values$correlation_matrix <- list(
        data = plot_data,
        labels = colnames(plot_data)
      )
      p
    } else {
      # Show message if only one gene set is selected
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "Select at least two gene sets to view correlation matrix",
                 size = 4) +
        theme_void()
      
      # Store the plot
      plot_values$correlation_matrix <- p
      p
    }
  })
  
  ##### string plot #####
  output$string_plot <- renderPlot({
    if(input$no_isoforms){
      ph <- ph_filt
    }
    
    # Add error handling
    tryCatch({
      # Get genes for current phylostrata
      current_genes <- ph %>% 
        filter(is_of_interest & mrca_name %in% input$ps) %>% 
        pull(gene) %>% 
        unique()
      
      # Check cache first
      cached_results <- get_cached_string_results(current_genes, input$ps)
      
      if (!is.null(cached_results)) {
        # Use cached results
        string_db <- cached_results$string_db
        ph_mapped <- cached_results$ph_mapped
      } else {
        # Generate new results
        string_db <- STRINGdb::STRINGdb$new(species = 9606, version = "12", score_threshold = 200, input_directory = ".")
        ph$is_of_interest <- ph$gene %in% genes_of_interest_reactive()
        ph$annot_of_interest <-
          ifelse(ph$is_of_interest,
                 disease_name(),
                 "All Human Genes")
        
        withProgress(message = "Pulling data from STRINGdb...", value = 0, {
          ph_mapped <- string_db$map(ph %>% filter(is_of_interest & mrca_name %in% input$ps) %>% as.data.frame(), "gene", removeUnmappedRows = T)
          incProgress(amount = 0.5, message = "Generating plot...")
          
          # Cache the results
          save_to_cache(current_genes, input$ps, list(
            string_db = string_db,
            ph_mapped = ph_mapped
          ))
          
          string_db$plot_network(ph_mapped$STRING_id, required_score = 800)
          incProgress(amount = 0.5, message = "Done!")
        })
      }
    }, error = function(e) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, paste("Error generating STRING plot:", e$message), cex = 1.2)
    })
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
  
  # Add download handler
  output$download_plots <- downloadHandler(
    filename = function() {
      paste0("phylostrat_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      # Create PDF with appropriate dimensions
      pdf(file, width = 12, height = 16)
      
      # Create a layout with 3 rows
      layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(1, 1, 1.2))
      
      # Plot 1: Phylostrata plot
      print(plot_values$ph_plot)
      
      # Plot 2: Rate plot
      print(plot_values$rate_plot)
      
      # Plot 3: Correlation matrix
      if (!is.null(plot_values$correlation_matrix)) {
        if (is.list(plot_values$correlation_matrix)) {
          # Recreate the pairs plot
          pairs(plot_values$correlation_matrix$data,
                labels = plot_values$correlation_matrix$labels,
                main = "Scatter Plot Matrix of Normalized Gene Emergence Rates",
                pch = 16,
                col = rgb(0, 0, 0, 0.5),
                cex.labels = 1.5,
                cex.axis = 1.2,
                upper.panel = function(x, y, ...) {
                  cor_test <- cor.test(x, y)
                  p_val <- cor_test$p.value
                  if(p_val < 0.05) {
                    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
                         col = rgb(0.9, 0.9, 1, 0.3))
                  }
                  points(x, y, ...)
                  abline(lm(y ~ x), col = "red", lty = 2)
                  cor_val <- cor_test$estimate
                  p_text <- ifelse(p_val < 0.001, 
                                 sprintf("p = %.1e", p_val),
                                 sprintf("p = %.3f", p_val))
                  text(mean(range(x)), max(y), 
                       sprintf("r = %.2f\n%s", cor_val, p_text),
                       cex = 1.2, col = "red", font = 2)
                },
                lower.panel = function(x, y, ...) {
                  cor_test <- cor.test(x, y)
                  p_val <- cor_test$p.value
                  if(p_val < 0.05) {
                    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
                         col = rgb(0.9, 0.9, 1, 0.3))
                  }
                  points(x, y, ...)
                  abline(lm(y ~ x), col = "red", lty = 2)
                  cor_val <- cor_test$estimate
                  p_text <- ifelse(p_val < 0.001, 
                                 sprintf("p = %.1e", p_val),
                                 sprintf("p = %.3f", p_val))
                  text(mean(range(x)), max(y), 
                       sprintf("r = %.2f\n%s", cor_val, p_text),
                       cex = 1.2, col = "red", font = 2)
                },
                diag.panel = function(x, ...) {
                  usr <- par("usr")
                  on.exit(par(usr))
                  par(usr = c(usr[1:2], 0, 1.5))
                  h <- hist(x, plot = FALSE)
                  breaks <- h$breaks
                  nB <- length(breaks)
                  y <- h$counts
                  y <- y/max(y)
                  rect(breaks[-nB], 0, breaks[-1], y, col = "grey80")
                })
        } else {
          # Print the ggplot object
          print(plot_values$correlation_matrix)
        }
      }
      
      dev.off()
    },
    contentType = "application/pdf"  # Explicitly set content type to PDF
  )
}

##### RUN #####
shinyApp(ui = ui, server = server)
