library(shiny)
library(shinyWidgets)
library(DT)
source("global_v2.R")

ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    titlePanel("Phensim parameters"),
    
    pickerInput(
      inputId = "selected_organism",
      label = "Organism:",
      multiple = FALSE,
      choices = NULL,
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE,
        noneSelectedText = "Nothing selected"
      )
    ),
    
    selectizeInput(
      inputId = "selected_pubmedid",
      label = "PubMed Id:",
      choices = NULL,
      options = list(
        placeholder = "Nothing selected",
        allowEmptyOption = TRUE
      )
    ),
    
    pickerInput(
      inputId = "selected_sample",
      label = "Sample:",
      multiple = FALSE,
      choices = NULL,
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE,
        noneSelectedText = "Nothing selected"
      )
    ),
    
    pickerInput(
      inputId = "selected_overexpressed_miRNA",
      label = "OverExpressed MiRNA:",
      multiple = TRUE,
      choices = NULL,
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE,
        noneSelectedText = "Nothing selected"
      )
    ),
    
    pickerInput(
      inputId = "selected_underexpressed_miRNA",
      label = "UnderExpressed MiRNA:",
      multiple = TRUE,
      choices = NULL,
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE,
        noneSelectedText = "Nothing selected"
      )
    ),
    
    selectInput(
      "fdr",
      label = "FDR method:",
      choices = c(
        "Benjamini & Hochberg" = "BH",
        "Q-value (Storey et al.)" = "QV",
        "Local FDR (Efron et al.)" = "LOC"
      ),
      selected = "BH"
    ),
    
    numericInput(
      "epsilon",
      label = "Epsilon:",
      value = 0.001,
      min = 0,
      step = 0.001,
      width = "300px"
    ),
    
    numericInput(
      "random_seed",
      label = "Random Seed:",
      value = 780,
      width = "300px"
    ),
    
    pickerInput(
      inputId = "selected_miRNAsEvidence",
      label = "miRNAsEvidence",
      choices = c("STRONG", "WEAK", "PREDICTION"),
      multiple = FALSE,
      selected = "STRONG",
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE
      )
    ),
    
    pickerInput(
      inputId = "selected_submitType",
      label = "Simulation type",
      choices = c("COMPACT", "SPLITTED"),
      multiple = FALSE,
      selected = "COMPACT",
      options = list(
        `actions-box` = TRUE,
        liveSearch = TRUE,
        liveSearchPlaceholder = TRUE
      )
    ),
    fileInput("input_parameters", "Input Parameters"),
    textInput(
      inputId = "description",
      label = "Description"
    ),
    actionButton(inputId = "run_query_simulation", label = "Submit simulation")
  ),
  mainPanel(
    titlePanel("Phensim simulations"),
    DTOutput('data'),
    titlePanel("Pathway results"),
    DTOutput('result'),
    actionButton(inputId = "showGenes", label = "Show expressed genes"),
    actionButton(inputId = "showParameters", label = "Show parameters"),
    downloadButton(outputId = "downloadResult", label = "Download results"),
    actionButton(inputId = "open_website", label = "More details"),
    titlePanel("Merge results"),
    textInput("job_id_range", "Insert Job IDs (es. 1001-1005,1008):", ""),
    downloadButton("download_combined_tsv", "Merge & download")
  ),
))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  global_values <-
    reactiveValues(titles = NULL,
                   url = NULL,
                   selected_id = NULL,
                   miRNAs = NULL,)
  
  observe({
    updatePickerInput(session,
                      "selected_organism",
                      choices = sort(unique(combined_DB$organism), decreasing = TRUE))
  })
  
  observeEvent(input$selected_organism, {
    organism <- input$selected_organism
    
    db <- filter_mirandola_DB(organism = organism)
    
    db <- unique(db[, c('PubMed_ID', 'disease_or_cell_line')])
    
    db <- db[order(db$PubMed_ID),]
    
    pubmedids <- db$PubMed_ID
    
    disease_or_cell_lines <- db$disease_or_cell_line
    
    titles <-
      paste(disease_or_cell_lines, " (", pubmedids, ")", sep = "")
    
    names(titles) <- pubmedids
    
    global_values$titles <- titles
    
    names(pubmedids) <- unlist(titles)
    
    updateSelectizeInput(session, "selected_pubmedid", choices = pubmedids)
    
  })

  
  observeEvent(input$selected_pubmedid, {
    organism <- input$selected_organism
    
    #sample <- input$selected_sample
    
    pubmedid <- input$selected_pubmedid
    
    if(pubmedid == ""){
      pubmedid <- NULL
    }
    
    #print(pubmedid)
    
    db <- filter_mirandola_DB(organism = organism, pubmedid = pubmedid)
    
    samples <- sort(unique(db$sample), decreasing = TRUE)
    
    updatePickerInput(session, "selected_sample", choices = samples)
    
  })
  
  observeEvent(input$selected_sample, {
    organism <- input$selected_organism
    
    sample <- input$selected_sample
    
    pubmedid <- input$selected_pubmedid
    
    db <- filter_mirandola_DB(organism = organism, sample = sample, pubmedid = pubmedid)
    
    #print(unique(db$organism))
    
    global_values$miRNAs <-
      sort(unique(
        filter_mirandola_DB(
          organism = organism,
          sample = sample,
          pubmedid = pubmedid
        )$miRBase_Last_Version
      ))
    
    updatePickerInput(session, "selected_overexpressed_miRNA", choices = global_values$miRNAs)
    updatePickerInput(session, "selected_underexpressed_miRNA", choices = global_values$miRNAs)
    
  })
  
  observe({
    overexpressed <- input$selected_overexpressed_miRNA
    underexpressed <- input$selected_underexpressed_miRNA
    
    remaining_miRNAs <- setdiff(global_values$miRNAs, union(overexpressed, underexpressed))
    
    updatePickerInput(session, "selected_overexpressed_miRNA", choices = c(overexpressed, remaining_miRNAs), selected = overexpressed)
    updatePickerInput(session, "selected_underexpressed_miRNA", choices = c(underexpressed, remaining_miRNAs), selected = underexpressed)
  })
  
  observeEvent(input$run_query_simulation, {
    if (input$selected_pubmedid == "") {
      showModal(modalDialog(
        title = "Enter Title",
        textInput("name", "Title:", ""),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("submit_text", "Submit")
        )
      ))
    } else {
      if (input$selected_submitType == "COMPACT"){
        result <- run_simulation(input, 
                       NULL,
                       input$selected_sample,
                       global_values$titles[input$selected_pubmedid],
                       input$description,
                       input$selected_overexpressed_miRNA,
                       input$selected_underexpressed_miRNA,
                       input$selected_organism,
                       input$epsilon,
                       input$randomSeed,
                       input$fdr,
                       input$selected_miRNAsEvidence,
                       input$input_parameters)
        
        showModal(modalDialog(title = "Simulation Started",
                              result,
                              easyClose = TRUE))
        
        output$data <- renderDT({
          df <- list_simulations()$data[, c('id', 'name', 'readable_status')]
          datatable(
            df,
            options = list(
              lengthChange = FALSE,
              rowId = ~ id,
              order = list(list(0, 'desc')),
              pageLength = 5
            ),
            selection = "single",
            rownames = FALSE
          )})
          
        print("Submitted compact simulation")
      }
      else{
        withProgress(message = 'Starting splitted simulations...', value = 0, {
          total_mirnas <- length(input$selected_overexpressed_miRNA) + length(input$selected_underexpressed_miRNA)
          
          # Ciclo per miRNA sovraespressi
          for (mirna in input$selected_overexpressed_miRNA) {
            run_simulation(input, 
                           NULL,
                           input$selected_sample,
                           global_values$titles[input$selected_pubmedid],
                           input$description,
                           mirna,
                           NULL,
                           input$selected_organism,
                           input$epsilon,
                           input$randomSeed,
                           input$fdr,
                           input$selected_miRNAsEvidence,
                           input$input_parameters)
            
            # Incrementa la barra di progresso
            incProgress(1 / total_mirnas)
          }
          
          # Ciclo per miRNA sottoespressi
          for (mirna in input$selected_underexpressed_miRNA) {
            run_simulation(input, 
                           NULL,
                           input$selected_sample,
                           global_values$titles[input$selected_pubmedid],
                           input$description,
                           NULL,
                           mirna,
                           input$selected_organism,
                           input$epsilon,
                           input$randomSeed,
                           input$fdr,
                           input$selected_miRNAsEvidence,
                           input$input_parameters)
            
            # Incrementa la barra di progresso
            incProgress(1 / total_mirnas)
          }
        })
        
        showModal(modalDialog(title = "Simulations Started",
                              paste0(
                                paste0(
                                  "Overexpressed miRNAs: ",
                                  paste0(input$selected_overexpressed_miRNA, collapse = "- ")
                                ),
                                paste0(
                                  "Underexpressed miRNAs: ",
                                  paste0(input$selected_underexpressed_miRNA, collapse = "- ")
                                ),
                                collapse = "\n"
                              ),
                              easyClose = TRUE))
        
        output$data <- renderDT({
          df <- list_simulations()$data[, c('id', 'name', 'readable_status')]
          datatable(
            df,
            options = list(
              lengthChange = FALSE,
              rowId = ~ id,
              order = list(list(0, 'desc')),
              pageLength = 5
            ),
            selection = "single",
            rownames = FALSE
          )
        })
        
        print("Submitted splitted run")
      } 
    }
  })
  
  observeEvent(input$submit_text, {
    name <- input$name
    if (!is.null(name) && name != "") {
      removeModal()
      result <- run_simulation(input, 
                     name,
                     input$selected_sample,
                     global_values$titles[input$selected_pubmedid],
                     input$description,
                     input$selected_overexpressed_miRNA,
                     input$selected_underexpressed_miRNA,
                     input$selected_organism,
                     input$epsilon,
                     input$randomSeed,
                     input$fdr,
                     input$selected_miRNAsEvidence,
                     input$input_parameters)
      
      showModal(modalDialog(title = "Simulation Started",
                            result,
                            easyClose = TRUE))
      
      output$data <- renderDT({
        df <- list_simulations()$data[, c('id', 'name', 'readable_status')]
        datatable(
          df,
          options = list(
            lengthChange = FALSE,
            rowId = ~ id,
            order = list(list(0, 'desc')),
            pageLength = 5
          ),
          selection = "single",
          rownames = FALSE
        )
      })
    }
  })
  
  
  output$data <- renderDT({
    df <- list_simulations()$data[, c('id', 'name', 'readable_status')]
    datatable(
      #df[order(df$id, decreasing = TRUE), ],
      df,
      options = list(
        lengthChange = FALSE,
        rowId = df$id,
        order = list(list(0, 'desc')),
        pageLength = 5
      ),
      selection = "single",
      rownames = FALSE
    )
  })
  
  
  output$result <- renderDT({
    dummy <- data.frame(
      "X..Pathway.Id" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
      "Pathway.Name" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
      "Pathway.Activity.Score" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
      "Pathway.p.value" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
      "Average.Pathway.Perturbation" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
      stringsAsFactors = FALSE
    )
    
    datatable(
      dummy,
      options = list(lengthChange = FALSE,
                     pageLength = 5),
      rownames = FALSE
    )
  })
  
  observeEvent(input$data_rows_selected, {

    selected_row <- input$data_rows_selected
    
    if (length(selected_row) == 0) {
      output$result <- renderDT({
        dummy <- data.frame(
          "X..Pathway.Id" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
          "Pathway.Name" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
          "Pathway.Activity.Score" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
          "Pathway.p.value" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
          "Average.Pathway.Perturbation" = c("N/a", "N/a", "N/a", "N/a", "N/a"),
          stringsAsFactors = FALSE
        )
        
        datatable(
          dummy,
          options = list(lengthChange = FALSE,
                         pageLength = 5),
          rownames = FALSE
        )
      })
    } else{
      simulation_id <- list_simulations()$data[selected_row, "id"]
      
      global_values$selected_id <- simulation_id
      
      print(simulation_id)
      
      result <- get_simulation_result(simulation_id)
      
      if (length(result) == 0) {
        global_values$url <- NULL
      } else{
        global_values$url <-
          glue("https://phensim.tech/simulations/{simulation_id}")
      }
      
      output$result <- renderDT({
        datatable(
          result,
          options = list(lengthChange = FALSE,
                         pageLength = 5),
          rownames = FALSE
        )
      })
      
      output$expressed_genes <- renderDT({
        result <- get_simulation_input(global_values$selected_id)
        
        result <-
          read.table(textConnection(result),
                     header = FALSE,
                     sep = "\t")
        
        names(result) <- c('Gene', 'Expression')
        
        
        datatable(
          result,
          options = list(lengthChange = FALSE,
                         pageLength = 5),
          rownames = FALSE
        )
      })
      
      output$parameters <- renderDT({
        parameters <- get_simulation(global_values$selected_id)
        parameters <- parameters$data$parameters
        
        parameters <- data.frame(
          Parameter = c(
            "fast",
            "fdr",
            "epsilon",
            "seed",
            "reactome",
            "enrichMiRNAs",
            "miRNAsEvidence"
          ),
          Value = c(
            ifelse(is.null(parameters$fast), "N/a", parameters$fast),
            ifelse(is.null(parameters$fdr), "N/a", parameters$fdr),
            ifelse(
              is.null(parameters$epsilon),
              "N/a",
              parameters$epsilon
            ),
            ifelse(is.null(parameters$seed), "N/a", parameters$seed),
            ifelse(
              is.null(parameters$reactome),
              "N/a",
              parameters$reactome
            ),
            ifelse(
              is.null(parameters$enrichMiRNAs),
              "N/a",
              parameters$enrichMiRNAs
            ),
            ifelse(
              is.null(parameters$miRNAsEvidence),
              "N/a",
              parameters$miRNAsEvidence
            )
          )
        )
        datatable(parameters, rownames = FALSE)
      })
    }
    
  })
  
  observeEvent(input$showGenes, {
    if (!is.null(global_values$selected_id)) {
      showModal(
        modalDialog(
          title = "Simulation expressed genes",
          DT::dataTableOutput("expressed_genes"),
          easyClose = TRUE
        )
      )
      
    } else{
      showModal(
        modalDialog(
          title = "Simulation expressed genes",
          "Nessuna simulazione selezionata",
          easyClose = TRUE
        )
      )
    }
  })
  
  observeEvent(input$showParameters, {
    if (!is.null(global_values$selected_id)) {
      showModal(
        modalDialog(
          title = "Simulation parameters",
          DT::dataTableOutput("parameters"),
          easyClose = TRUE
        )
      )
      
    } else{
      showModal(
        modalDialog(
          title = "Simulation parameters",
          "Select at least one simulation",
          easyClose = TRUE
        )
      )
    }
  })
  
  output$downloadResult <- downloadHandler(
    # La funzione filename non avvia il download se il job id non è selezionato
    filename = function() {
      if (!is.null(global_values$selected_id)) {
        paste("simulation_", global_values$selected_id, "_results.tsv", sep = "")
      } else {
        showModal(
          modalDialog(
            title = "Failed",
            "Select at least one simulation",
            easyClose = TRUE
          )
        )
        return(NULL)  # Interrompe l'azione se non c'è un job id
      }
    },
    
    content = function(file) {
      if (!is.null(global_values$selected_id)) {
        # Aggiunta della barra di progresso
        withProgress(message = 'Download in progress...', value = 0, {
          
          # Incremento la barra di progresso al 25%
          incProgress(0.25)
          
          # Simulazione del processo di combinazione dei risultati (potrebbe essere un'operazione costosa)
          combined_df <- get_combined_simulation_results(global_values$selected_id)
          
          # Incremento la barra di progresso al 75%
          incProgress(0.75)
          
          # Scrivo i risultati nel file
          write.table(combined_df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
          
          # Incremento la barra di progresso al 100%
          incProgress(1)
        })
      } else {
        showModal(
          modalDialog(
            title = "Failed",
            "Select at least one simulation",
            easyClose = TRUE
          )
        )
        return()  # Ferma l'azione se non c'è un job id selezionato
      }
    }
  )
  

  
  
  observeEvent(input$open_website, {
    hyperlink <-
      HTML(glue("<a href='{global_values$url}'>{global_values$url}</a>"))
    showModal(modalDialog(title = "Simulation URL",
                          hyperlink,
                          easyClose = TRUE))
  })
  
  parse_job_id_range <- function(job_id_range_str) {
    job_ids <- unlist(strsplit(job_id_range_str, ","))
    job_ids <- unlist(lapply(job_ids, function(x) {
      if (grepl("-", x)) {
        range_vals <- as.numeric(unlist(strsplit(x, "-")))
        return(seq(range_vals[1], range_vals[2]))
      } else {
        return(as.numeric(x))
      }
    }))
    return(job_ids)
  }

  output$download_combined_tsv <- downloadHandler(
    filename = function() {
      paste("combined_results_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Download in progress...', value = 0, {
        job_ids <- parse_job_id_range(input$job_id_range)
        incProgress(0.25)
        
        combined_df <- get_combined_simulation_results(job_ids)
        incProgress(0.50)
        
        if (nrow(combined_df) == 0) {
          stop("Result is empty.")
        }
        
        write.table(combined_df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        incProgress(1)
      })
    }
  )
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
