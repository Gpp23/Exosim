library(shiny)
library(shinyWidgets)
library(DT)
source("global.R")

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(titlePanel("Phensim parameters"),
                    
                    pickerInput(
                      inputId = "selected_organism",
                      label = "Organism:",
                      multiple = FALSE,
                      choices = sort(unique(mirandola_DB$organism), decreasing = TRUE),
                      options = list(
                        `actions-box` = TRUE,
                        liveSearch = TRUE,
                        liveSearchPlaceholder = TRUE,
                        noneSelectedText = "Nothing selected"
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
                      inputId = "selected_pubmedid",
                      label = "PubMed Id:",
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
                      inputId = "selected_miRNA",
                      label = "MiRNA:",
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
                      inputId = "selected_expression",
                      label = "Genes expression:",
                      multiple = FALSE,
                      choices = c("IN (-)" = 'UNDEREXPRESSION', "OUT (+)" = 'OVEREXPRESSION'),
                      options = list(
                        `actions-box` = TRUE,
                        liveSearch = TRUE,
                        liveSearchPlaceholder = TRUE,
                        noneSelectedText = "Nothing selected"
                      )
                    ),
                    
                    selectInput("fdr", label = "FDR method:",
                                choices = c("Benjamini & Hochberg" = "BH",
                                            "Q-value (Storey et al.)" = "QV",
                                            "Local FDR (Efron et al.)" = "LOC"),
                                selected = "BH"),
                    
                    numericInput("epsilon", label = "Epsilon:", 
                                 value = 0.001, min = 0, step = 0.001, width = "300px"),
                    
                    numericInput("random_seed", label = "Random Seed:", value = NULL, width = "300px"),
                    
                    checkboxInput(inputId = "selected_metaPathway", label = "metaPathway", value = TRUE),
                    
                    actionButton(inputId = "run_query_simulation", label = "Submit"),),
    mainPanel(titlePanel("Phensim simulations"),
      DTOutput('data'))
  )
  
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  global_values <- reactiveValues(titles = NULL)
  
  
  observe({
    organism <- input$selected_organism

    samples <- sort(unique(filter_mirandola_DB(organism = organism)$sample))
    
    updatePickerInput(session, "selected_sample", choices = samples)
  })
  
  observe({
    organism <- input$selected_organism
    
    sample <- input$selected_sample
    
    db <- (unique(filter_mirandola_DB(organism = organism, sample = sample)[, c('PubMed_ID', 'description')]))
    
    db <- db[order(db$PubMed_ID), ]
    
    pubmedids <- db$PubMed_ID
    descriptions <- db$description
    
    titles <- paste(descriptions, " (", pubmedids, ")", sep = "")
    
    names(titles) <- pubmedids
    
    global_values$titles <- titles
    
    names(pubmedids) <- unlist(titles)
    
    #print(titles)
    
    updatePickerInput(session, "selected_pubmedid", choices = pubmedids)
  })
  
  observe({
    organism <- input$selected_organism
    
    sample <- input$selected_sample
    
    pubmedid <- input$selected_pubmedid
    
    miRNAs <- sort(unique(filter_mirandola_DB(organism = organism, sample = sample, pubmedid = pubmedid)$miRBase_Last_Version))
    
    updatePickerInput(session, "selected_miRNA", choices = miRNAs)
  })
  
  observeEvent(input$run_query_simulation, {

    simulation_parameters <- set_simulation_parameters(name = global_values$titles[input$selected_pubmedid],
                                                       organism = unlist(organisms[input$selected_organism])[1],
                                                       simulationParams = set_simulation_input(input$selected_miRNA, input$selected_expression),
                                                       fdr = input$fdr, 
                                                       epsilon = input$epsilon, 
                                                       seed = input$random_seed,
                                                       metaPathway = input$selected_metaPathway)
    
    #query_simulation(simulation_parameters)
    simulation <- run_simulation(simulation_parameters = simulation_parameters)
    
    showModal(modalDialog(
      title = "Simulation Result",
      simulation,
      easyClose = TRUE
    ))
    
    
    output$data <- renderDT({
      datatable(list_simulations()[, c('id', 'job_name', 'job_status', 'uri')], options = list(lengthChange = FALSE))
    })
    
  })
  
  output$data <- renderDT({
    datatable(list_simulations()[, c('id', 'job_name', 'job_status', 'uri')], options = list(lengthChange = FALSE))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
