library(shiny)
library(ComplexHeatmap)
library(colorRamp2)
library(plotly)

# Define the heatmap UI and server as a module
heatmapAppUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file1"), "Upload the TSV file", accept = c(".tsv", ".txt")),
    sliderInput(ns("mean_threshold"), "Minimum Mean Absolute Value", 
                min = 0, max = 10, value = 2, step = 0.01),
    plotlyOutput(ns("heatmap"))
  )
}

heatmapAppServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    data_reactive <- reactive({
      req(input$file1)
      data <- read.delim(input$file1$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      pathway_names <- data$Pathway.Name
      if (anyDuplicated(pathway_names) > 0) {
        pathway_names <- make.names(pathway_names, unique = TRUE)
        data$Pathway.Name <- pathway_names
      }
      data <- data[, -c(1, 2)]
      rownames(data) <- pathway_names
      return(data)
    })
    
    filtered_data_reactive <- reactive({
      data <- data_reactive()
      req(!is.null(data) && nrow(data) > 0 && ncol(data) > 0)
      pathway_means <- rowMeans(abs(data), na.rm = TRUE)
      data[pathway_means >= input$mean_threshold, ]
    })
    
    output$heatmap <- renderPlotly({
      data <- filtered_data_reactive()
      req(nrow(data) > 0)
      
      z <- as.matrix(data)
      
      plot_ly(
        y = rownames(data),
        x = colnames(data),
        z = z,
        colors = colorRamp(c("blue", "white", "red")),
        type = "heatmap",
        showscale = TRUE,
        colorbar = list(title = "Perturbation"),
        zmin = min(z), zmax = max(z),
        hovertemplate = "Pathway: %{y}<br>Experiment: %{x}<br>Value: %{z}<extra></extra>"
      ) %>%
        layout(
          xaxis = list(title = "Experiments"),
          yaxis = list(title = "Pathways"),
          hovermode = "closest"
        )
    })
  })
}