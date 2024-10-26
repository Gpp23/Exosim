library(httr)
library(jsonlite)
library(glue)

sample_DB <- read.table("GSE147507_RawReadCounts_Human.tsv", header = TRUE)

# URL del file remoto
mirandola_url <- "http://mirandola.iit.cnr.it/download/miRandola_version_02_2017.txt"

# Path del file locale
local_file_path <- "local_mirandola_DB.tsv"

# Funzione per gestire la lettura e l'unione dei dati
load_and_merge_data <- function(remote_url, local_path) {
  tryCatch({
    # Leggi il file remoto
    mirandola_DB <- read.table(url(remote_url), header = TRUE, sep = "\t", fill = TRUE)
    
    # Filtra i dati per "exRNA_type"
    mirandola_DB <- mirandola_DB[mirandola_DB$exRNA_type == "exosome", ]
    
    # Controlla se il file locale esiste
    if (file.exists(local_path)) {
      # Leggi il file locale
      local_DB <- read.table(local_path, header = TRUE, sep = "\t", fill = TRUE)
      
      # Unisci i due data frame
      combined_DB <- rbind(mirandola_DB, local_DB)
    } else {
      # Se il file non esiste, crealo con lo stesso header di mirandola_DB
      write.table(mirandola_DB[0, ], file = local_path, sep = "\t", row.names = FALSE, col.names = TRUE)
      combined_DB <- mirandola_DB
    }
    
    return(combined_DB)
  }, error = function(e) {
    cat("Si è verificato un errore durante il tentativo di leggere il file:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Carica e unisci i dati
combined_DB <- load_and_merge_data(mirandola_url, local_file_path)


filter_mirandola_DB <- function(organism = NULL, sample = NULL, pubmedid = NULL) {
  filtered <- combined_DB
  
  if (!is.null(organism)) {
    filtered <- filtered[filtered$organism == organism, ]
  }
  
  if (!is.null(sample)) {
    filtered <- filtered[filtered$sample == sample, ]
  }
  
  if (!is.null(pubmedid)) {
    filtered <- filtered[filtered$PubMed_ID == pubmedid, ]
  }
  
  return(filtered)
}

organisms <- list("Homo sapiens" = "hsa")

read_response <- function(response) {
  content <- content(response, "text")
  
  data <- fromJSON(content)
  
  return(data)
}

api_key <- "DVZ0AGbNbsVWwyDWVm2VxB03BVCPvW2Iu4HB2RwF"

headers <- c(
  "Accept" = "application/json",
  "Authorization" = glue("Bearer {api_key}")
)

set_simulation_input <- function(overexpressed_genes, underexpressed_genes, input_parameter_datapath = NULL) {
  
  if(length(overexpressed_genes) > 0 | length(underexpressed_genes) > 0 ){
    
    overexpressed_input <- NULL
    underexpressed_input <- NULL
    
    if(length(overexpressed_genes) > 0){
      overexpressed_input <- paste(overexpressed_genes, "OVEREXPRESSION", sep = "\t")
      
      overexpressed_input <- paste(overexpressed_input, collapse = "\n")
    }
    
    if(length(underexpressed_genes) > 0){
      underexpressed_input <- paste(underexpressed_genes, "UNDEREXPRESSION", sep = "\t")
      
      underexpressed_input <- paste(underexpressed_input, collapse = "\n")
    }
    
    if (!is.null(input_parameter_datapath)) {
      input_parameter_content <- readLines(input_parameter_datapath)
      
      writeLines(c(overexpressed_input, underexpressed_input, input_parameter_content), temp_file)
    }else{
      temp_file <- tempfile(fileext = ".tsv")
      writeLines(c(overexpressed_input, underexpressed_input), temp_file)
    }
    
    #print(readLines(temp_file))
    
    return(upload_file(temp_file))
  }
  
  return(upload_file(input_parameter_datapath))
}

set_nonExpressedNodes <- function(nonExpressedGenes){
  if(!is.null(nonExpressedGenes)){
    input_nonExpressedGenes <- paste(nonExpressedGenes, collapse = "\n")
    
    print(input_nonExpressedGenes)
    
    temp_file <- tempfile(fileext = ".tsv")
    writeLines(input_nonExpressedGenes, temp_file)
    
    return(upload_file(temp_file))
  }
  return(NULL)
}

set_simulation_parameters <- function(name, organism, epsilon = NULL, seed = NULL, fdr = NULL, reactome = NULL, fast = NULL, miRNAs = NULL, miRNAsEvidence = "STRONG", submit = NULL, nodes_overExpressed = NULL, nodes_underExpressed = NULL, nodes_nonExpressed = NULL, nodes_knockout = NULL, simulationParametersFile = NULL, enrichmentDatabaseFile = NULL, filter = NULL, nonExpressedNodesFile = NULL, knockoutNodesFile = NULL, customNodeTypesFile = NULL, customEdgeTypesFile = NULL, customEdgeSubtypesFile = NULL) {
  # Costruisci i parametri del corpo della richiesta
  simulation_parameters <- list(
    "name" = name,
    "organism" = organism,
    "epsilon" = epsilon,
    "seed" = seed,
    "fdr" = fdr,
    "reactome" = reactome,
    "fast" = fast,
    "miRNAs" = miRNAs,
    "miRNAsEvidence" = miRNAsEvidence,
    "submit" = submit,
    "nodes.overExpressed" = nodes_overExpressed,
    "nodes.underExpressed" = nodes_underExpressed,
    "nodes.nonExpressed" = nodes_nonExpressed,
    "nodes.knockout" = nodes_knockout,
    "simulationParametersFile" = simulationParametersFile,
    "enrichmentDatabaseFile" = enrichmentDatabaseFile,
    "filter" = filter,
    "nonExpressedNodesFile" = nonExpressedNodesFile,
    "knockoutNodesFile" = knockoutNodesFile,
    "customNodeTypesFile" = customNodeTypesFile,
    "customEdgeTypesFile" = customEdgeTypesFile,
    "customEdgeSubtypesFile" = customEdgeSubtypesFile
  )
  
  return(simulation_parameters)
}

create_simulation <- function(simulation_parameters) {
  url <- "https://phensim.tech/api/v1/simulations"
  # Invia la richiesta POST
  response <- POST(
    url,
    body = simulation_parameters,
    add_headers(headers),
    encode = "multipart"
  )
  
  # Verifica la risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    return(read_response(response))
    
  } else {
    return(read_response(response))
  }
}

submit_simulation <- function(SIMULATION_ID){
  url <- glue("https://phensim.tech/api/v1/simulations/{SIMULATION_ID}/submit")
 
   response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    
    return(read_response(response))
  }  
}

list_simulations <- function() {
  
  url <- "https://phensim.tech/api/v1/simulations"
  
  response <- GET(url, 
                  query = list("per_page" = 1000), 
                  add_headers(headers),
                  encode = "json")
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    
    simulations = read_response(response)
    
    return(simulations)
  }
  
}

get_simulation <- function(JOB_ID) {
  
  url <- glue("https://phensim.tech/api/v1/simulations/{JOB_ID}")
  
  response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    #print("Error:", content(response, "text"))
    
    return(NULL)
  } else {
    
    return(read_response(response))
  }
  
}

get_simulation_result <- function(JOB_ID) {
  
  url <- glue("https://phensim.tech/api/v1/simulations/{JOB_ID}/download/output")
  
  response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    return(data.frame(
      "X..Pathway.Id" = character(),
      "Pathway.Name" = character(),
      "Pathway.Activity.Score" = numeric(),
      "Pathway.p.value" = numeric(),
      "Average.Pathway.Perturbation" = numeric(),
      stringsAsFactors = FALSE
    ))
  } else {
    
    result <- read.delim(text = content(response, "text"), sep = "\t")
    
    result <- result[order(result$Pathway.Activity.Score, decreasing = TRUE), ]
    
    return(unique(result[, c("X..Pathway.Id", "Pathway.Name", "Pathway.Activity.Score", "Pathway.p.value", "Average.Pathway.Perturbation")]))
    
  }
  
}

get_simulation_input <- function(JOB_ID) {
  
  url <- glue("https://phensim.tech/api/v1/simulations/{JOB_ID}/download/input_parameters")
  
  response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    
  } else {
    
    result <- content(response, "text")
    
    return(result)
    
  }
  
}

get_combined_simulation_results <- function(job_ids) {
  
  results_list <- list()
  
  for (job_id in job_ids) {
    
    print(paste0("Getting ", job_id, " results"))
    
    data <- get_simulation(job_id)
    
    if(!is.null(data)){
    
      status <- data$data$status
    
      if(status == 3){
        
        input <- get_simulation_input(job_id)
        
        #print(paste0(job_id, " ", input))
        
        if(!is.null(input)){
          
          input <- gsub("\n", "", input)
        
          miRNAs <- strsplit(input, "\r")[[1]]
          
          if(length(miRNAs) > 1){
            name <- get_simulation(job_id)$data$name
            
            name_array <- strsplit(name, " ")[[1]]
            
            group <- name_array[1]
            
            sample <- name_array[5]
            
            result <- get_simulation_result(job_id)
            
            if (nrow(result) > 0) {
              result <- result[, c("X..Pathway.Id", "Pathway.Name", "Average.Pathway.Perturbation")]
              colnames(result)[3] <- paste0(group, "-", sample)
              results_list <- append(results_list, list(result))
            }
            
          }else{
            
            tmp <- strsplit(miRNAs, "\t")[[1]]
            
            miRNA <- tmp[1]
            
            expression <- tmp[2]
            
            if(grepl("UNDEREXPRESSION", expression)){
              expression <- '(-)'
            } else {
              expression <- '(+)'
            }
            
            result <- get_simulation_result(job_id)
            
            if (nrow(result) > 0) {
              result <- result[, c("X..Pathway.Id", "Pathway.Name", "Average.Pathway.Perturbation")]
              colnames(result)[3] <- paste0(miRNA, expression)
              results_list <- append(results_list, list(result))
            }
          }
        }
      }
    }else{
      print(paste0("Skipping: Simulation ",job_id," is not present"))
    }
  }
  
  # Combina i risultati
  if (length(results_list) > 1) {
    print("Merging results")
    combined_results <- Reduce(function(x, y) merge(x, y, by = c("X..Pathway.Id", "Pathway.Name"), all = TRUE), results_list)
    print("Merging completed")
  } else if (length(results_list) == 1) {
    print("Getting result")
    combined_results <- results_list[[1]]
  } else {
    print("Results are empty")
    combined_results <- data.frame()  # Ritorna un data frame vuoto se non ci sono risultati
  }
  
  # Ordina le colonne dei miRNA in ordine lessicografico
  if (ncol(combined_results) > 2) {
    # Ottieni i nomi delle colonne miRNA (escludendo le prime due)
    mirna_columns <- colnames(combined_results)[-(1:2)]
    
    # Ordina i nomi delle colonne in ordine lessicografico
    sorted_columns <- sort(mirna_columns)
    
    # Ricostruisci il data frame con le colonne ordinate
    combined_results <- combined_results[, c(colnames(combined_results)[1:2], sorted_columns)]
  }
  
  return(combined_results)
}




delete_simulation <- function(JOB_ID) {
  
  url <- glue("https://phensim.tech/api/v1/simulations/{JOB_ID}")
  
  response <- DELETE(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    
    return(read_response(response))
  }
}

run_simulation <- function(input, name, sample, titles, description, overExpressedMiRNA, underExpressedMiRNA, organism, epsilon, randomSeed, fdr, miRNAsEvidence, inputParameters) {
  if (is.null(name)) {
    name <- paste(titles, sample, description, sep = " - ")
  }
  
  # Inizializzare nonExpressedGenes
  nonExpressedGenes <- NULL
  
  # Controllare se input$selected_sample è una colonna di sample_DB
  if (gsub("-", ".",sample) %in% colnames(sample_DB)) {
    
    # Filtrare i geni con espressione < 10 nella colonna specificata da input$selected_sample
    nonExpressedGenes <- sample_DB$gene[sample_DB[[sample]] == 1]
  }
  
  if (length(overExpressedMiRNA) + length(underExpressedMiRNA) == 1){
    tmpMiRNA <- NULL
    if (!is.null(overExpressedMiRNA)) {
      tmpMiRNA <- overExpressedMiRNA
      name <- paste(name, tmpMiRNA, "(+) - ")
    }
    else {
      tmpMiRNA <- underExpressedMiRNA
      name <- paste(name, tmpMiRNA, "(-) - ")
      }
  }
  
  simulation_parameters <-
    set_simulation_parameters(
      name = name,
      organism = unlist(organisms[organism])[1],
      epsilon = epsilon,
      seed = randomSeed,
      fdr = fdr,
      reactome = 1,
      fast = NULL,
      miRNAs = NULL,
      miRNAsEvidence = miRNAsEvidence,
      submit = NULL,
      simulationParametersFile = set_simulation_input(overExpressedMiRNA, underExpressedMiRNA, inputParameters$datapath),
      nonExpressedNodesFile = set_nonExpressedNodes(nonExpressedGenes)
    )
  
  simulation <-
    create_simulation(simulation_parameters = simulation_parameters)
  
  if (!is.null(simulation$data)) {
    result <- submit_simulation(simulation$data$id)
    
    return(result)
    
  } else{
    showModal(modalDialog(title = "Simulation Failed",
                          simulation,
                          easyClose = TRUE))
    return(NULL)
  }
}




