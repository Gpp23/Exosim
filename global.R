library(httr)
library(jsonlite)
library(glue)

tryCatch({
  
  mirandola_DB <- read.table(url("http://mirandola.iit.cnr.it/download/miRandola_version_02_2017.txt"), header = TRUE, sep = "\t", fill = TRUE)
  
  mirandola_DB <- mirandola_DB[mirandola_DB$exRNA_type == "exosome", ]
}, error = function(e) {
  cat("Si è verificato un errore durante il tentativo di leggere il file:", conditionMessage(e), "\n")
})

filter_mirandola_DB <- function(organism = NULL, sample = NULL, pubmedid = NULL) {
  filtered <- mirandola_DB
  
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

headers <- c(
  "Accept" = "application/json",
  "Authorization" = "Bearer DVZ0AGbNbsVWwyDWVm2VxB03BVCPvW2Iu4HB2RwF"
)

set_simulation_input <- function(genes, expression) {
  
  simulation_input <- paste(genes, expression, sep = "\t")
  
  simulation_input <- paste(simulation_input, collapse = "\n")
  
  return(simulation_input)
}

set_simulation_parameters <- function(name, organism, simulationParams, fdr = "BH", epsilon = 0.1, seed = 780, metaPathway = TRUE) {
  simulation_parameters <- list(
    name = name,
    organism = organism,
    "simulation-input" = simulationParams,
    "enrich-db" = NULL,
    "db-filter" = NULL,
    "node-types" = NULL,
    "custom-edge-types" = NULL,
    "custom-edge-subtypes" = NULL,
    "fdr" = fdr,
    epsilon = epsilon,
    "random-seed" = seed,
    "enrich-mirnas" = 'on',
    metaPathway = metaPathway
  )
  
  return(simulation_parameters)
}

query_simulation <- function(query_simulation_parameters) {
  simulations <- list_simulations()
  
  for (i in 1:nrow(simulations)) {
    simulation <- simulations[i, ]
    
    simulation_parameters <- get_simulation_parameters(simulation$id)
    
    print(simulation_parameters)
  }
}


run_simulation <- function(simulation_parameters) {
  
  url <- "https://phensim.tech/api/v1/simulations"
  
  response <- POST(url, body = simulation_parameters, add_headers(headers))
  
  if (http_type(response) == "application/json" && http_error(response)) {
    return(content(response, "text"))
  } else {
    return(paste("Simulation started with success with id:", read_response(response)$id))
  }
}

list_simulations <- function() {
  
  url <- "https://phensim.tech/api/v1/simulations"
  
  response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    
    simulations = read_response(response)
    
    return(data.frame(simulations))
  }
  
}

get_simulation_details <- function(JOB_ID) {
  
  url <- glue("https://phensim.tech/api/v1/simulations/{JOB_ID}")
  
  response <- GET(url, add_headers(headers))
  
  # Controlla lo stato della risposta
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    
    return(read_response(response))
  }
  
}

get_simulation_parameters <- function(JOB_ID) {
  
  url <- glue("https://phensim.atlas.dmi.unict.it/api/simulations/{JOB_ID}/parameters")
  
  response <- GET(url, add_headers(headers))
  
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    return(read_response(response))
  }
}

get_raw_simulation_results <- function(JOB_ID) {
  
  url <- glue("https://phensim.atlas.dmi.unict.it/api/simulations/<JOB_ID>/results/raw")
  
  response <- GET(url, add_headers(headers))
  
  if (http_type(response) == "application/json" && http_error(response)) {
    stop("Error:", content(response, "text"))
  } else {
    return(read_response(response))
  }
  
}




