#'@name Nexus
#'@title Nexus
#'@description Ancestor class for time series online event detection
#'@return Nexus object
#'@examples onelin_detector <- nexus(datasource, fbiad())
#'@export
nexus <- function(datasource, detector, warm_size = 30, batch_size = 30, mem_batches = 0) {
  obj <- dal_base()
  obj$detection <- NULL
  obj$stable_detection <- NULL
  obj$datasource <- datasource
  obj$detector <- detector
  obj$warm_size <- warm_size
  obj$batch_size <- batch_size
  obj$mem_size <- batch_size * mem_batches
  obj$base_pos <- 1
  obj$serie <- NULL
  obj$stable_serie <- NULL
  class(obj) <- append("nexus", class(obj))
  return(obj)
}

#'@name nex_datasource
#'@title nex_datasource
#'@description Function to create an object with series data
#'@return Dal object with data
#'@examples onelin_detector <- nexus(datasource, fbiad())
#'@export
nex_datasource <- function(name) {
  obj <- dal_base()
  obj$name <- name
  class(obj) <- append("nex_datasource", class(obj))
  return(obj)
}

#'@name getnext
#'@title getnext
#'@export
getnext <- function(obj, ...) {
  UseMethod("getnext")
}

#'@export
getnext.default <- function(obj) {
  return(obj)
}

#'@name fetch
#'@title fetch
#'@export
fetch <- function(obj, ...) {
  UseMethod("fetch")
}

#'@export
fetch.default <- function(obj) {
  return(NA)
}

#'@name nex_simulated_datasource
#'@title nex_simulated_datasource
#'@export
nex_simulated_datasource <- function(name, serie) {
  obj <- nex_datasource(name)
  obj$pos <- 0
  obj$serie <- serie
  obj$name <- name
  class(obj) <- append("nex_simulated_datasource", class(obj))
  return(obj)
}

#'@export
getnext.nex_simulated_datasource <- function(obj) {
  result <- NULL
  if (obj$pos < length(obj$serie)) {
    obj$pos <- obj$pos + 1
    result <- obj
  }
  return(result)
}

#'@export
fetch.nex_simulated_datasource <- function(obj) {
  return(result <- obj$serie[obj$pos])
}

#'@name warmup
#'@title warmup
#'@export
warmup <- function(obj, ...) {
  UseMethod("warmup")
}

#'@export
warmup.default <- function(obj) {
  return(obj)
}

#'@export
warmup.nexus <- function(obj) {
  while (length(obj$serie) < obj$warm_size - 1) {
    obj$datasource <- getnext(obj$datasource)
    if (is.null(obj$datasource))
      break
    obj$serie <- c(obj$serie, fetch(obj$datasource))
  }
  return(obj)
}

#'@name adjust_memory
#'@title adjust_memory
#'#'@description adjusts how batch memory will occur
#'@export
adjust_memory <- function(obj) {
  if (obj$mem_size == 0)
    return(obj) # full memory
  n <- length(obj$serie)
  if (n >= obj$mem_size) {
    if(!is.null(obj$detection) & length(obj$detection)>0){
      obj$stable_detection <-rbind(obj$stable_detection, obj$detection[1:obj$batch_size, ])
      obj$stable_detection$idx <- 1:nrow(obj$stable_detection)
    }
    obj$stable_serie <- c(obj$stable_serie, obj$serie[1:obj$batch_size])
    obj$serie <- obj$serie[(obj$batch_size+1):n]
  }
  return(obj)
}

#'@name detect
#'@title detect
#'#'@description Slides through time series for online detection of events.
#'@export
detect.nexus <- function(obj) {
  obj$datasource <- getnext(obj$datasource)
  if (!is.null(obj$datasource)) {
    obj$serie <- c(obj$serie, fetch(obj$datasource))
    if ((length(obj$serie) %% obj$batch_size) == 0) {
      obj$detector <- fit(obj$detector, obj$serie)
      obj <- adjust_memory(obj)
    }
    idxref <- 0
    if (!is.null(obj$stable_detection))
      idxref <- nrow(obj$stable_detection)
    detection <- detect(obj$detector, obj$serie)
    detection$idx <- detection$idx + idxref
    detection$event <- as.integer(detection$event)
    ##Under development
    #Probability of the observation being an event
    detection$ef <- 0 #Events frequency of a series point
    detection$bf <- 0 #Batch frequency of a series point
    detection$pe <- detection$ef / detection$bf
    ##Under development
    #Detection temporal lag
    detection$first_batch <- 0 #First batch where the event was detected
    detection$last_batch <- 0 #Last batch where the event was detected
    
    obj$detection <- rbind(obj$stable_detection, detection)
  }
  return(obj)
}
