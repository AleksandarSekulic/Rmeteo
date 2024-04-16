get_meteo <- function(loc,
                     var = "tmean",
                     agg_level = "agg",
                     time_scale = "day",
                     time,
                     from,
                     to,
                     api_key) {
  # check loc - sf, SpatVector or c(lon, lat)
  if (any(class(loc) == "SpatVector")) {
    loc <- crds(loc, df=T)
  } else if (any(class(loc) == "sf")) {
    loc <- st_coordinates(loc)
  } else if (any(class(loc) == "data.frame" | any(class(loc) == "matrix"))) {
    # do nothing
  } else if (any(class(loc) == "numeric" | any(class(loc) == "integer"))) { # vector (numeric or integer)
    loc <- matrix(loc, nrow = 1)
  } else {
    stop('The argument loc must be of sf, SpatVector, data.frame, matrix or vector (numeric or integer) class!') # "STSDF"
  }
  if (missing(var)) {
    stop("The argument var is missing")
  } else {
    if (! var %in% c("tmax", "tmin", "tmean", "prcp", "slp")) {
      stop('The argument var must be "tmax", "tmin", "tmean", "prcp", or "slp"')
    }
  }
  if (missing(agg_level)) {
    stop("The argument agg_level is missing")
  } else {
    if (! agg_level %in% c("agg", "ltm")) {
      stop('The argument agg_level must be "agg" or "ltm"')
    }
  }
  if (missing(time_scale)) {
    stop("The argument time_scale is missing")
  } else {
    if (! time_scale %in% c("day", "mon", "ann")) {
      stop('The argument time_scale must be "day", "mon", or "ann"')
    }
  }
  if (missing(time) & (missing(from) | missing(to))) {
    stop("The argument time or arguments from/to are missing")
  }
  # Construct the URL for the dailymeteo.com API endpoint
  url <- "https://api.dailymeteo.com/meteo/pq/"
  # Set query parameters
  
  data <- c()
  
  for (r in 1:nrow(loc)) {
    print(paste0("Location ", r))
    
    query <- list(
      var = var,
      agg_level = agg_level,
      time_scale = time_scale,
      lat = loc[r, 2],
      lon = loc[r, 1],
      api_key=api_key
    )
    if (missing(time)) {
      query <- append(query, list(from=from, to=to))
    } else {
      if (!is.null(time)){
        query <- append(query, list(time=paste(time, collapse = ",")))
      } else {
        query <- append(query, list(from=from, to=to))
      }
    }
    # Construct the query string
    query_string <- paste0(names(query), "=", unlist(query), collapse = "&")
    full_url <- paste0(url, "?", query_string)
    # Make GET request to the API
    response <- url(full_url)
    # Check if request was successful
    if (inherits(response, "error")) {
      # Print error message if request was not successful
      print("Error: Failed to retrieve data from the API.")
      break
    } else {
      # Read the response and parse JSON
      response_content <- readLines(response, warn = FALSE)
      row_data <- fromJSON(response_content)
      # Close the connection
      close(response)
      if (length(row_data) > 0) {
        row_data <- cbind(loc=r, row_data)
        data <- rbind(data, row_data)
      } else {
        print("Error: No data returned from the API.")
        break
      }
    }
  }
  return(data)
}