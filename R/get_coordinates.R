get_coordinates <- function(location_name = "Belgrade") {
  # Construct the URL for the Nominatim API endpoint
  nominatim_url <- "https://nominatim.openstreetmap.org/search"
  # Set query parameters
  query <- list(
    q = location_name,
    format = "json",
    limit = 1
  )
  # Make GET request to the Nominatim API
  # response <- GET(url, query = query)
  query_string <- paste0(names(query), "=", unlist(query), collapse = "&")
  full_url <- paste0(nominatim_url, "?", query_string)
  response <- url(full_url)
  # Check if request was successful
  if (inherits(response, "error")) {
    stop("Error: Failed to retrieve coordinates.")
  } else {
    # Parse JSON response
    response_content <- readLines(response, warn = FALSE)
    data <- fromJSON(response_content)
    close(response)
    # Extract latitude and longitude from the response
    if (length(data) > 0) {
      lat <- as.numeric(data[1, ]$lat)
      lon <- as.numeric(data[1, ]$lon)
      return(c(lon, lat))
    } else {
      stop("Error: Location not found.")
    }
  }
}
