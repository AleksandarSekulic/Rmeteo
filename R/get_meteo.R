get_meteo <- function(var = "tmean", agg_level = "agg", time_scale = "day", from=NULL, to=NULL, time=NULL, lat, lon, api_key=NULL) {
  if (is.null(api_key)) {stop("api_key parameter is missing")}
  if (is.null(var)) {stop("var parameter is missing")}
  if (is.null(agg_level)) {stop("agg_level parameter is missing")}
  # if (is.null(time_scale)) {stop("time_scale parameter is missing")}
  if (is.null(time) & (is.null(from) | is.null(to))) {stop("time or from/to parameters are missing")}
  if (is.null(lat) | is.null(lon)) {stop("lat/lon parameters are missing")}
  # Construct the URL for the dailymeteo.com API endpoint
  url <- "https://api.dailymeteo.com/meteo/pq/"
  # Set query parameters
  query <- list(
    var = var,
    agg_level = agg_level,
    time_scale = time_scale,
    lat = lat,
    lon = lon,
    api_key=api_key)
  if (is.null(time)) {
    query <- append(query, list(from=from, to=to))
  } else {
    query <- append(query, list(time=time))
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
    return(NULL)
  } else {
    # Read the response and parse JSON
    response_content <- readLines(response, warn = FALSE)
    data <- fromJSON(response_content)
    data <- data[, c("timestamp", "value")]
    # Close the connection
    close(response)
    # Convert data to data frame
    if (length(data) > 0) {
      # df <- data.frame(date = as.character(names(unlist(data))),
      #                  value = unlist(data))
      return(data)
    } else {
      print("Error: No data returned from the API.")
      return(NULL)
    }
  }
}