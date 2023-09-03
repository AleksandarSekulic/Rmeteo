get_meteo <- function(loc,
                      dates,
                      var = "tmean", # "tmax" "tmin" "prcp" "slp"
                      source = "MeteoEurope1km"){
  # check loc - sf, SpatVector or c(lon, lat)
  if (any(class(loc) == "SpatVector")) {
    # do nothing
  } else if (any(class(loc) == "sf")) {
    loc <- vect(loc)
  } else if (any(class(loc) == "data.frame" | any(class(loc) == "matrix"))) {
    loc <- vect(as.matrix(loc), crs="+proj=longlat +datum=WGS84")
  } else if (any(class(loc) == "numeric" | any(class(loc) == "integer"))) { # vector (numeric or integer)
    loc <- vect(matrix(loc, nrow = 1), crs="+proj=longlat +datum=WGS84")
  } else {
    stop('The argument loc must be of sf, SpatVector, data.frame, matrix or vector (numeric or integer) class!') # "STSDF"
  }
  # to 3035
  loc <- terra::project(loc, "EPSG:3035")
  
  # check var
  if (! var %in% c("tmax", "tmin", "tmean", "prcp", "slp")) {
    stop('The argument var must be "tmax", "tmin", "tmean", "prcp", or "slp"')
  }
  
  # check dates
  if(any(class(dates) == "Date")) {
    dates <- as.character(dates)
  } else {
    test_date <- tryCatch(!is.na(as.Date(dates)),  
                          error = function(err) {FALSE})
    if (all(test_date)) {
      dates <- as.Date(dates)
    } else {
      stop('The argument dates must be of Date or character (YYYY-MM-DD) class')
    }
  }
  dates1 <- gsub("-", "", dates, fixed=T)
  
  # create links
  # http://dailymeteo.com/cog/europe/tmean/tmean_day_19910101_3035.tif
  tif.cog <- sapply(dates1, function(d) paste0("/vsicurl/https://dailymeteo.com/cog/europe/", var, "/", var, "_day_", d, "_3035.tif"))
  # terra query
  r = rast(tif.cog)
  result <- extract(r, loc)
  result[, 2:length(result)] <- result[,2:length(result)] / 10
  names(result) <- c("ID", as.character(dates))
  
  # results
  return(result)
}