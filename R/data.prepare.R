data.prepare <- function (data,
                          data.staid.x.y.z=NULL,
                          obs.col=NULL,
                          s.crs=NA
) {
  # terra
  if (any(class(data) == "SpatVector")) {
    data.df = as.data.frame(data)
    coords = crds(data, df=T)
    data.df = cbind(data.df, coords)
    x.y.loc <- match(dimnames(coords)[[2]], names(data.df))
    ############## OVO JE ZA PREDICTION, TUNE i CV
    if (is.null(data.staid.x.y.z)) {
      if ("staid" %in% names(data.df)) {
        staid.loc <- match("staid", names(data.df))
      } else {
        data.df$staid <- 1:nrow(data.df)
        staid.loc <- length(data.df)
      }
      data.staid.x.y.z <- c(staid.loc,x.y.loc,NA)
    } else {
      if (!is.numeric(data.staid.x.y.z)) {
        data.staid.x.y.z <- match(data.staid.x.y.z, names(data.df)) # sapply(data.staid.x.y.z, function(i) index(names(data))[names(data) == i])
      }
      data.staid.x.y.z <- c(data.staid.x.y.z[1],x.y.loc,data.staid.x.y.z[4]) # NA)
    }
    if (!is.na(st_crs(data))) {
      s.crs <- st_crs(data)
    }
  } else if (any(class(data) == "SpatRaster")) { # new data only
    data.df = crds(data, df=T, na.rm=F)
    # if (is.null(data.staid.x.y.z)) {
    data.df$staid <- 1:nrow(data.df)
    data.staid.x.y.z <- c(3,1,2,NA)
    # } else {
    #   if (!is.numeric(data.staid.x.y.z)) {
    #     data.staid.x.y.z <- match(data.staid.x.y.z, names(data.df)) # sapply(data.staid.x.y.z, function(i) index(names(data))[names(data) == i])
    #   }
    # }
    if (!is.na(st_crs(data))) {
      s.crs <- st_crs(data)
    }
  } else if (any(class(data) == "sftime")) {
    data.df = st_drop_time(data)
    geom.col.check = sapply(data.df, function (x) any(class(x)=="sfc"))
    geom.col = names(geom.col.check[geom.col.check])[1]
    coords = st_coordinates(data.df[, geom.col])
    data.df = cbind(st_drop_geometry(data.df), coords, time=st_time(data))
    x.y.loc <- match(dimnames(coords)[[2]], names(data.df))
    ############## OVO JE ZA PREDICTION, TUNE i CV
    if (is.null(data.staid.x.y.z)) {
      if ("staid" %in% names(data.df)) {
        staid.loc <- match("staid", names(data.df))
        data.staid.x.y.z <- c(staid.loc,x.y.loc,length(data.df))
      } else {
        data.df$staid = cumsum(!duplicated(coords))
        staid.loc <- length(data.df)
        data.staid.x.y.z <- c(staid.loc,x.y.loc,length(data.df)-1)
      }
    } else {
      if (!is.numeric(data.staid.x.y.z)) {
        data.staid.x.y.z <- match(data.staid.x.y.z, names(data.df)) # sapply(data.staid.x.y.z, function(i) index(names(data))[names(data) == i])
      }
      data.staid.x.y.z <- c(data.staid.x.y.z[1],x.y.loc,length(data.df))
    }
    if (!is.na(st_crs(data))) {
      s.crs <- st_crs(data)
    }
  } else if (any(class(data) == "sf")) {
    coords = st_coordinates(data)
    data.df = st_drop_geometry(data)# as.data.frame(data)
    data.df = cbind(data.df, coords)
    x.y.loc <- match(dimnames(coords)[[2]], names(data.df))
    ############## OVO JE ZA PREDICTION, TUNE i CV
    if (is.null(data.staid.x.y.z)) {
      if ("staid" %in% names(data.df)) {
        staid.loc <- match("staid", names(data.df))
      } else {
        data.df$staid <- 1:nrow(data.df)
        staid.loc <- length(data.df)
      }
      data.staid.x.y.z <- c(staid.loc,x.y.loc,NA)
    } else {
      if (!is.numeric(data.staid.x.y.z)) {
        data.staid.x.y.z <- match(data.staid.x.y.z, names(data.df)) # sapply(data.staid.x.y.z, function(i) index(names(data))[names(data) == i])
      }
      data.staid.x.y.z <- c(data.staid.x.y.z[1],x.y.loc,data.staid.x.y.z[4]) # NA)
    }
    if (!is.na(st_crs(data))) {
      s.crs <- st_crs(data)
    }
  } else if (class(data)[1] == "data.frame") {
    if (is.null(data.staid.x.y.z)) {
      stop('The argument (new)data.staid.x.y.z must not be empty if class(data) = data.frame!')
    }
    # if data.staid.x.y.z is character
    if (!is.numeric(data.staid.x.y.z)) {
      data.staid.x.y.z <- match(data.staid.x.y.z, names(data)) # sapply(data.staid.x.y.z, function(i) index(names(data))[names(data) == i])
    }
    data.df = data
  } else {
    stop('The argument data/newdata must be of sf, sftime, SpatVector, SpatRaster (as newdata only), data.frame class!') # "STSDF"
  }
  if (!is.null(obs.col) & is.numeric(obs.col)) {
      obs.col <- names(data)[obs.col]
  }
  return(list(data.df=data.df, data.staid.x.y.z=data.staid.x.y.z, s.crs=s.crs, obs.col=obs.col))
}