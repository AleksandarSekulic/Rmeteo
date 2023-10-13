cv.rfsi <- function (formula, # without nearest obs
                     data, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                     data.staid.x.y.z = NULL, # if data.frame
                     use.idw = FALSE,
                     # avg = FALSE,
                     # increment, # avg(nearest point dist)
                     # range, # bbox smaller(a, b) / 2
                     # direct = FALSE,
                     s.crs = NA, # brisi
                     p.crs = NA, # brisi
                     tgrid, # caret tune grid - by default random (min.node.size, mtry, no, sample.fraction, ntree, splitrule)
                     tgrid.n = 10,
                     tune.type = "LLO", # type of cv - LLO for now, after LTO, LLTO - CAST
                     k = 5, # number of folds
                     seed = 42,
                     out.folds, # if user want to create outer folds or column name
                     in.folds, # if user want to create innner folds or column name
                     acc.metric, # for tuning on inner folds
                     output.format = "data.frame", #"STFDF", # brisi - koristi kao data
                     cpus=detectCores()-1,
                     progress = 1,
                     soil3d = FALSE, # soil RFSI
                     no.obs = 'increase', # exactly
                     ...){ # fixed ranger parameters
  
  # check the input
  if (progress %in% 1:3) print('Preparing data ...')
  if ((missing(data)) | missing(formula) | missing(tgrid)) {
    stop('The arguments data, formula and tgrid must not be empty!')
  }
  formula <- as.formula(formula)
  all_vars <- all.vars(formula)
  obs.col.name <- all_vars[1]
  names_covar <- all_vars[-1]
  
  # prepare data
  data.prep <- data.prepare(data=data, data.staid.x.y.z=data.staid.x.y.z, s.crs=s.crs)
  data.df <- data.prep[["data.df"]]
  data.staid.x.y.z <- data.prep[["data.staid.x.y.z"]]
  s.crs <- data.prep[["s.crs"]]
  
  # Criteria accuracy parameter
  if (is.factor(data.df[, obs.col.name])) {
    # classification
    if (missing(acc.metric) || !(acc.metric %in% c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue"))) {
      acc.metric <- "Kappa"
      warning("Criteria accuracy parameter is missing or is not valid for classification task. Using Kappa acccuracy parameter!")
    }
  } else {
    # regression
    if (missing(acc.metric) || !(acc.metric %in% c("RMSE","NRMSE","MAE","NMAE","ME","R2","CCC"))) {
      acc.metric <- "RMSE"
      warning("Criteria accuracy parameter is missing or is not valid for regression task. Using RMSE acccuracy parameter!")
    }
  }
  
  # set out.folds and in.folds
  if (missing(out.folds)) {
    # create out.folds
    
    # check if time?
    
    # space_id <- rep(1:length(data@sp), length(time))
    # time_id <- rep(1:length(time), each = length(data@sp))
    # st_df <- cbind(space_id, time_id)
    # if (type == "LLO") {
    spacevar <- names(data.df)[data.staid.x.y.z[1]]# "space_id"
    timevar <- NA
    # TO DO LTO and LLTO
    # } else if (type == "LTO") {
    #   spacevar <- NA
    #   timevar <- "time_id"
    # } else if (type == "LLTO") {
    #   spacevar <- "space_id"
    #   timevar <- "time_id"
    # }
    indices <- CreateSpacetimeFolds(data.df, spacevar = spacevar, timevar = timevar,
                                    k=k, seed = seed)
    out.folds <- c()
    for (f in 1:length(indices$indexOut)) {
      out.folds[indices$indexOut[[f]]] <- f
    }
    data.df$out.folds <- out.folds
    out.fold.column <- "out.folds"
  } else if (class(out.folds) %in% c("numeric", "character", "integer")) {
    # outer folds
    if (length(out.folds) == 1) { # column
      out.fold.column <- out.folds
      if (class(out.fold.column) %in% c("numeric", "integer")) {
        out.fold.column <- names(data)[out.fold.column]
      } else if (inherits(out.fold.column, "character")) {
        if (!out.fold.column %in% names(data)){
          stop(paste0('Colum with name "', out.fold.column, '" does not exist in data'))
        }
      }
      print(paste0("Outer fold column: ", out.fold.column))
    } else if (length(out.folds) == nrow(data.df)){ # vector
      data.df$out.folds <- out.folds
      out.fold.column <- "out.folds"
    } else {
      stop('length(out.folds) != nrow(data).')
    }
    # inner folds
    if (missing(in.folds)) {
      in.fold.column <- NA
    } else if (class(in.folds) %in% c("numeric", "character", "integer")) {
      if (length(in.folds) == 1) { # column
        in.fold.column <- in.folds
        if (class(in.fold.column) %in% c("numeric", "integer")) {
          in.fold.column <- names(data)[in.fold.column]
        } else if (inherits(in.fold.column, "character")) {
          if (!in.fold.column %in% names(data)){
            stop(paste0('Colum with name "', in.fold.column, '" does not exist in data'))
          }
        }
        print(paste0("Inner fold column: ", out.fold.column))
      } else if (length(in.folds) == nrow(data.df)){ # vector
        data.df$in.folds <- in.folds
        in.fold.column <- "in.folds"
      } else {
        stop('length(in.folds) != nrow(data).')
      }
    } else {
      stop('The argument in.folds must numeric, integer or character.')
    }
    
  } else {
    stop('The argument out.folds must numeric, integer or character.')
  }
  
  # do CV
  if (progress %in% 1:3) print('Doing CV ...')
  pred <- c()
  # for val_fold in outer folds
  for (val_fold in sort(unique(data.df[, out.fold.column]))) {
    if (progress %in% 1:3) print(paste('### Main fold ', val_fold, " ###", sep=""))
    # tune RFSI model
    dev.df <- data.df[data.df[, out.fold.column] != val_fold, ]
    val.df <- data.df[data.df[, out.fold.column] == val_fold, ]
    if (progress %in% 1:3) print('Tuning RFSI model ...')
    # tune.rfsi
    tuned_model <- tune.rfsi(formula, # without nearest obs
                             data=dev.df, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                             data.staid.x.y.z = data.staid.x.y.z, # if data.frame
                             # n.obs=n.obs, # nearest obs
                             # time.nmax, # use all if not specified
                             s.crs=s.crs,
                             p.crs=p.crs,
                             # use.tps=use.tps,
                             # tps.df=tps.df,
                             use.idw=use.idw,
                             # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                             tgrid=tgrid, # caret tune grid - by default random (min.node.size, mtry, no, sample.fraction, ntree, splitrule)
                             tgrid.n=tgrid.n,
                             tune.type = tune.type, # type of cv - LLO for now, after LTO, LLTO - CAST
                             k = k, # number of folds
                             folds=in.folds, # if user wants to create folds or column name
                             # fold.column=in.fold.column, # by which column
                             acc.metric = acc.metric,
                             cpus=cpus, # for near.obs
                             progress=ifelse(progress==2, T, F),
                             fit.final.model=TRUE,
                             soil3d = soil3d, # soil RFSI
                             no.obs = no.obs, # exactly
                             # quantreg = quantreg,
                             # seed = seed,
                             ...)
    # validate
    if (progress %in% 1:3) print('Validation ...')
    fold_prediction <- pred.rfsi(tuned_model$final.model, # RFSI model iz rfsi ili tune rfsi funkcije
                                 data=dev.df, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                 obs.col=obs.col.name,
                                 data.staid.x.y.z = data.staid.x.y.z, # if data.frame
                                 newdata=val.df, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                 newdata.staid.x.y.z = data.staid.x.y.z, # if data.frame
                                 output.format = "data.frame",
                                 # n.obs=10, # nearest obs 3 vidi iz modela
                                 # time.nmax, # use all if not specified
                                 s.crs=s.crs,
                                 newdata.s.crs = s.crs,
                                 p.crs=p.crs,
                                 # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                                 cpus=cpus,
                                 progress=ifelse(progress==3, T, F),
                                 soil3d = soil3d, # soil RFSI
                                 depth.range = tuned_model$tuned.parameters$depth.range, # in units of depth
                                 no.obs = no.obs, # exactly
                                 ...)
    
    if (is.na(data.staid.x.y.z[4])) {
      fold_prediction <- fold_prediction[, c(names(val.df)[data.staid.x.y.z[1:3]], "pred")]
      val.df <- val.df[, c(names(val.df)[data.staid.x.y.z[1:3]], out.fold.column, obs.col.name)]
    } else {
      fold_prediction <- fold_prediction[, c(names(val.df)[data.staid.x.y.z], "pred")]
      val.df <- val.df[, c(names(val.df)[data.staid.x.y.z], out.fold.column, obs.col.name)]
    }
    val.df <- join(val.df, fold_prediction)
    pred <- rbind(pred, val.df)
    if (progress %in% 1:3) print(paste(acc.metric, ": ", acc.metric.fun(val.df[, obs.col.name], val.df$pred, acc.metric), sep=""))
  }
  names(pred)[names(pred) == obs.col.name] <- "obs"
  pred <- pred[complete.cases(pred), ]
  
  # return
  if (output.format == "sftime") {
    sf <- pred
    names(sf)[4] <- "time"
    sf <- st_as_sf(sf, coords = names(sf)[2:3], crs = s.crs, agr = "constant")
    # sf$time <- Sys.time()
    sftime <- st_sftime(sf)
    if (progress %in% 1:3) print("Done!")
    return(sftime)
  } else if (output.format == "sf") {
    sf <- st_as_sf(pred, coords = names(pred)[2:3], crs = s.crs, agr = "constant")
    if (progress %in% 1:3) print("Done!")
    return(sf)
  } else if (output.format == "SpatVector") {
    if (!is.na(s.crs)) {s.crs <- s.crs$wkt}
    sv <- terra::vect(pred, geom = names(pred)[2:3], crs = s.crs)
    if (progress %in% 1:3) print("Done!")
    return(sv)
  } else { #  (output.format == "df")
    if (progress %in% 1:3) print("Done!")
    return(pred)
  }
}

# list - final RF model, parameters - n.obs, mtry, min.node.size, splitrule