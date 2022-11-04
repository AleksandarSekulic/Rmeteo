tune.rfsi <- function (formula, # without nearest obs
                       data, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                       data.staid.x.y.z=  NULL, # if data.frame
                       zero.tol=0,
                       # avg = FALSE,
                       # increment, # avg(nearest point dist)
                       # range, # bbox smaller(a, b) / 2
                       # direct = FALSE,
                       use.idw = FALSE,
                       s.crs = NA,
                       p.crs = NA,
                       tgrid, # caret tune grid (min.node.size, mtry, no, sample.fraction, ntree, splitrule, idw.p, depth.range)
                       tgrid.n = 10,
                       tune.type = "LLO", # type of cv - LLO for now, after LTO, LLTO - CAST
                       k = 5, # number of folds
                       seed = 42,
                       folds, # if user want to create folds
                       fold.column, # by which column
                       acc.metric, # caret parameters - RMSE, MAE, R2, ...
                       fit.final.model = TRUE,
                       cpus=detectCores()-1,
                       progress = TRUE,
                       soil3d = FALSE, # soil RFSI
                       no.obs = 'increase', # exactly
                       ...){ # ranger parameters
  
  # check the input
  if (progress) print('Preparing data ...')
  if ((missing(data)) | missing(formula) | missing(tgrid)) {
    stop('The arguments data, formula and tgrid must not be empty!')
  }
  formula <- as.formula(formula)
  all_vars <- all.vars(formula)
  obs.col.name <- all_vars[1]
  names_covar <- all_vars[-1]
  
  # prepare data
  data.prep <- data.prepare(data = data, data.staid.x.y.z)
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
    if (missing(acc.metric) || !(acc.metric %in% c("RMSE","NRMSE", "MAE", "NMAE","ME","R2","CCC"))) {
      acc.metric <- "RMSE"
      warning("Criteria accuracy parameter is missing or is not valid for regression task. Using RMSE acccuracy parameter!")
    }
  }
  
  # set folds
  if (missing(fold.column)){
    if (missing(folds)) {
      # create folds
      
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
      folds <- c()
      for (f in 1:length(indices$indexOut)) {
        folds[indices$indexOut[[f]]] <- f
      }
    }
    data.df$folds <- folds
    fold.column <- "folds"
  }
  
  if (length(tgrid) == 0) {
    stop('The paremeter tgrid is empty!')
  }
  # check tgrid parameters, if missing add default
  params.diff <- setdiff(c("min.node.size", "num.trees", "mtry", "n.obs", "sample.fraction", "idw.p"), names(tgrid))
  if (!use.idw) {
    params.diff <- params.diff[!(params.diff %in% "idw.p")]
  }
  if (!soil3d) {
    params.diff <- params.diff[!(params.diff %in% "depth.range")]
  }
  if (!identical(params.diff, character(0))) {
    warning(paste('The paremeter(s) ', paste(params.diff, collapse = ", "), ' - missing from tgrid. Setting default values!', sep = ""))
    for (pd in params.diff) {
      if (pd == "min.node.size") {
        tgrid$min.node.size <- NULL
      } else if (pd == "num.trees") {
        tgrid$num.trees <- 500
      } else if (pd == "mtry") {
        tgrid$mtry <- NULL
      } else if (pd == "n.obs") {
        tgrid$n.obs <- 10
      } else if (pd == "sample.fraction") {
        tgrid$sample.fraction <- 1
      } else if (pd == "splirule") {
        tgrid$splirule <- NULL
      } else if (pd == "idw.p") {
        tgrid$idw.p <- 2
      } else if (pd == "depth.range") {
        tgrid$depth.range <- 0.1
      }
    }
  }

  # set combinations
  if ("mtry" %in% names(tgrid)) {
    tgrid <- tgrid[tgrid$mtry < (length(names_covar)+2*tgrid$n.obs-1), ]
  }
  if(tgrid.n > nrow(tgrid)) {
    warning(paste('The argument tgrid.n = ', tgrid.n, ' is larger than nrow(tgrid) = ', nrow(tgrid), '! -> tgrid.n is set to ', nrow(tgrid), '.', sep = ""))
    tgrid.n <- nrow(tgrid)
  }
  tgrid <- tgrid[sample(nrow(tgrid), tgrid.n),]
  tgrid <- tgrid[order(tgrid$n.obs),]
  acc.metric.vector <- rep(NA, nrow(tgrid))
  
  for (tg in 1:nrow(tgrid)) {
    comb <- tgrid[tg, ]
    if (progress) print(paste("combination ", tg, ": ", sep=""))
    if (progress) print(comb)
    
    # parameters
    min.node.size=comb$min.node.size
    num.trees=comb$num.trees
    mtry=comb$mtry
    n.obs=comb$n.obs
    sample.fraction=comb$sample.fraction
    splitrule=comb$splitrule
    idw.p=comb$idw.p
    depth.range=comb$depth.range
    
    # fold_obs <- list()
    fold_pred <- c()
    for (val_fold in sort(unique(data.df[, fold.column]))) {
      if (progress) print(paste('Fold ', val_fold, sep=""))
      # fit RF model
      dev.df <- data.df[data.df[, fold.column] != val_fold, ]
      val.df <- data.df[data.df[, fold.column] == val_fold, ]
      if (progress) print('Fitting RFSI model ...')
      fold_model <- rfsi(formula, # without nearest obs
                         data=dev.df, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                         data.staid.x.y.z = data.staid.x.y.z, # if data.frame
                         zero.tol=zero.tol,
                         n.obs=n.obs, # nearest obs
                         # time.nmax, # use all if not specified
                         s.crs=s.crs,
                         p.crs=p.crs,
                         # use.tps=use.tps,
                         # tps.df=tps.df,
                         use.idw=use.idw,
                         idw.p=idw.p,
                         # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                         cpus=cpus, # for near.obs
                         progress=progress,
                         soil3d = soil3d,
                         depth.range = depth.range,
                         no.obs = no.obs,
                         num.trees = num.trees, mtry = mtry, splitrule = splitrule,
                         min.node.size = min.node.size,  sample.fraction = sample.fraction,
                         ...)
      # predict 
      if (progress) print('Prediction ...')
      fold_prediction <- pred.rfsi(model=fold_model, # RFSI model iz rfsi ili tune rfsi funkcije
                                   data=dev.df, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                   obs.col=obs.col.name,
                                   data.staid.x.y.z = data.staid.x.y.z, # if data.frame
                                   newdata=val.df, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                                   newdata.staid.x.y.z = data.staid.x.y.z, # if data.frame
                                   output.format = "data.frame",
                                   zero.tol=zero.tol,
                                   # n.obs=10, # nearest obs 3 vidi iz modela
                                   # time.nmax, # use all if not specified
                                   s.crs=s.crs,
                                   newdata.s.crs=s.crs,
                                   p.crs=p.crs,
                                   # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                                   cpus=cpus,
                                   progress=progress,
                                   soil3d = soil3d, # soil RFSI
                                   depth.range = depth.range, # in units of depth
                                   no.obs = no.obs, # exactly
                                   ...)
      
      # fold_obs <- c(fold_obs, val.df[, obs.col.name])
      # fold_pred <- c(fold_pred, fold_prediction$pred)
      # fold_obs[[val_fold]] <- val.df[, obs.col.name]
      # fold_pred[[val_fold]] <- fold_prediction$pred
      if (is.na(data.staid.x.y.z[4])) {
        val.df <- join(val.df[, c(names(val.df)[data.staid.x.y.z[1:3]], obs.col.name)], fold_prediction[, c(names(val.df)[data.staid.x.y.z[1:3]], "pred")])
      } else {
        val.df <- join(val.df[, c(names(val.df)[data.staid.x.y.z], obs.col.name)], fold_prediction[, c(names(val.df)[data.staid.x.y.z], "pred")])
      }
      
      val.df <- val.df[, c(obs.col.name, "pred")]
      fold_pred <- rbind(fold_pred, val.df)
    }
    fold_pred <- fold_pred[complete.cases(fold_pred), ]
    acc.metric.vector[tg] <- acc.metric.fun(fold_pred[, obs.col.name], fold_pred[, "pred"], acc.metric)
    if (progress) print(paste(acc.metric, ": ", acc.metric.vector[tg], sep=""))
    if (progress) print("")
    gc(); gc()
  }
  tgrid <- cbind(tgrid, acc.metric.vector)
  names(tgrid)[length(tgrid)] <- acc.metric
  
  if (progress) print('Done!')
  
  if (acc.metric %in% c("RMSE", "NRMSE", "MAE", "NMAE", "ME", "AccuracyPValue", "McnemarPValue")) {
    dev_parameters <- tgrid[which.min(acc.metric.vector), ]
  } else {
    dev_parameters <- tgrid[which.max(acc.metric.vector), ]
  }
  if (progress) print('Final parameters: ')
  if (progress) print(dev_parameters)
  
  results <- list(combinations=tgrid, tuned.parameters=dev_parameters)
  if (fit.final.model) {
    final.model <- rfsi(formula, # without nearest obs
                        data=data.df, # data.frame(x,y,obs,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                        data.staid.x.y.z = data.staid.x.y.z, # if data.frame
                        zero.tol=zero.tol,
                        n.obs=dev_parameters$n.obs, # nearest obs
                        # time.nmax, # use all if not specified
                        s.crs=s.crs,
                        p.crs=p.crs,
                        # use.tps=use.tps,
                        # tps.df=tps.df,
                        use.idw=use.idw,
                        idw.p=dev_parameters$idw.p,
                        # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                        cpus=cpus, # for near.obs
                        progress=progress,
                        soil3d = soil3d,
                        depth.range = dev_parameters$depth.range,
                        no.obs = no.obs,
                        num.trees = dev_parameters$num.trees,
                        mtry = dev_parameters$mtry,
                        splitrule = dev_parameters$splitrule,
                        min.node.size = dev_parameters$min.node.size,
                        sample.fraction = dev_parameters$sample.fraction,
                        ...)#,
                        # quantreg = quantreg,
                        # seed = seed,
                        # ...)
    results[["final.model"]] <- final.model
  }
  return(results)
  
}

# list - final RF model, parameters - n.obs, mtry, min.node.size, splitrule