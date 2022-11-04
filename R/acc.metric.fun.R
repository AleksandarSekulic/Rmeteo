acc.metric.fun <-  function(obs, pred, acc.m) {
  # classification
  if (acc.m %in% c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue")) {
    cm <- confusionMatrix(pred, obs)
    return(cm$overall[acc.m])
  # regression
  } else if (acc.m == "R2") {
    return(1 - (t(obs - pred) %*% (obs - pred)) / (t(obs - mean(obs)) %*% (obs - mean(obs))))
  } else if (acc.m == "MAE") {
    return(mean(abs(obs - pred), na.rm=TRUE)/mean(obs, na.rm=TRUE))
  } else if (acc.m == "NMAE") {
    return(sqrt(mean((obs - pred)^2, na.rm=TRUE)))
  } else if (acc.m == "CCC") {
    return(DescTools::CCC(obs, pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c)
  }  else if (acc.m == "ME") {
    return(mean((obs - pred), na.rm=TRUE))
  } else if (acc.m == "RMSE") {
    return(sqrt(mean((obs - pred)^2, na.rm=TRUE)))
  } else if (acc.m == "NRMSE") {
    return(sqrt(mean((obs - pred)^2, na.rm=TRUE))/mean(obs, na.rm=TRUE))
  } else {
    stop(paste("Accuracy paremeter ", acc.m, " is not valid!", sep=""))
  }
  
}