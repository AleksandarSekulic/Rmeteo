near.obs <- function(
  locations,
  locations.x.y = c(1,2),
  observations,
  observations.x.y = c(1,2),
  obs.col = 3,
  n.obs = 10,
  rm.dupl = TRUE,
  avg = FALSE,
  increment = 10000, # avg(nearest point dist)
  range = 50000, # bbox smaller(a, b) / 2
  quadrant = FALSE,
  idw=FALSE,
  idw.p=2
)
{
  ### PREPARE DATA ###
  if (any(class(locations) %in% c("SpatRaster", "SpatVector"))) {
    locations <- crds(locations)
  } else if (any(class(locations) == "sf")) {
    locations <- st_coordinates(locations)
  } else {
    locations <- locations[, locations.x.y]
  }
  if (any(class(observations) == "SpatVector")) {
    variable <- observations[[obs.col]] # as.data.frame(observations)[, obs.col]
    observations <- crds(observations)
  } else if (any(class(observations) == "sf")) {
    variable <- observations[, obs.col, drop=T]
    observations <- st_coordinates(observations)
  } else {
    variable <- observations[, obs.col]
    observations <- observations[, observations.x.y]
  }
  
  if (nrow(observations) < (n.obs+1)) {
    # return NA
    stop("Number of observations must be > n.obs")
    # nl_df <- matrix(NA, nrow = nrow(locations), ncol = (2*n.obs))
  } else {
    if (rm.dupl){
      knn1 <- nabor::knn(observations, locations, k=n.obs+1)
      knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.idx[round(knn1$nn.dists[, 1]) == 0, -1]
      knn1$nn.idx <- as.matrix(knn1$nn.idx[, -(n.obs+1)])
      knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, 1:n.obs] <- knn1$nn.dists[round(knn1$nn.dists[, 1]) == 0, -1]
      knn1$nn.dists <- as.matrix(knn1$nn.dists[, -(n.obs+1)])
    } else {
      knn1 <- nabor::knn(observations, locations, k=n.obs)
    }
    if(any(class(knn1$nn.idx)!='integer')) {
      near_o1 <- apply(knn1$nn.idx, 2, function(x) {variable[x]})
      if (nrow(locations) == 1) {
        near_o1 <- matrix(near_o1, nrow = 1)
        knn1$nn.dists <- matrix(knn1$nn.dists, nrow = 1)
      }
      near_o1 <- cbind(near_o1)
      nl_df <- cbind(knn1$nn.dists, near_o1)
    } else {
      near_o1 <- variable[knn1$nn.idx]
      nl_df <- t(c(knn1$nn.dists, near_o1))
    }
  }
  name1 <- c()
  name2 <- c()
  for (i in 1:n.obs) {
    name1 <- c(name1, paste("dist", i, sep = ""))
    name2 <- c(name2, paste("obs", i, sep = ""))
  }
  all_names <- c(name1, name2)
  
  ### IDW ###
  if (idw) {
    for (ip in idw.p) {
      if (nrow(observations) < (n.obs+1)) {
        nl_df <- cbind(nl_df, rep(NA, nrow(locations)))
      } else {
        if (any(class(knn1$nn.dists) != "numeric")) {
          if (nrow(locations) == 1) {
            idw.w <- apply(knn1$nn.dists, 1, function(x) (1/(x)^ip) / sum(1/(x)^ip) )
          } else {
            idw.w <- t(apply(knn1$nn.dists, 1, function(x) (1/(x)^ip) / sum(1/(x)^ip) ))
          }
          # wi = 1/d^2 / sum(1/d^2)
          # pred = sum(wi*obs)
          nl_df <- cbind(nl_df, apply(idw.w*near_o1, 1, sum))
        } else {
          idw.w <- (1/(knn1$nn.dists)^ip) / sum(1/(knn1$nn.dists)^ip) 
          nl_df <- cbind(nl_df, sum(idw.w*near_o1))
        }
      }
      all_names <- c(all_names, paste("idw", ip, sep="_"))
    }
  }
  
  ### AVG - BAND ###
  if(avg) {
    knn2 <- nabor::knn(observations, locations, k = nrow(observations), radius = range)
    if (rm.dupl){
      knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.idx <- as.matrix(knn2$nn.idx[, -(nrow(observations))])
      knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.dists <- as.matrix(knn2$nn.dists[, -(nrow(observations))])
    } else {
      knn1 <- nabor::knn(observations, locations, k = nrow(observations), radius = range)
    }
    knn2$nn.idx <- ifelse(knn2$nn.idx==0,NA,knn2$nn.idx)
    near_o2 <- apply(knn2$nn.idx, 2, function(x) {variable[x]})
    near_o2 <- cbind(near_o2)
    nc <- ceiling(range/increment)
    avg_st <- matrix(nrow = nrow(locations), ncol = nc)
    for (inc in 1:nc) {
      avg_st[, inc] <- apply(ifelse((knn2$nn.dists <= inc * increment) & (knn2$nn.dists >= (inc-1) * increment), near_o2, NA), 1, "mean", na.rm = T)
      # avg_st[, inc] <- apply(ifelse((knn2$nn.dists <= inc * increment), near_o2, NA), 1, "mean", na.rm = T)
    }
    avg_st[, nc] <- ifelse(is.na(avg_st[, nc]), mean(variable), avg_st[, nc])
    for (inc in (nc-1):1){
      avg_st[, inc] <- ifelse(is.na(avg_st[, inc]), avg_st[, (inc+1)], avg_st[, inc])
    }
    nl_df <- cbind(nl_df, avg_st)
    name3 <- paste("avg", trimws(format(seq(increment, range + increment-1, increment), scientific = F)), sep="")
    all_names <- c(all_names, name3)
  }
  
  ### Quadrants ###
  if(quadrant) {
    knn2 <- nabor::knn(observations, locations, k = nrow(observations))
    if (rm.dupl){
      knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.idx[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.idx <- as.matrix(knn2$nn.idx[, -(nrow(observations))])
      knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, 1:(nrow(observations)-1)] <- knn2$nn.dists[round(knn2$nn.dists[, 1]) == 0, -1]
      knn2$nn.dists <- as.matrix(knn2$nn.dists[, -(nrow(observations))])
    }
    near_o2 <- apply(knn2$nn.idx, 2, function(x) {variable[x]})
    near_o2 <- cbind(near_o2)
    
    quadrant_obs <- matrix(nrow = nrow(locations), ncol = 4)
    quadrant_obs[, 1] <- sapply(1:nrow(locations),
                              function(x) near_o2[x, (observations[knn2$nn.idx[x, ], 1]>locations[x, 1] & observations[knn2$nn.idx[x, ], 2]>locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_obs[, 2] <- sapply(1:nrow(locations),
                              function(x) near_o2[x, (observations[knn2$nn.idx[x, ], 1]>locations[x, 1] & observations[knn2$nn.idx[x, ], 2]<locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_obs[, 3] <- sapply(1:nrow(locations),
                              function(x) near_o2[x, (observations[knn2$nn.idx[x, ], 1]<locations[x, 1] & observations[knn2$nn.idx[x, ], 2]<locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_obs[, 4] <- sapply(1:nrow(locations),
                              function(x) near_o2[x, (observations[knn2$nn.idx[x, ], 1]<locations[x, 1] & observations[knn2$nn.idx[x, ], 2]>locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_obs <- ifelse(is.na(quadrant_obs), 999999, quadrant_obs)
    nl_df <- cbind(nl_df, quadrant_obs)
    name3 <- paste("dir", seq(1, 4, 1), sep="")
    all_names <- c(all_names, name3)
    
    quadrant_dist <- matrix(nrow = nrow(locations), ncol = 4)
    quadrant_dist[, 1] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], 1]>locations[x, 1] & observations[knn2$nn.idx[x, ], 2]>locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_dist[, 2] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], 1]>locations[x, 1] & observations[knn2$nn.idx[x, ], 2]<locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_dist[, 3] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], 1]<locations[x, 1] & observations[knn2$nn.idx[x, ], 2]<locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_dist[, 4] <- sapply(1:nrow(locations),
                               function(x) knn2$nn.dists[x, (observations[knn2$nn.idx[x, ], 1]<locations[x, 1] & observations[knn2$nn.idx[x, ], 2]>locations[x, 2])][1]) # knn2$nn.dists[x, observations[knn2$nn.idx, ]]
    quadrant_dist <- ifelse(is.na(quadrant_dist), 999999, quadrant_dist)
    nl_df <- cbind(nl_df, quadrant_dist)
    name3 <- paste("dir_dist", seq(1, 4, 1), sep="")
    all_names <- c(all_names, name3)
  }
  nl_df <- as.data.frame(nl_df)
  names(nl_df) <- all_names
  return(nl_df)
}









