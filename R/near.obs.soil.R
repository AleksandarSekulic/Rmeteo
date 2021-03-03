near.obs.soil <- function(
  locations,
  locations.x.y.md = c(1,2,3),
  observations,
  observations.x.y.md = c(1,2,3),
  zcol = 4,
  n.obs = 5,
  depth.range = 0.1, # in units of depth
  no.obs = 'increase' # exactly
)
{
  ### input data preparation
  if (class(locations) == "SpatialPoints" ||
      class(locations) == "SpatialPointsDataFrame" ||
      class(locations) == "SpatialPixelsDataFrame") {
    mid.depth <- locations[, locations.x.y.md[3]]
    locations <- coordinates(locations)
    locations <- cbind(locations, mid.depth)
  } else {
    locations <- locations[, locations.x.y.md]
  }
  if (class(observations) == "SpatialPoints" || class(observations) == "SpatialPointsDataFrame") {
    mid.depth <- observations[, locations.x.y.md[3]]
    variable <- observations[[zcol]]
    observations <- coordinates(observations)
    observations <- cbind(observations, mid.depth, variable)
  } else {
    observations <- observations[, c(observations.x.y.md, zcol)]
  }
  
  ############# probably the best solution - to find the obs in the profile in which horizon is the midpoint
  ############### i.e. the midpoint is in the lower-upper range of the obs in the profile ###
  #### if there is no obs -> error("")
  # or
  #### use only horizons that have obs in the specified depth.range!!!
  
  near_o1 <- matrix(NA, nrow = nrow(locations), ncol = n.obs)
  nn.dists <- matrix(NA, nrow = nrow(locations), ncol = n.obs)
  ### loop through locations
  # add foreach in parallel - add
  for (l in 1:nrow(locations)) {
    loc <- locations[l, ]
    ### remove observations at the same location as loc
    obs.dupl <- observations[which(observations[, 1] != loc[, 1] & observations[, 2] != loc[, 2]), ]
    ### find observations at location depth +-depth.range
    obs.depth.range <- obs.dupl[obs.dupl[, observations.x.y.md[3]] >= (loc[, locations.x.y.md[3]] - depth.range) &
                                obs.dupl[, observations.x.y.md[3]] <= (loc[, locations.x.y.md[3]] + depth.range), ]
    
    ### find n nearest observations from +-depth.range
    if (nrow(obs.depth.range) < n.obs) {
      ### return error
      if (no.obs == 'exactly') {
        stop(paste(paste('There are only ', nrow(obs.depth.range), ' observations (n.obs = ', n.obs, ') within +-', depth.range, ' from location: ', sep=""),
                   paste('x = ', loc[, 1], sep=""),
                   paste('y = ', loc[, 2], sep=""),
                   paste('depth = ', loc[, 3], sep=""),
                   paste('Please, increase the depth.range.', sep=""), sep="\n"))
      } else { # increase
        multi <- 2
        while (nrow(obs.depth.range) < n.obs) {
          depth.range2 <- depth.range*multi
          ### find observations at location depth +-depth.range
          obs.depth.range <- obs.dupl[obs.dupl[, observations.x.y.md[3]] >= (loc[, locations.x.y.md[3]] - depth.range2) &
                                      obs.dupl[, observations.x.y.md[3]] <= (loc[, locations.x.y.md[3]] + depth.range2), ]
          multi <- multi + 1
        }
        warning(paste('The depth.range for location:', 
                      paste('x = ', loc[, 1], sep=""),
                      paste('y = ', loc[, 2], sep=""),
                      paste('depth = ', loc[, 3], sep=""),
                      paste('was increased to ', depth.range2, ' because there were no ', n.obs, ' nearest observations for +-', depth.range, ' depth.range.', sep=""),
                      sep="\n"))
      }
    }
    ### find n.obs nearest observations and distances to them
    knn1 <- nabor::knn(obs.depth.range[, 1:2], loc[, 1:2], k=n.obs)
    near_o1[l, ] <- obs.depth.range[knn1$nn.idx, 4]
    nn.dists[l, ] <- knn1$nn.dists
  }
  
  nl_df <- cbind(nn.dists, near_o1) # rbind
  
  name1 <- c()
  name2 <- c()
  for (i in 1:n.obs) {
    name1 <- c(name1, paste("dist", i, sep = ""))
    name2 <- c(name2, paste("obs", i, sep = ""))
  }
  all_names <- c(name1, name2)
  
  nl_df <- as.data.frame(nl_df)
  names(nl_df) <- all_names
  
  return(nl_df)

}








