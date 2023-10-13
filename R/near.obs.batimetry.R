near.obs.batimetry <- function(
  locations,
  locations.x.y = c(1,2),
  observations,
  observations.x.y = c(1,2),
  zcol = 3,
  IHP="IHP",
  n.obs = 10,
  rm.dupl = TRUE
){
  if (inherits(locations, "SpatialPoints") ||
      inherits(locations, "SpatialPointsDataFrame") ||
      inherits(locations, "SpatialPixelsDataFrame")) {
    locations <- coordinates(locations)
  } else {
    locations <- locations[, c(locations.x.y, IHP, zcol)]
  }
  if (inherits(observations, "SpatialPoints") || inherits(observations, "SpatialPointsDataFrame")) {
    variable <- observations[[IHP]]
    observations <- coordinates(observations)
  } else {
    # variable <- observations[, IHP]
    observations <- observations[, c(observations.x.y, IHP, zcol)]
  } 
  
  locations1 <- locations
  locations2 <- observations
  locations1_sf <-sf::st_as_sf(locations1, coords = locations.x.y)
  locations2_sf <-sf::st_as_sf(locations2, coords = observations.x.y)
  
  #create matrix with euclidian distances between stations
  mat <- sf::st_distance(locations1_sf, locations2_sf, by_element=FALSE)
  locations_df <- as.data.frame(mat)%>%
  units::drop_units() # get rid of "m" 
  locations_df[locations_df==0] <- NA # distance to station itself replaced by NA
  
  # create matrix with differences in bathymetric depth between stations
  mat_bat <-outer(-locations1[[IHP]], locations2[[IHP]], "+")
  mat_bat[mat_bat==0] <- NA
  
  #compute real distances: taking into account bathymetry
  distmat <- sqrt((mat_bat^2)+(locations_df^2)) # pythagoras
  
  # creating and labeling of nearmat dataframe (contains distances to n nearest stations and the TOC value)
  counter=0
  nplus <- n.obs+1
  cols=2*n.obs
  myrange <- 1:(2*n.obs)
  myrange2 <- 1:n.obs
  nearmat <- data.frame(matrix(NA, nrow = nrow(locations1), ncol = cols)) #create matrix equal to output of near.obs

  for(i in myrange) { #name matrix with the colnames needed (obs1 to obs n, dist1 to distn)
    if (i <= n.obs){
      distcol <- paste("dist", i, sep = "")
      colnames(nearmat)[i] <- distcol
    } else {
      k <- i-n.obs
      #print(k)
      obscol <- paste("obs", k, sep = "")
      colnames(nearmat)[i] <- obscol  }
  }
  
  
  locations2$ind <- 1:nrow(locations2) # creates index column that can be used to get toc values for nearmat later
  
  #filling of nearmat dataframe
  for (r1 in 1:nrow(locations1)){
    #print(paste("r:",r1))
    
    ref <-locations1[[IHP]][r1]
    print(paste("refrence depth:",ref))
    mask <- dplyr::filter(locations2, IHP > ref) #selecting all the stations that lie "higher" than station
    #first version line above:mask <- dplyr::filter(locations2, IHP >= ref): problem with the >= it chooses
    #the station itself as well (and we get NAs from the distance matrix), and I wasn't able to get rid of the station itself in locations2 before doing filtering (
    #only possible when locations1 and lcoations2 and indices are the same). problem with this approach: if two stations have by chance the same bathymetric depth, 
    #a close station doesn't get chosen although it could be very close by. But the bigger as makes also sense fi we assume
    #during particle transport some downward settling motion and exclude the possibility of horizontal transfer from one
    #height to the same height....
    
    if (nrow(mask) >= n.obs) { # enough observations in mask in order to fill nearmat for respective station
      #print(paste("obs mask:", nrow(mask), "nn:", n.obs))
      print("deep water")
      vect<- mask$ind
      dis <- distmat[r1,vect]
      dis2 <-dis %>% 
        select(where(~!any(is.na(.))))      
      tempmat <- rbind(dis2,vect)
      tempmat2 <-data.table::data.table(t(tempmat))
      colnames(tempmat2) <- c("V1", "V2")
      ordered <- tempmat2[order(tempmat2$V1, na.last=TRUE)]
      chosen <- utils::head(ordered, n.obs) #gets first n items (smallest n distances)
      chosen2 <- data.table::data.table(chosen)
      vect2 <- chosen2$V2
      distances_insert <- chosen2$V1
      observation_insert <-locations2[, zcol][vect2] # $toc_station[vect2]
      #browser()
      nearmat[r1,1:n.obs] <- distances_insert #get distances into first n columns of near matrix
      nearmat[r1,nplus:cols] <- observation_insert # get TOC values of nearest stations
      #browser()
    }
    
    else { # chances are high station is close to coast in shallow water and there is not 
      #enough stations for the set n.obs that lie higher then the respective station
      #print(paste("obs mask:", nrow(mask), "nn:", n.obs))
      print("shallow water")
      dis <- distmat[r1,]
      dis2 <-dis %>% 
        select(where(~!any(is.na(.))))
      #browser()
      index <- 1:ncol(dis2)
      tempmat <- rbind(dis2,index)
      tempmat2 <-data.table::data.table(t(tempmat))
      colnames(tempmat2) <- c("V1", "V2")
      ordered <- tempmat2[order(tempmat2$V1, na.last=TRUE)]
      chosen <- utils::head(ordered, n.obs)
      chosen2 <- data.table::data.table(chosen)
      vect<-chosen2[["V2"]]
      nearmat[r1,1:n.obs] <- chosen2$V1 #get distances into first n columns of near matrix
      nearmat[r1,nplus:cols] <- locations2[, zcol][vect] # $toc_station[vect]
      #browser()
    }
  }
  return(nearmat)
  #}
}
