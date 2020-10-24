\name{rm.dupl}
\alias{rm.dupl}
\title{Find point pairs with equal spatial coordinates from \link[spacetime]{STFDF-class} object.
}
\description{ 
This function finds point pairs with equal spatial coordinates from \link[spacetime]{STFDF-class} object and remove locations with less observations.
}
\usage{rm.dupl(obj, zcol = 1, zero.tol = 0)}
\arguments{
  \item{obj}{ \link[spacetime]{STFDF-class} object}
  \item{zcol}{ variable column name, or column number, from obj@data }
  \item{zero.tol}{distance values less than or equal to this threshold value are considered as duplicates; units are those of the coordinates for projected data or unknown projection, or km if coordinates are defined to be longitute/latitude}
}
\value{ \link[spacetime]{STFDF-class} object with removed duplicate locations. Stations with less observation is removed, if number of observation is the same for two stations the first is removed.
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}
}
\seealso{
\code{\link{tgeom2STFDF}}
\code{\link{pred.strk}}
}
\examples{
# load observation - data frame
data(dtempc) 
# load stations - data frame
data(stations)

str(dtempc)
str(stations)
# create STFDF
temp <- meteo2STFDF(dtempc,stations) # create STFDF from 2 data frames

nrow(temp@sp) # number of stations before removing dupl.

temp2 <-rm.dupl(temp)
nrow(temp2@sp) # number of stations after
}


