\name{nlmodis20110712}
\alias{nlmodis20110712}
\docType{data}
\title{MODIS LST 8 day images image for the Netherlands ('2011-07-12') 
}
\description{The original 8 day MODIS LST images were also converted from Kelvin to degrees Celsius using the formula indicated in the MODIS user's manual.
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}}
\references{Wan, Z., Y. Zhang, Q. Zhang, and Z.-L. Li (2004), Quality assessment and validation of the MODIS global land surface temperature, Int. J. Remote Sens., 25(1), 261-274}
\usage{data(nlmodis20110712)}

\examples{
library(sp)
data(nlmodis20110712)
\donttest{
spplot(nlmodis20110712)
}
}
\keyword{datasets}
