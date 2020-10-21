\name{dtemp_minc}
\docType{data}
\encoding{utf8}
\alias{dtemp_minc}
\title{Minimum daily temperature in degrees Celsius for July 2011}
\description{Sample data set showing values of merged minimum daily temperature measurements from the Global Surface Summary of Day data (GSOD) with European Climate Assessment & Data set (ECA&D) for the month July 2011.}
\usage{data(dtemp_minc)}
\format{
The \code{dtemp_minc} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from GSOD or ECA&D data set}
  \item{\code{time}}{Date; day of the measurement}
  \item{\code{temp_minc}}{numeric; minimum daily temperature in degrees Celsius}
}
}
\note{The data summaries provided here are based on data exchanged under the World Meteorological Organization (WMO) World Weather Watch Program. To prepare a point map, merge with the \code{\link{stations}} table containing stations' coordinates.}
\author{Milan Kilibarda and Tomislav Hengl}
\references{
\itemize{
\item Global Surface Summary of the day data (\url{ftp://ftp.ncdc.noaa.gov/pub/data/gsod/}) 
\item European Climate Assessment & Dataset (\url{http://eca.knmi.nl/dailydata/predefinedseries.php}) 
}
}
\examples{
## load data
data(dtemp_minc)
str(dtemp_minc)
}
\keyword{datasets}