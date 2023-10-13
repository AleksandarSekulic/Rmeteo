\name{dslp}
\docType{data}
\encoding{utf8}
\alias{dslp}
\title{Mean sea level pressure in hPa for July 2011}
\description{Sample data set showing values of merged mean sea level pressure measurements from the Global Surface Summary of Day data (GSOD) with European Climate Assessment & Data set (ECA&D) for the month July 2011.}
\usage{data(dslp)}
\format{
The \code{dslp} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from GSOD or ECA&D data set}
  \item{\code{time}}{Date; day of the measurement}
  \item{\code{slp}}{numeric; mean sea level pressure amount in hPa}
}
}
\note{The data summaries provided here are based on data exchanged under the World Meteorological Organization (WMO) World Weather Watch Program. To prepare a point map, merge with the \code{\link{stations}} table containing stations' coordinates.}
\author{Milan Kilibarda and Tomislav Hengl}
\references{
\itemize{
\item Global Surface Summary of the day data (\url{ftp://ftp.ncdc.noaa.gov/pub/data/gsod/}) 
\item European Climate Assessment & Dataset (\url{https://www.ecad.eu/dailydata/predefinedseries.php}) 
}
}
\examples{
## load data
data(dslp)
str(dslp)
}
\keyword{datasets}