\name{dsndp}
\docType{data}
\encoding{utf8}
\alias{dsndp}
\title{Daily snow depth in cm for July 2011}
\description{Sample data set showing values of merged daily snow depth measurements from the Global Surface Summary of Day data (GSOD) with European Climate Assessment & Data set (ECA&D) for the month July 2011.}
\usage{data(dsndp)}
\format{
The \code{dsndp} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from GSOD or ECA&D data set}
  \item{\code{time}}{Date; day of the measurement}
  \item{\code{sndp}}{numeric; daily snow depth in cm}
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
# load data 
data(dsndp)
str(dsndp)
}
\keyword{datasets}