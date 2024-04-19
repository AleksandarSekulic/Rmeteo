\name{dwdsp}
\docType{data}
\encoding{utf8}
\alias{dwdsp}
\title{Daily mean wind speed in m/s for July 2011}
\description{Sample data set showing values of merged daily mean wind speed measurements from the Global Surface Summary of Day data (GSOD) with European Climate Assessment & Data set (ECA&D) for the month July 2011.}
\usage{data(dwdsp)}
\format{
The \code{dwdsp} contains the following columns:
  \describe{
  \item{\code{staid}}{character; station ID from GSOD or ECA&D data set}
  \item{\code{time}}{Date; day of the measurement}
  \item{\code{wdsp}}{numeric; daily mean wind speed in m/s}
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
data(dwdsp)
str(dwdsp)
}
\keyword{datasets}