\name{get_coordinates}
\alias{get_coordinates}
\title{Get lon/lat coordinates for a specific location name.}
\description{The function gives back lon/lat coordinates for a specific location name using Nominatim (https://nominatim.org/) service.}
\usage{
get_coordinates(location_name = "Belgrade")
}
\arguments{
\item{location_name}{\code{character}; Location name (e.g. "Belgrade").}
}

\value{\code{numeric} vector with lon/lat values for a specific location name.}
\author{Aleksandar Sekulić \email{asekulic@grf.bg.ac.rs}}

\seealso{
\code{\link{get_meteo}}
}
\examples{
# coords <- get_coordinates("Belgrade")
# coords
}

