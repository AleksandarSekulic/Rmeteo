\name{get_coordinates}
\alias{get_coordinates}
\title{Get lat/lon coordinates for a specific location name.}
\description{The function gives back lat/lon coordinates for a specific location name using Nominatim (https://nominatim.org/) service.}
\usage{
get_coordinates(location_name = "Belgrade")
}
\arguments{
\item{location_name}{\code{character}; Location name (e.g. "Belgrade").}
}

\value{\code{numeric} vector with latitude and longitude (WGS84) for a specific location name, in that order.}
\author{Aleksandar Sekulić \email{asekulic@grf.bg.ac.rs}}

\seealso{
\code{\link{get_meteo}}
}
\examples{
\donttest{
coords <- get_coordinates("Belgrade")
str(coords)
}
}

