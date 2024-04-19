\name{get_meteo}
\alias{get_meteo}
\title{Get daily, monthly, or annual; aggregated or long-term means meteorological data for specific location(s) and date(s).}
\description{The function gives back daily, monthly, or annual and aggregated or long-term means metorological data for specific location(s) and date(s) for entire World for 1961-2020 period from DailyMeteo portal (https://app.dailymeteo.com/).}
\usage{
get_meteo(loc,
          var = "tmean",
          agg_level = "agg",
          time_scale = "day",
          time,
          from,
          to,
          api_key)
}
\arguments{
\item{loc}{\code{sf}, \code{SpatVector}, \link[base]{data.frame}, \code{matrix}, \code{numeric} or \code{integer}; Location(s) in WGS84 (EPSG:4326). If \link[base]{data.frame} or \code{matrix} columns are lon/lat. See \code{\link{get_coordinates}}.}
\item{var}{\code{character}; Meteorological variable. Possible values are: \code{"tmean"} (mean temperature, default), \code{"tmax"} (maximum temperature), \code{"tmin"} (minimum temperature), \code{"prcp"} (total precipitation), or \code{"slp"} (sea-level pressure).}
\item{agg_level}{\code{character}; Aggregation level. Possible values are: \code{"agg"} (aggregated, default) or \code{"ltm"} (long-term mean).}
\item{time_scale}{\code{character}; Time scale. Possible values are: \code{"day"} (daily, default), \code{"mon"} (monthly), or \code{"ann"} (annual).}
\item{time}{\code{character}; Vector of time references. If specified, ignore \code{from} and \code{to}. If \code{time_scale="day"} then format is \code{'YYYY-MM-DD'}. If \code{time_scale="mon"} then the format is \code{'YYYY-MM'}. If \code{time_scale="ann"} then format is \code{'YYYY'} (e.g. if \code{time_scale="day"} then \code{c('2020-01-01,2020-01-02')}; if \code{time_scale="mon"} then \code{c('2020-01,2020-02')}; if \code{time_scale="ann"} then \code{c('2020,2019')}).}
\item{from}{\code{character}; Time reference from. Format is the same as for the parameter \code{time}}
\item{to}{\code{character}; Time reference to. Format is the same as for the parameter \code{time}}
\item{api_key}{\code{character}; API key. Check DailyMeteo portal (https://app.dailymeteo.com/)}
}

\value{\link[base]{data.frame} object with columns \code{loc} (location index from \code{loc} argument), \code{timestamp} (time reference), and \code{value} (meteorological value).}
\author{Aleksandar Sekulić \email{asekulic@grf.bg.ac.rs}}

\seealso{
\code{\link{get_coordinates}}
}
\examples{
\dontrun{
loc <- get_coordinates("Belgrade")
loc <- rbind(loc, get_coordinates("Kopaonik"))
loc
api_key <- "" # get API key from DailyMeteo portal (https://app.dailymeteo.com/)
result <- get_meteo(loc = loc,
                   var = "tmean",
                   agg_level = "agg",
                   time_scale = "day",
                   from = "2020-01-01",
                   to = "2020-01-02", # or use time = c("2020-01-01", "2020-01-02"),
                   api_key = api_key)
# result
#   loc  timestamp value
# 1   1 2020-01-01   0.7
# 2   1 2020-01-02   1.0
# 3   2 2020-01-01  -9.2
# 4   2 2020-01-02  -8.6
# 5   3 2020-01-01  -9.2
# 6   3 2020-01-02  -8.6
}
}

