\name{get_meteoeurope1km}
\alias{get_meteoeurope1km}
\alias{meteoeurope1km_files}
\title{Download the MeteoEurope1km gridded meteorological dataset from Zenodo.}
\description{\code{get_meteoeurope1km} downloads archives of the MeteoEurope1km dataset (daily gridded meteorological data for Europe at a 1-km spatial resolution, 1961-2020; Sekulić and Kilibarda, 2026) from Zenodo and verifies their MD5 checksums. \code{meteoeurope1km_files} only lists the matching archives, with their size and download URL.

The dataset is published as one Zenodo record per variable and decade with daily data (one ZIP archive per year), and one record per variable with monthly and annual aggregates and climate normals for 1961-1990 and 1991-2020. Records are looked up through the catalog in the umbrella record (\doi{10.5281/zenodo.23264155}), so the latest published version is always used.

Rasters are Cloud-Optimized GeoTIFFs in the Equi7 Europe projection (EPSG:27704), 1000 m, Int16 in tenths of the unit (°C, mm, hPa), NoData -9999. Aggregates are means for temperatures and sums for precipitation. For \code{"slp"} only daily data are available.}
\usage{
meteoeurope1km_files(var = "tmean",
                     product = "day",
                     years = NULL,
                     period = NULL)

get_meteoeurope1km(var = "tmean",
                   product = "day",
                   years = NULL,
                   period = NULL,
                   dest_dir = ".",
                   unzip = FALSE,
                   overwrite = FALSE)
}
\arguments{
\item{var}{\code{character}; Meteorological variable. Possible values are: \code{"tmean"} (mean temperature, default), \code{"tmax"} (maximum temperature), \code{"tmin"} (minimum temperature), \code{"prcp"} (total precipitation), or \code{"slp"} (sea-level pressure).}
\item{product}{\code{character}; Possible values are: \code{"day"} (daily data, one archive per year, default), \code{"mon"} (monthly aggregates, one archive per decade), \code{"ann"} (annual aggregates 1961-2020, one archive), \code{"ltm_day"} (daily climate normals), \code{"ltm_mon_ann"} (monthly and annual climate normals), or \code{"ltm_diff"} (difference of the normals, 1991-2020 minus 1961-1990).}
\item{years}{\code{numeric}; Years between 1961 and 2020. Required for \code{product = "day"} and \code{product = "mon"}; for \code{"mon"} the archives of the decades containing these years are returned.}
\item{period}{\code{character}; Climate normal period(s) for \code{product = "ltm_day"} and \code{product = "ltm_mon_ann"}: \code{"1961-1990"} and/or \code{"1991-2020"}. If \code{NULL} (default), both periods.}
\item{dest_dir}{\code{character}; Directory for the downloaded archives; created if it does not exist.}
\item{unzip}{\code{logical}; If \code{TRUE}, archives are unzipped into \code{dest_dir}.}
\item{overwrite}{\code{logical}; If \code{FALSE} (default), archives already in \code{dest_dir} with a matching MD5 checksum are not downloaded again.}
}

\value{\link[base]{data.frame} with columns \code{record} (record name), \code{doi} (DOI of the Zenodo record), \code{file} (archive name), \code{size} (bytes), \code{md5}, and \code{url}. \code{get_meteoeurope1km} also adds \code{path} (local path) and returns the data frame invisibly.}
\references{
Sekulić, A., & Kilibarda, M. (2026). MeteoEurope1km (1961-2020): daily gridded meteorological dataset for Europe at a 1-km spatial resolution. Zenodo. \doi{10.5281/zenodo.23264155}
}
\author{Aleksandar Sekulić \email{asekulic@grf.bg.ac.rs}}

\seealso{
\code{\link{get_meteo}}
}
\examples{
\donttest{
# list the daily mean temperature archives for 2001 and 2002
files <- meteoeurope1km_files(var = "tmean", product = "day", years = 2001:2002)
files[, c("file", "size", "doi")]
}
\dontrun{
# download (about 3.5 GB per year) and unzip them
get_meteoeurope1km(var = "tmean", product = "day", years = 2001:2002,
                   dest_dir = "meteoeurope1km", unzip = TRUE)
# difference of the annual precipitation normals
get_meteoeurope1km(var = "prcp", product = "ltm_diff", dest_dir = "meteoeurope1km")
}
}
