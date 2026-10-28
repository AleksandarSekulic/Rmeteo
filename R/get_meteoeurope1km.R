.meteoeurope1km_concept <- "23264155" # concept record of the MeteoEurope1km umbrella; always resolves to the latest version

.meteoeurope1km_vars <- c("tmax", "tmin", "tmean", "prcp", "slp")
.meteoeurope1km_products <- c("day", "mon", "ann", "ltm_day", "ltm_mon_ann", "ltm_diff")
.meteoeurope1km_periods <- c("1961-1990", "1991-2020")

.zenodo_record <- function(id) {
  fromJSON(paste0("https://zenodo.org/api/records/", id))
}

.meteoeurope1km_catalog <- function() {
  umbrella <- .zenodo_record(.meteoeurope1km_concept)
  files <- umbrella$files
  url <- files$links$self[files$key == "catalog.csv"]
  if (length(url) != 1) {stop("catalog.csv not found in the MeteoEurope1km umbrella record")}
  catalog <- utils::read.csv(url, stringsAsFactors = FALSE, encoding = "UTF-8")
  catalog$id <- sub(".*zenodo\\.", "", catalog$doi)
  catalog
}

meteoeurope1km_files <- function(var = "tmean", product = "day", years = NULL, period = NULL) {
  var <- match.arg(var, .meteoeurope1km_vars)
  product <- match.arg(product, .meteoeurope1km_products)
  if (var == "slp" && product != "day") {stop("only daily data (product = \"day\") are available for slp")}
  if (product %in% c("day", "mon")) {
    if (is.null(years)) {stop("years parameter is missing")}
    years <- sort(unique(as.integer(years)))
    if (any(is.na(years) | years < 1961 | years > 2020)) {stop("years must be between 1961 and 2020")}
  }
  if (product %in% c("ltm_day", "ltm_mon_ann")) {
    period <- if (is.null(period)) .meteoeurope1km_periods else match.arg(period, .meteoeurope1km_periods, several.ok = TRUE)
  }
  decade <- function(y) {d <- 1961 + 10 * ((y - 1961) %/% 10); paste0(d, "_", d + 9)}
  per <- function(p) sub("-", "", p)
  if (product == "day") {
    record <- paste0(var, "_", decade(years))
    key <- paste0(var, "_", years, "_equi7.zip")
  } else {
    key <- switch(product,
                  mon = paste0(var, "_mon_", unique(decade(years)), "_equi7.zip"),
                  ann = paste0(var, "_ann_1961_2020_equi7.zip"),
                  ltm_day = paste0(var, "_ltm_day_", per(period), "_equi7.zip"),
                  ltm_mon_ann = paste0(var, "_ltm_mon_ann_", per(period), "_equi7.zip"),
                  ltm_diff = paste0(var, "_ltm_diff_equi7.zip"))
    record <- rep(paste0(var, "_agg_ltm_1961_2020"), length(key))
  }
  catalog <- .meteoeurope1km_catalog()
  out <- NULL
  for (rec in unique(record)) {
    id <- catalog$id[catalog$record == rec]
    if (length(id) != 1) {stop(paste0("record ", rec, " not found in the MeteoEurope1km catalog"))}
    files <- .zenodo_record(id)$files
    want <- key[record == rec]
    m <- match(want, files$key)
    if (any(is.na(m))) {stop(paste0("file(s) not found in record ", rec, ": ", paste(want[is.na(m)], collapse = ", ")))}
    out <- rbind(out, data.frame(record = rec,
                                 doi = catalog$doi[catalog$record == rec],
                                 file = files$key[m],
                                 size = files$size[m],
                                 md5 = sub("^md5:", "", files$checksum[m]),
                                 url = files$links$self[m],
                                 stringsAsFactors = FALSE))
  }
  out
}

get_meteoeurope1km <- function(var = "tmean", product = "day", years = NULL, period = NULL,
                               dest_dir = ".", unzip = FALSE, overwrite = FALSE) {
  files <- meteoeurope1km_files(var = var, product = product, years = years, period = period)
  if (!dir.exists(dest_dir)) {dir.create(dest_dir, recursive = TRUE)}
  # daily archives are 1-4 GB; the default 60 s timeout is far too short
  old <- options(timeout = max(getOption("timeout"), 3 * 3600))
  on.exit(options(old))
  files$path <- file.path(dest_dir, files$file)
  for (i in seq_len(nrow(files))) {
    path <- files$path[i]
    if (file.exists(path) && !overwrite && unname(tools::md5sum(path)) == files$md5[i]) {
      message(files$file[i], ": already downloaded")
    } else {
      message(files$file[i], ": downloading ", round(files$size[i] / 1e6), " MB")
      tmp <- paste0(path, ".part")
      status <- utils::download.file(files$url[i], tmp, mode = "wb", quiet = TRUE)
      if (status != 0) {unlink(tmp); stop("download failed: ", files$file[i])}
      if (unname(tools::md5sum(tmp)) != files$md5[i]) {unlink(tmp); stop("MD5 checksum mismatch: ", files$file[i])}
      file.rename(tmp, path)
    }
    if (unzip) {utils::unzip(path, exdir = dest_dir)}
  }
  invisible(files)
}
