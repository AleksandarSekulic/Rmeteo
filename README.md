# meteo
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/meteo)](https://cran.r-project.org/package=meteo)
[![cran checks](https://cranchecks.info/badges/worst/meteo)](https://cran.r-project.org/web/checks/check_results_meteo.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/meteo?color=brightgreen)](http://www.r-pkg.org/pkg/meteo)

## Overview
meteo is an R package for spatio-temporal geostatistical mapping of meteorological data.
Main functions:
* `pred.strk` - spatio-temporal regression kriging prediction (Kilibarda et al. 2014)
* `cv.strk` - k-fold cross-validation for spatio-temporal regression kriging (Kilibarda et al. 2014)
* `rfsi` - Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020)
* `tune.rfsi` - tuning of RFSI model (Sekulić et al. 2020)
* `pred.rfsi` - RFSI prediction (Sekulić et al. 2020)
* `cv.rfsi` - nested k-fold cross-validation for RFSI (Sekulić et al. 2020)


## Repositories
* [R-forge](http://meteo.r-forge.r-project.org/)
* [Github](https://github.com/AleksandarSekulic/Rmeteo)
* [CRAN](https://cran.r-project.org/package=meteo)

*Note that the latest version are on R-forge and Github repository. CRAN repository will updated only with stabile versions.*

## Citation

* Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, https://doi.org/10.1002/2013JD020803.
* Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).

