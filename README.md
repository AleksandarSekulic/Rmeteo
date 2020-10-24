# meteo: RFSI and spacetime geostatistical interpolations for meteorological and other enviromental varialbles
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/meteo)](https://cran.r-project.org/package=meteo)
[![cran checks](https://cranchecks.info/badges/worst/meteo)](https://cran.r-project.org/web/checks/check_results_meteo.html)
[![Downloads](http://cranlogs.r-pkg.org/badges/meteo?color=brightgreen)](http://www.r-pkg.org/pkg/meteo)

## Installing
Install development versions (the most recent version) from github with
```
library(devtools)
install_github("https://github.com/AleksandarSekulic/Rmeteo")
```
or from R-forge
```
install.packages("meteo", repos="http://R-Forge.R-project.org")
```

## Overview
meteo is an R package for RFSI ang spacetime geostatistical interpolations for meteorological and other enviromental varialbles.
Main functions:
* `pred.strk` - spatio-temporal regression kriging prediction (Kilibarda et al. 2014)
* `cv.strk` - k-fold cross-validation for spatio-temporal regression kriging (Kilibarda et al. 2014)
* `rfsi` - Random Forest Spatial Interpolation (RFSI) model (Sekulić et al. 2020b)
* `tune.rfsi` - tuning of RFSI model (Sekulić et al. 2020b)
* `pred.rfsi` - RFSI prediction (Sekulić et al. 2020b)
* `cv.rfsi` - nested k-fold cross-validation for RFSI (Sekulić et al. 2020b)


## Repositories
* [R-forge](http://meteo.r-forge.r-project.org/)
* [Github](https://github.com/AleksandarSekulic/Rmeteo)
* [CRAN](https://cran.r-project.org/package=meteo)

*Note that the latest version is on R-forge and Github repository. CRAN repository will updated only with stabile version.*

## Citation

* Kilibarda, M., T. Hengl, G. B. M. Heuvelink, B. Graeler, E. Pebesma, M. Percec Tadic, and B. Bajat (2014), Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution, J. Geophys. Res. Atmos., 119, 2294-2313, https://doi.org/10.1002/2013JD020803.
* Kilibarda, M., Tadić, M. P., Hengl, T., Luković, J., & Bajat, B. (2015). Global geographic and feature space coverage of temperature data in the context of spatio-temporal interpolation. Spatial Statistics, 14, 22–38. https://doi.org/10.1016/j.spasta.2015.04.005
* Sekulić, A., Kilibarda, M., Protić, D., Tadić, M. P., & Bajat, B. (2020a). Spatio-temporal regression kriging model of mean daily temperature for Croatia. Theoretical and Applied Climatology, 140(1–2), 101–114. https://doi.org/10.1007/s00704-019-03077-3
* Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020b).

