## Overview

Created originally as a companion of`rsplash`, `splashTools` is collection of useful functions and wrappers designed to download and pre-process datasets for ecohydrological analysis or simulation frameworks 

## Features
- aggRaster: computes daily to monthly aggregation in large grid datasets in paralell
- clipGrid: crops and masks a raster object with a polygon
- clipPointGrid: subset a sp dataframe using a raster value as conditional
- cloud2Rad: computes shortwave radiation from cloudiness
- downscaleKrige: downscales grid-type datasets using kriging
- downscaleRLM:  downscales grid-type datasets using linear models
- downscaleSolar: downscales shortwave radiation datasets using the theoretical effect terrain features
- getForcingNRCS: downloads forcing data from the US NRCS database
- getGLMchelsa: Automates the download of CHELSA downscaled datasets of the CMIP5
- getSoilGrid: downloads grid datasets, from soilgrids.org, computes a weigthed average acording to the depth intervals
- getSoilNRCS: search and download measured soil data close to US SNOTEL stations
- getSoilSite: wrapper of the REST API to query data from soilgrids.org
- krigeForcing: automates the kriging of timeseries in parallel
- readFluxdata: reads subdaily data from fluxnet or ameriflux, makes xts dataframe with daytime daily timeseries

## Installation
To install the development release  of the `splashTools` package please run: 
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "dsval/splashTools")
```


