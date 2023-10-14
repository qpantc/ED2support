rm(list = ls())

library(ncdf4)
library(dplyr)
library(lubridate)
library(tibble)
library(glue)
library(utils)
library(udunits2)
library(hdf5r)

## working dir "./EDsupport/R/"
source("./function/download.CRUNCEP.R")
source("./function/convert.CRUNCEP.R")
source("./function/other.functions.R")

output.folder <- "../run/BCI_histo/inputs/driver/" #"/data/gent/vo/000/gvo00074/vsc44253/R/climate.site/"
all.years <- 1981:2010
lat <- 9.25
lon <- -79.75

#######################################################################
# No change beyond this point

dir.create(output.folder,showWarnings = FALSE)

site.folder <- download.CRUNCEP(output.folder,years = all.years,lat.in = lat,lon.in = lon,verbose = FALSE)

convert.CRUNCEP(in.path = file.path(site.folder),
                in.prefix = "CRUNCEP",
                outfolder = file.path(site.folder,"ED2"),
                start_date = paste0(min(all.years),"/01/01"),
                end_date = paste0(max(all.years),"/12/31"),
                lst = 0,
                lat = lat,
                lon = lon,
                overwrite = TRUE,
                verbose = FALSE,
                leap_year = TRUE)


