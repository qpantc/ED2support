rm(list = ls())

######################################################################################################
# Step 0) to do, only once

# These packages need to be installed first:
# need_packages <- c("abind", "agricolae", "akima", "beanplot", "boot", "callr", "car",
#                    "caTools", "chron", "cluster", "compiler", "data.table", "devtools",
#                    "FAdist", "fields", "gbm", "gdalUtils", "geoR", "gpclib", "grDevices",
#                    "gstat", "Hmisc", "klaR", "kriging", "leaps", "maps", "mapdata",
#                    "maptools", "MASS", "MCMCpack", "nlme", "numDeriv", "onls", "PBSmapping",
#                    "plotrix", "pls", "proto", "raster", "rgdal", "rgeos", "rlas", "robustbase",
#                    "rworldmap", "RSEIS", "R.utils","smatr","VoxR","survival")
#
# for (package in need_packages) {
#   if(!package %in% rownames(installed.packages())){install.packages(package)}
# }

#
# rhdf5 is a specific package that needs to be installed separately
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")


## ==========================================================================247
# Some default r files also need to be edited:
# load.everything.r (line 247)
# loaded.package[["hdf5"       ]] = discreet.require(hdf5        )
# must become:
# loaded.package[["rhdf5"       ]] = discreet.require(rhdf5        )

# =============================================================================
# comment some lines in load.everything.r
# # loaded.package[["gpclib"      ]] = discreet.require(gpclib      )

# #envir = as.environment("package:survival")
# #try(unlockBinding("pspline",envir),silent=TRUE)

# #envir = as.environment("package:ggplot2")
# #try(unlockBinding("theme",envir),silent=TRUE)

# #envir = as.environment("package:forecast")
# #try(unlockBinding("%>%",envir),silent=TRUE)

## ============================================================================
# functions operator.r magma.r and inferno.r: replace <<- by <-
# for operator only first part

## ===========================================================================38
# read.q.files:  replace:  soilcp     = datum$soil.prop$soilcp
#                     to:  soilcp     = rep(datum$soil.prop$soilcp,nzg)

# ED2/R-utils/read.q.files.r, line 113
# "mymont    = hdf5load(file=h5file,load=FALSE,verbosity=0,tidy=TRUE)" needs to be replaced with:
# mymont    = lapply(h5read_opt(h5file),FUN=aperm)
# names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".")



## ===========================================================================14
# ED2/R-utils/monthly.template.r, line 14
# "mymont    = hdf5load(file=h5first,load=FALSE,verbosity=0,tidy=TRUE)" needs to be replaced with:
# mymont    = lapply(h5read_opt(h5first),FUN=aperm)
# names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".")



# Then we read the outputs files with the read_and_plot_ED2.2_all_tspft function
## working dir "./EDsupport/R/"

source("../../../R/h5read_opt.r")
source("../../../R/read_and_save_ED2.2.R")
source("../../../R/read_save_plot_ED2.2.R")
ED_utils_dir = "../../../R-utils"

# /scratch/gent/vo/000/gvo00074/vsc44253/ED2.2/EDsupport/outputs/Results/BCI/analy
# /Users/tiacc/projects/ED2.2/EDsupport/outputs/Results/BCI/analy
# "../output/Results/BCI/analy"
read_and_save_ED2.2(there = '/scratch/gent/vo/000/gvo00074/vsc44253/ED2.2/EDsupport/run/BCI_histo/outputs/BCI/analy', # path to the analy outputs (Q files)
                              place = 'analysis',                                              # output name
                              yeara = '1901/01/01',                                                 # first year/month to process
                              yearz = '2011/01/01',                                                 # last year/month (+ 1) to process
                              ED2srcdir = ED_utils_dir)

# read_and_save_ED2.2(there = '/Users/tiacc/projects/ED2.2/EDsupport/outputs/Results/BCI/analy', # path to the analy outputs (Q files)
#                     place = 'analysis',                                              # output name
#                     yeara = '1901/01/01',                                                 # first year/month to process
#                     yearz = '1903/01/01',                                                 # last year/month (+ 1) to process
#                     ED2srcdir = "../R-utils")

# read_and_plot_ED2_Q2R(there = '/Users/tiacc/projects/ED2.2/EDsupport/outputs/Results/BCI/analy', # path to the analy outputs (Q files)
#                     place = 'analysis',                                              # output name
#                     yeara = '1901/01/01',                                                 # first year/month to process
#                     yearz = '1903/01/01',                                                 # last year/month (+ 1) to process
#                     ED2srcdir = "../R-utils")

# ED2srcdir is the location of the R-utils library from the github ED2 repository
# read_and_plot_ED2.2_all_tspft function generates a .RData file + figures

# We then load the output files
#load("./outputs/Japan_default.RData")

# Plot something
#matplot(datum$szpft$gpp[,12,c(2,3,4,18)],type = "l")

