#rm(list = ls())

# To be run on the cluster

# devtools::install_github("pecanproject/pecan", subdir = 'base/db')
# devtools::install_github("pecanproject/pecan", subdir = 'base/logger')
# devtools::install_github("pecanproject/pecan", subdir = 'base/qaqc')
# devtools::install_github("pecanproject/pecan", subdir = 'base/remote')
# devtools::install_github("pecanproject/pecan", subdir = 'base/settings')
# devtools::install_github("pecanproject/pecan", subdir = 'base/utils')
# devtools::install_github("pecanproject/pecan", subdir = 'base/visualization')

# devtools::install_github("pecanproject/pecan", subdir = 'modules/data.atmosphere')
# devtools::install_github("pecanproject/pecan", subdir = 'modules/emulator')
# devtools::install_github("pecanproject/pecan", subdir = 'modules/priors')
# devtools::install_github("pecanproject/pecan", subdir = 'modules/uncertainty')


#not needed (yet?)
# devtools::install_github("pecanproject/pecan", subdir = 'base/workflow')

#for workflow: 'PEcAn.data.land', 'PEcAn.settings', 'PEcAn.uncertainty'
#for data.land: 'datapack', 'dataone', 'redland'


# devtools::install_github("pecanproject/pecan", subdir = 'modules/data.land')





#PEcAn.DB
# 'datapack', 'dataone', 'redland' are not available for package 'PEcAn.data.land'
# devtools::install_github("pecanproject/pecan", subdir = 'base/utils')
# devtools::install_github("Kacoppie/pecan", subdir = 'models/ed')
# devtools::install_github("femeunier/ED2scenarios")

#PEcAn.data.atmosphere', 'PEcAn.data.land', 'PEcAn.settings' are not available for package 'PEcAn.ED2'
#'PEcAn.visualization', 'PEcAn.settings', 'PEcAn.benchmark'

#first script to be run for the SA

# Libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ED2scenarios)
library(PEcAn.ED2)
library(BayesianTools)

# Directories
ref_dir <- "/data/gent/vo/000/gvo00074/ED2kasper/run"
ed2in <- read_ed2in(file.path(ref_dir,"ED2IN_ref_short_control"))  # reference ED2IN file

rundir <- "/data/gent/vo/000/gvo00074/ED2kasper/run/SA"   # Directory for the run folders
outdir <- "/data/gent/vo/000/gvo00074/ED2kasper/out/SA"                       # Directory for the outfolders folders

# Parameters d
# Priors
#define the lower and upper limits of the parameters you want to test. for now, we assume a uniform distribution, but in a later stage we can adjust that distribution.
pft_lowers <- c(stomatal_slope = 2,
                D0 = 0.005,
                Vm0 = 5,
                Delta_Vm0 = 0,
                vm_q10 = 1.8,
                Vm_high_temp = 40,
                Vm_decay_ehigh = 0.2)

pft_uppers <- c(stomatal_slope = 16,
                D0 = 0.03,
                Vm0 = 35,
                Delta_Vm0 = 10,
                vm_q10 = 3,
                Vm_high_temp = 50,
                Vm_decay_ehigh = 0.6)

# Global thresholds
global_min <- c(Vm0 = 5)
global_max <- c(Vm0 = 35)

#create the uniform distribution
prior <- map2(pft_lowers,pft_uppers,createUniformPrior)
param_names <- names(pft_lowers)
Nparam <- length(param_names)

# Multi jobs --> to group the simulations in a single job file
Nsimuperjob = 4
isimu = 0

# Config file
PREFIX_XML <- "<?xml version=\"1.0\"?>\n<!DOCTYPE config SYSTEM \"ed.dtd\">\n"
defaults <- list_dir <- list()

# Default settings
settings <- list(model = list(revision = "git",
                              config.header = NULL),
                 pfts = list(pft = list(num = 2,
                                        ed2_pft_number = 2,
                                        name = "Early"),
                             pft = list(num = 3,
                                        ed2_pft_number = 3,
                                        name = "Mid"),
                             pft = list(num = 4,
                                        ed2_pft_number = 4,
                                        name = "Late")))

# Default config
config <- list()

#unlist simplifies the list structure and produces a vector which contains all the atomic components that occur in the initial list.
#It results in a named num.
#here, we give all the reference parameters of the pft, which are not necessarily the parameters we want to test the sensitivity of.
config[["Early"]] <- unlist(list(num = 2,
                                 Vm0 = 20,
                                 wood_psi50 = -400,
                                 stoma_psi_b = -350,
                                 phenology = 5,
                                 b1Rd = -4,
                                 mort1 = 0.0005,
                                 leaf_psi_tlp = 200))

config[["Mid"]] <- unlist(list(num = 3,
                               Vm0 = 15,
                               wood_psi50 = -450,
                               stoma_psi_b = -400,
                               phenology = 5,
                               b1Rd = -4,
                               mort1 = 0.0005,
                               leaf_psi_tlp = 200))


config[["Late"]] <- unlist(list(num = 4,
                                Vm0 = 10,
                                wood_psi50 = -500,
                                stoma_psi_b = -450,
                                phenology = 5,
                                b1Rd = -4,
                                mort1 = 0.0005,
                                leaf_psi_tlp = 200,
                                wood_water_cap = 0.2))

# Which quantiles to be run (in addition to the median)
quantiles <- c(0.025,0.16,0.25,0.75,0.84,0.975)

#################################################################
# Reference run

# Directories
run_name <- "SA_reference"
isimu = isimu + 1

#file.path: Construct the path to a file from components in a platform-independent way.
run_ref <- file.path(rundir,run_name)
out_ref <- file.path(outdir,run_name)

#check whether or not the directory exists, and if not, create it.
if(!dir.exists(run_ref)) dir.create(run_ref)
if(!dir.exists(out_ref)) dir.create(out_ref)
if(!dir.exists(file.path(out_ref,"analy"))) dir.create(file.path(out_ref,"analy"))
if(!dir.exists(file.path(out_ref,"histo"))) dir.create(file.path(out_ref,"histo"))

# ED2IN file
ed2in_scenar <- ed2in
ed2in_scenar$IEDCNFGF <- file.path(run_ref,"config.xml") #set the path to the right folder for various extra parameter files.
ed2in_scenar$FFILOUT = file.path(out_ref,"analy","analysis")
ed2in_scenar$SFILOUT = file.path(out_ref,"histo","history")

write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN")) #writes the ED2IN file based on the inputs of ed2in_scenar

# Config file
config_simu <- config

#for loop that goes over all the parameters. If the parameter is Vm0 (the maximum capacity of Rubisco to perform the carboxylase function),
#a separate loop is executed. 
for (iparam in seq(1,Nparam)){
  if (param_names[iparam] == "Vm0"){ # Pft specific case

    param_name = param_names[iparam]
    param0 <- median(prior[[iparam]]$sample(n = 100000)) #results in a value close to 20
    delta_param <- median(prior[[paste0("Delta_",param_name)]]$sample(n = 100000)) #results in a value close to 5

    params <- param0 - delta_param*c(0,1,2) #results in three values, 20, 15 and 10
    params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name]) #checks if a value lies above or below the global treshold

  } else if (param_names[iparam] == "Delta_Vm0"){ # We skip the delta, but actually we don't, if we want to skip it we can add some extra code before next()
    next()
  } else { # Default case
    param_name <- param_names[iparam]
    samples <- prior[[iparam]]$sample(n = 100000)
    params_actual <- rep(median(samples),3) #for default parameters, all pfts get the same value
  }

  for (ipft in seq(1,length(params_actual))){ 
    config_simu[[ipft]][param_name] <- params_actual[ipft]
    #in this case this is equal to three, as we have three pfts defined, but this has to be changed manually,
    #for all parameters except for the Vm0, the pfts get the same value, for Vm0, the first pft gets the highest value, and the next are that value
    #substracted by one and two times the Delta_Vm0 value (median of sampling, close to five in this case). All other parameters that were defined in
    #advance can differ between the different pfts.
  }
}

config_reference <- config_simu

# Write config file
xml <- write.config.xml.ED2(defaults = defaults,
                            settings = settings,
                            trait.values = config_simu) #writes the config file for ED2

XML::saveXML(xml, file = file.path(run_ref,"config.xml"), indent = TRUE,
             prefix = PREFIX_XML)

# job.sh
if (isimu == 1){
  isfirstjob = TRUE
  dir_joblauncher = run_ref
  list_dir[[run_name]] = run_ref
} else{
  isfirstjob = FALSE
}

# Write job file
write_joblauncher(file = file.path(dir_joblauncher,"job.sh"), # job file name
                  nodes = 1,ppn = 18,mem = 16,walltime = 3,   # job config
                  prerun = "ml purge; ml UDUNITS/2.2.26-intel-2018a R/3.4.4-intel-2018a-X11-20180131 HDF5/1.10.1-intel-2018a; ulimit -s unlimited",  # modules to load
                  CD = run_ref,  # path to the run folder
                  ed_exec = "/data/gent/vo/000/gvo00074/ED2kasper/ED2/ED/build/ed_2.1-opt-master-2bb6872", # ED2 executable name
                  ED2IN = "ED2IN", # ED2IN name
                  Rplot_function = '/data/gent/vo/000/gvo00074/ED2kasper/R/read_and_plot_ED2_Q2R_tspft.r',   # postprocessing R file
                  firstjob = isfirstjob,
                  clean = TRUE  # Do you want to clean some of the outputs
)

if (isimu == Nsimuperjob){
  isimu = 0
}

#################################################################
# Main loops for the SA runs (similar to above but one parameter is changed at the time)
# SA
for (iparam in seq(1,Nparam)){
  for (iquantile in seq(1,length(quantiles))){

    # Directories
    run_name <- paste0("SA_",param_names[iparam],"_",quantiles[iquantile])
    isimu = isimu + 1

    run_ref <- file.path(rundir,run_name)
    out_ref <- file.path(outdir,run_name)

    if(!dir.exists(run_ref)) dir.create(run_ref)
    if(!dir.exists(out_ref)) dir.create(out_ref)
    if(!dir.exists(file.path(out_ref,"analy"))) dir.create(file.path(out_ref,"analy"))
    if(!dir.exists(file.path(out_ref,"histo"))) dir.create(file.path(out_ref,"histo"))

    # ED2IN
    ed2in_scenar <- ed2in
    ed2in_scenar$IEDCNFGF <- file.path(run_ref,"config.xml")
    ed2in_scenar$FFILOUT = file.path(out_ref,"analy","analysis")
    ed2in_scenar$SFILOUT = file.path(out_ref,"histo","history")

    write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN"))
    
    
    # Config
    config_simu <- config_reference
    
    for (iparam2 in seq(1,Nparam)){
      
      param_name = param_names[iparam2]
      
      if (param_name == param_names[iparam]){ #only execute if the parameter is the same as above: we only want to change one parameter at a time.
        if (param_names[iparam] == "Vm0"){ # Pft specific case
          
          param0 <- quantile(prior[[iparam]]$sample(n = 100000),quantiles[iquantile]) #takes the respective quantile from the distribution.
          delta_param <- median(prior[[paste0("Delta_",param_name)]]$sample(n = 100000)) #here we don't take into account the quantile.
          
          params <- param0 - delta_param*c(0,1,2)
          params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name])
          
        } else if (param_names[iparam] == "Delta_Vm0"){ # We skip the delta
          # does this part override the first calculation of the Vm0 calculation?
          
          param0 <- median(prior[["Vm0"]]$sample(n = 100000))
          delta_param <- quantile(prior[[paste0("Delta_Vm0")]]$sample(n = 100000),quantiles[iquantile])
          
          params <- param0 - delta_param*c(0,1,2)
          param_name <- "Vm0"
          params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name])
          
          
        } else { # Default case
          samples <- prior[[iparam]]$sample(n = 100000)
          params_actual <- rep(quantile(samples,quantiles[iquantile]),3) #here we get the value of the parameter for each quantile value that is specified.
        }
        
        for (ipft in seq(1,length(params_actual))){
          config_simu[[ipft]][param_name] <- params_actual[ipft]
        }
      }
    }
    
    xml <- write.config.xml.ED2(defaults = defaults,
                                settings = settings,
                                trait.values = config_simu)
    
    XML::saveXML(xml, file = file.path(run_ref,"config.xml"), indent = TRUE,
                 prefix = PREFIX_XML)

    # job.sh

    if (isimu == 1){
      isfirstjob = TRUE
      dir_joblauncher = run_ref
      list_dir[[run_name]] = run_ref
    } else{
      isfirstjob = FALSE
    }

    write_joblauncher(file = file.path(dir_joblauncher,"job.sh"), # job file name
                      nodes = 1,ppn = 18,mem = 16,walltime = 3,   # job config
                      prerun = "ml purge; ml UDUNITS/2.2.26-intel-2018a R/3.4.4-intel-2018a-X11-20180131 HDF5/1.10.1-intel-2018a; ulimit -s unlimited",  # modules to load
                      CD = run_ref,  # path to the run folder
                      ed_exec = "/data/gent/vo/000/gvo00074/ED2kasper/ED2/ED/build/ed_2.1-opt-master-2bb6872", # ED2 executable name
                      ED2IN = "ED2IN", # ED2IN name
                      Rplot_function = '/data/gent/vo/000/gvo00074/ED2kasper/R/read_and_plot_ED2_Q2R_tspft.r',   # postprocessing R file
                      firstjob = isfirstjob,
                      clean = TRUE ) # Do you want to clean some of the outputs

    if (isimu >= Nsimuperjob){
      isimu = 0
    }
  }
}

dumb <- write_bash_submission(file = file.path(rundir,"all_jobs_SA.sh"),
                              list_files = list_dir,
                              job_name = "job.sh")


# To transfer the files
# scp /home/femeunier/Documents/projects/YGB/scripts/generate_SA_YGB.R hpc:/data/gent/vo/000/gvo00074/felicien/R
