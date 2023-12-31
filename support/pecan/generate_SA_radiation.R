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

rundir <- "/data/gent/vo/000/gvo00074/ED2kasper/run/SA/Radiation"   # Directory for the run folders
outdir <- "/data/gent/vo/000/gvo00074/ED2kasper/out/SA/Radiation"                       # Directory for the outfolders folders

# Parameters d
# Priors
#define the lower and upper limits of the parameters you want to test. for now, we assume a uniform distribution, but in a later stage we can adjust that distribution.
#I leave root respiration out as I don't find the correct parameter name 
#All parameters
pft_lowers <- c(stomatal_slope = 2,
                D0 = 0.005,
                Vm0 = 5,
                vm_q10 = 1.8,
                Vm_high_temp = 40,
                Vm_decay_ehigh = 0.2,
                b1Rd = -2, #median different than the ed default
                b2Rd = 0.05, #median different than the ed default
                q = 0.5,
                root_beta = 0.0001, #median different than the ed default
                SRA = 24,
                stoma_psi_c = 1,
                #root_resp = 0.14, #correct parameter name?
                orient_factor = -0.5, 
                clumping_factor = 0.4,
                #from here, the normal distributions
                wood_Kexp = 2, #median different than the ed default
                stoma_psi_b = -160, #median different than the ed default
                b1Bl = 0.049,
                b2Bl = 1.89,
                
                #lognormal distributions
                wood_Kmax = -3.0, # my estimation is -2.803619, 
                wood_water_cap = -4.906275, #used to be 2, but then the median was wrong. This value is calculated under the assumption that the sdlog is correct in the paper. formula: median = exp^mu, and solve for mu #median different than the ed default
                leaf_water_cap = -7.195437, #Same as above for wood_water_cap #median different than the ed default
                leaf_psi_tlp = 5.42, #median different than the ed default
                
                #radiation related:
                leaf_reflect_vis = -2.63242631741772,
                leaf_reflect_nir = -0.684674273668841,
                leaf_trans_vis = -4.97609552720594,
                leaf_trans_nir = -1.423545451587,
                #weibull distribution
                root_turnover_rate = 1.6,
                # Vm0 = 40, to apply in a later stage if Félicien answers
                
                #gamma distributions
                wood_psi50 = 2.93045092,
                rho = 14.49104 #from own dataset
                )

pft_uppers <- c(stomatal_slope = 16,
                D0 = 0.03,
                Vm0 = 35,
                vm_q10 = 3,
                Vm_high_temp = 50,
                Vm_decay_ehigh = 0.6,
                b1Rd = -0.1,
                b2Rd = 0.6,
                q = 1.5,
                root_beta = 0.1, 
                SRA = 72,
                stoma_psi_c = 5,
                #root_resp = 0.42, #correct parameter name?
                orient_factor = 0.5,
                clumping_factor = 1,
                #from here, the normal distributions
                wood_Kexp = 0.5,
                stoma_psi_b = 40,
                b1Bl = 0.010,
                b2Bl = 0.2,
                #lognormal distributions
                wood_Kmax = 0.75, # my estimation is 1.22863, 
                wood_water_cap = 0.5, 
                leaf_water_cap = 0.76,
                leaf_psi_tlp = 0.53,
                
                #radiation related:
                leaf_reflect_vis = 0.333003856292686,
                leaf_reflect_nir = 0.236135184459279,
                leaf_trans_vis = 0.923579621142062,
                leaf_trans_nir = 0.355592014323824,
                #weibull
                root_turnover_rate = 1.6,
                #Vm0 = 1.35, to apply in a later stage if Félicien answers
                #gamma distributions
                wood_psi50 = 0.01884666,
                rho = 30.53080)

# Global thresholds
global_min <- c(Vm0 = 5, leaf_reflect_vis = 0, leaf_reflect_nir = 0, leaf_trans_vis = 0, leaf_trans_nir = 0)
global_max <- c(Vm0 = 35, leaf_reflect_vis = 1, leaf_reflect_nir = 1, leaf_trans_vis = 1, leaf_trans_nir = 1)

#create the uniform distribution
param_names <- names(pft_lowers)
Nparam <- length(param_names)

prior <- list()
for (iparam in 1:Nparam){
  if (param_names[iparam] %in% c("stomatal_slope", "D0", "Vm0", "vm_q10", "Vm_high_temp", "Vm_decay_ehigh", "b1Rd",
                                 "b2Rd", "q", "root_beta", "SRA", "stoma_psi_c", "orient_factor", "clumping_factor")) {
    name <- param_names[iparam]
    prior[[name]]$sampler <- function(n = 1){
      d = runif(n, min = pft_lowers[iparam], max = pft_uppers[iparam])
      return(d)
    }
    
  } else if (param_names[iparam] %in% c("wood_Kexp", "stoma_psi_b", "b1Bl", "b2Bl")) {
    
    prior[[param_names[iparam]]]$sampler <- function(n = 1){
      d = rnorm(n, mean = pft_lowers[iparam], sd = pft_uppers[iparam])
      return(d)
    }
    
  } else if (param_names[iparam] %in% c("wood_water_cap", "wood_Kmax", "leaf_water_cap", "leaf_psi_tlp", "leaf_reflect_vis",
                                        "leaf_reflect_nir", "leaf_trans_vis", "leaf_trans_nir")) { #
    
    if (param_names[iparam] %in% c("leaf_psi_tlp")){
      
      prior[[param_names[iparam]]]$sampler <- function(n = 1){
        d = - rlnorm(n, meanlog = pft_lowers[iparam], sdlog = pft_uppers[iparam])
        return(d)
      }
    }else{
      prior[[param_names[iparam]]]$sampler <- function(n = 1){
        d =  rlnorm(n, meanlog = pft_lowers[iparam], sdlog = pft_uppers[iparam])
        return(d)
      }
    }
    
  } else if (param_names[iparam] %in% c("root_turnover_rate")){ #here, Vm0 was included
    prior[[param_names[iparam]]]$sampler <- function(n = 1){
      d = rweibull(n, shape = pft_lowers[iparam], scale = pft_uppers[iparam])
      return(d)
    }
    
  } else if (param_names[iparam] %in% c("wood_psi50", "rho")){
    
    if (param_names[iparam] %in% c("wood_psi50")){
      prior[[param_names[iparam]]]$sampler <- function(n = 1){
        d = - rgamma(n, shape = pft_lowers[iparam], rate = pft_uppers[iparam])
        return(d)
      }
    } else {
      prior[[param_names[iparam]]]$sampler <- function(n = 1){
        d = rgamma(n, shape = pft_lowers[iparam], rate = pft_uppers[iparam])
        return(d)
      }
    }
  }
}

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
                                        name = "Late"),
                             pft = list(num = 17,
                                        ed2_pft_number = 17,
                                        name = "Liana")))

# Default config
config <- list()

#unlist simplifies the list structure and produces a vector which contains all the atomic components that occur in the initial list.
#It results in a named num.
#here, we give all the reference parameters of the pft, which are not necessarily the parameters we want to test the sensitivity of.

#some parameters have to be changed, as the median of the distribution is not the same as the ed default value.
config[["Early"]] <- unlist(list(num = 2,
                                 Vm0 = 20,
                                 wood_psi50 = -400,
                                 stoma_psi_b = -350,
                                 phenology = 5,
                                 b1Rd = -4,
                                 mort1 = 0.0005,
                                 leaf_psi_tlp = -200))

config[["Mid"]] <- unlist(list(num = 3,
                               Vm0 = 15,
                               wood_psi50 = -450,
                               stoma_psi_b = -400,
                               phenology = 5,
                               b1Rd = -4,
                               mort1 = 0.0005,
                               leaf_psi_tlp = -200))


config[["Late"]] <- unlist(list(num = 4,
                                Vm0 = 10,
                                wood_psi50 = -500,
                                stoma_psi_b = -450,
                                phenology = 5,
                                b1Rd = -4,
                                mort1 = 0.0005,
                                leaf_psi_tlp = -200,
                                wood_water_cap = 0.2))

#aan te vullen
config[["Liana"]] <- unlist(list(num = 17,
                                 #phenology = 5, which phenology do lianas have?
                                 mort1 = 2,
                                 stomatal_slope = 9,
                                 D0 = 0.0175, #different from standard, which is 0.0160000008
                                 Vm0 = 20, #different then standard, and because it is in a uniform distribution now.
                                 vm_q10 = 2.4, #different from standard, which is 2.4000000954
                                 Vm_high_temp = 45, 
                                 Vm_decay_ehigh = 0.4,
                                 b1Rd = -1.05, #median different than the ed default
                                 b2Rd = 0.325, #median different than the ed default
                                 q = 1,
                                 root_beta = 0.05, #median different than the ed default which is 0.0010000000
                                 SRA = 48,
                                 stoma_psi_c = 3,
                                 #root_resp = 0.14, #correct parameter name?
                                 orient_factor = 0.0,
                                 clumping_factor = 0.7,
                                 #from here, the normal distributions
                                 wood_Kexp = 2, #median different than the ed default
                                 stoma_psi_b = -160, #median different than the ed default
                                 b1Bl = 0.049,
                                 b2Bl = 1.89,
                                 #lognormal distributions
                                 wood_Kmax = 0.05, # my estimation is 0.06059039, #default is 0.014 
                                 wood_water_cap = 0.0074, #median different than the ed default, which is  0.0057100886
                                 leaf_water_cap = 0.00075, #median different than the ed default, which is 0.0021480524
                                 leaf_psi_tlp = -225.8, #median different than the ed default
                                 
                                 #radiation related:
                                 leaf_reflect_vis = 0.07190379,
                                 leaf_reflect_nir = 0.5042545,
                                 leaf_trans_vis = 0.006900955,
                                 leaf_trans_nir = 0.2408585,
                                 #weibull
                                 root_turnover_rate = 1.27,
                                 #Vm0 = 21.47, #different from standard, which is 9.0970001221
                                 #gamma distributions
                                 wood_psi50 = -138.2052, #different from ed median, which is -116.4605331421
                                 rho = 0.4638128
))
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
ed2in_scenar$EXPNME = 'SA water traits test'                                                  #newly added


write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN")) #writes the ED2IN file based on the inputs of ed2in_scenar

# Config file
config_simu <- config

#for loop that goes over all the parameters. If the parameter is Vm0 (the maximum capacity of Rubisco to perform the carboxylase function),
#a separate loop is executed. 

#We only loop over the liana pft, so the others can remain constant

#deze for loop dan weg om enkel voor lianen aan te passen?

for (iparam in seq(1,Nparam)){
  param_name <- param_names[iparam]
  
  if(param_name == "Vm0"){
    samples <- prior[[iparam]]$sample(n = 100000)
    params <- median(samples) #for default parameters, all pfts get the same value
    params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name]) #checks if a value lies above or below the global threshold
    
  } else if (param_name %in% c("leaf_reflect_vis", "leaf_reflect_nir", "leaf_trans_vis", "leaf_trans_nir")){
    samples <- prior[[iparam]]$sample(n = 100000)
    params <- median(samples) #for default parameters, all pfts get the same value
    params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name]) #checks if a value lies above or below the global threshold
    
  } else {
    samples <- prior[[iparam]]$sample(n = 100000)
    params_actual <- median(samples) #for default parameters, all pfts get the same value
    
  }
  
  config_simu[[4]][param_name] <- params_actual
  
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



#hier ook dan? enkel lianen aanpassen? enkel de default loop?

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
    ed2in_scenar$EXPNME = 'SA water traits test'                                                   #newly added
    
    write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN"))
    
    
    # Config
    config_simu <- config_reference
    
    for (iparam2 in seq(1,Nparam)){
      
      param_name = param_names[iparam2]
      
      if (param_name == param_names[iparam]){ #only execute if the parameter is the same as above: we only want to change one parameter at a time.
        if (param_names[iparam] == "Vm0"){ # With global maximum and minimum
          
          params <- quantile(prior[[iparam]]$sample(n = 100000),quantiles[iquantile])
          params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name]) #checks if a value lies above or below the global treshold
          
        } else if (param_name %in% c("leaf_reflect_vis", "leaf_reflect_nir", "leaf_trans_vis", "leaf_trans_nir")){
          
          params <- quantile(prior[[iparam]]$sample(n = 100000),quantiles[iquantile])
          params_actual <- pmax(pmin(params,global_max[param_name]),global_min[param_name])
          
        } else {
          
          samples <- prior[[iparam]]$sample(n = 100000)
          params_actual <- quantile(samples,quantiles[iquantile]) #here we get the value of the parameter for each quantile value that is specified.
        }
        
        config_simu[[4]][param_name] <- params_actual
        
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


