#rm(list = ls())
#Second script to be run for the SA

library(dplyr)
library(LidarED)   # https://github.com/femeunier/LidarED/
library(purrr)
library(stringr)

outdir <- "/data/gent/vo/000/gvo00074/ED2kasper/out/SA"

df.seasonal.C.all <- df.param.all <- data.frame()

#under: have to be the same as in the first script.
params <- c("stomatal_slope","D0","Vm0","vm_q10","Vm_high_temp","Vm_decay_ehigh")
Nparam <- length(params)

quantiles <- c(0.025,0.16,0.25,0.75,0.84,0.975)

# SA runs
for (iparam in seq(1,Nparam)){
  for (iquantile in seq(1,length(quantiles))){

    param_name <- params[iparam]
    run_name <- paste0("SA_",param_name,"_",quantiles[iquantile]) #the naming from the first script.
    out_ref <- file.path(outdir,run_name)

    datum.file <- file.path(out_ref,"analy","analysis.RData")

    if (file.exists(datum.file)){
      load(datum.file)
      #get the information from the simulation: extract what you want.

      ED2.seasonal.C <- data.frame(year = datum$year,
                                   month = datum$month,
                                   GPP = datum$emean$gpp,
                                   Reco = datum$emean$reco,
                                   Rauto = datum$emean$plant.resp,
                                   Rhetero = datum$emean$het.resp,
                                   NEP = -datum$emean$nee/1000*1e-6*12*86400*365,
                                   NPP = datum$emean$npp)

      history.file <-  file.path(out_ref,"histo","history.xml")

      if (params[iparam] == "Vm0") {
        df.param <- data.frame(pft = c(2, 3, 4),
                               param = params[iparam],
                               value = unlist(map(c(2, 3, 4), function(x)
                                 get_ED_default_pft(history.file, params[iparam], pft_target = x)))) #get the value for each pft separately, as they differ.
      } else{
        df.param <-
          data.frame(pft = c(2, 3, 4),
                     param = params[iparam],
                     value = rep(get_ED_default_pft(history.file, params[iparam], pft_target = 2), 3)) #here it is always the same value,
        #so we repeat the value of the second pft three times.
      }

      df.param.all <- bind_rows(list(df.param.all,
                                     df.param %>% mutate(quantile = quantiles[iquantile],
                                                         param = params[iparam])))



      df.seasonal.C.all <- bind_rows(list(df.seasonal.C.all,
                                          ED2.seasonal.C %>% mutate(quantile = quantiles[iquantile],
                                                                    param = params[iparam])))
      #here: get the values from the datum file that correspond to the parameter and the quantile.
    }
  }
}

# Reference run
#same as above but for the reference run.
run_name <- "SA_reference"
out_ref <- file.path(outdir,run_name)

datum.file <- file.path(out_ref,"analy","analysis.RData")
history.file <-  file.path(out_ref,"histo","history.xml")

if (file.exists(datum.file)){
  load(datum.file)

  ED2.seasonal.C <- data.frame(year = datum$year,
                               month = datum$month,
                               GPP = datum$emean$gpp,
                               Reco = datum$emean$reco,
                               Rauto = datum$emean$plant.resp,
                               Rhetero = datum$emean$het.resp,
                               NEP = -datum$emean$nee/1000*1e-6*12*86400*365,
                               NPP = datum$emean$npp)

  for (iparam in seq(1,Nparam)){
    if (params[iparam] == "Vm0") {
      df.param <- data.frame(pft = c(2, 3, 4),
                             param = params[iparam],
                             value = unlist(map(c(2, 3, 4), function(x)
                               get_ED_default_pft(history.file, params[iparam], pft_target = x))))
    } else{
      df.param <-
        data.frame(pft = c(2, 3, 4),
                   param = params[iparam],
                   value = rep(get_ED_default_pft(history.file, params[iparam], pft_target = 2), 3))
    }
    df.param.all <- bind_rows(list(df.param.all,
                                   df.param %>% mutate(quantile = 0.5,
                                                       param = params[iparam])))
  }



  df.seasonal.C.all <- bind_rows(list(df.seasonal.C.all,
                                      ED2.seasonal.C %>% mutate(quantile = 0.5,
                                                                param = "reference"))) #not in the for loop here because the quantile and the parameter
  #are always the same.

}

saveRDS(object = df.seasonal.C.all,file = file.path('.',"df.seasonal.C.SA.RDS"))
saveRDS(object = df.param.all,file = file.path('.',"df_param_SA_YGB.RDS"))

# scp /home/femeunier/Documents/projects/YGB/scripts/analyze_SA.R hpc:/data/gent/vo/000/gvo00074/felicien/R

