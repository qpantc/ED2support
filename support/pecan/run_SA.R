#rm(list = ls())

setwd("~/ED2/Sensitivity analysis/Univariate SA")
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(BayesianTools)

source("sensitvity.analysis.R") #to get the functions that are listed in this file.

# Import outputs and parameter data.frame
#This
system2("rsync",paste("-avz","hpc:/data/gent/vo/000/gvo00074/ED2kasper/R/Sensitivity_analysis/df.seasonal.C.SA.RDS",
                      "./outputs/"))
system2("rsync",paste("-avz","hpc:/data/gent/vo/000/gvo00074/ED2kasper/R/Sensitivity_analysis/df_param_SA_YGB.RDS",
                      "./outputs/"))
#system2: invokes the OS command specified by command, rsync: rsync is a utility for efficiently transferring and synchronizing files between 
#a computer and an external hard drive and across networked computers by comparing the modification times and sizes of files
#or this
# system2("rsync",paste("/data/gent/vo/000/gvo00074/ED2kasper/R/df.seasonal.C.SA.RDS",
#                       "./outputs/"))
# system2("rsync",paste("/data/gent/vo/000/gvo00074/ED2kasper/R/df_param_SA_YGB.RDS",
#                       "./outputs/"))


# Priors
pft_lowers <- c(stomatal_slope = 2,
                D0 = 0.005,
                Vm0 = 5,
                vm_q10 = 1.8,
                Vm_high_temp = 40,
                Vm_decay_ehigh = 0.2)

pft_uppers <- c(stomatal_slope = 16,
                D0 = 0.03,
                Vm0 = 35,
                vm_q10 = 3,
                Vm_high_temp = 50,
                Vm_decay_ehigh = 0.6)

#same priors as in the first script, but does not seem to work...
# pft_lowers <- c(stomatal_slope = 2,
#                 D0 = 0.005,
#                 Vm0 = 5,
#                 Delta_Vm0 = 0,
#                 vm_q10 = 1.8,
#                 Vm_high_temp = 40,
#                 Vm_decay_ehigh = 0.2)
# 
# pft_uppers <- c(stomatal_slope = 16,
#                 D0 = 0.03,
#                 Vm0 = 35,
#                 Delta_Vm0 = 10,
#                 vm_q10 = 3,
#                 Vm_high_temp = 50,
#                 Vm_decay_ehigh = 0.6)


prior <- map2(pft_lowers,pft_uppers,createUniformPrior) # for each parameter, a lower, upper and best value is created
param_names <- names(pft_lowers)
Nparam <- length(param_names)

# Trait samples
trait.samples <- list()
for (iparam in seq(1,Nparam)){
  trait.samples[[param_names[iparam]]] <- as.vector(prior[[iparam]]$sample(n = 100000)) #100000 samples drawn from each distribution (uniform in this case).
} #for all parameters, gets 100000 samples from the distribution, needed for the sensitivity_analysis function

# sa.samples

df.samples <- readRDS("./outputs/df_param_SA_YGB.RDS") #for each PFT, parameter and quantile, a value is provided (for the model outputs that were defined in de previous script).
#for stomatal slope, the quantiles above 0.5 are missing, for D0, the first quantile is missing. NOT ANYMORE (but didn't change anything)

df.samples.current <- df.samples
df.samples.current.all <- df.samples.current %>% pivot_wider(names_from = param,
                                                             values_from = value) %>% arrange(quantile)
sa.samples <- df.samples.current.all #needed for the sensitivity analysis function

###############################
#change sa.samples#############
###############################
sa.samples2 <- sa.samples %>% group_by(quantile) %>% summarise(stomatal_slope = mean(stomatal_slope),
                                                               D0 = mean(D0),
                                                               Vm0 = mean(Vm0),
                                                               vm_q10 = mean(vm_q10),
                                                               Vm_high_temp = mean(Vm_high_temp),
                                                               Vm_decay_ehigh = mean(Vm_decay_ehigh))
##############################

# SA outputs original
df.SA <- readRDS("./outputs/df.seasonal.C.SA.RDS") #columns: year, month, GPP, Reco, Rauto, Rhetero, NEP, NPP, quantile, param
# a file with quantile values and then the 6 parameters, 48 values for each parameter per quantile, this is because there is one value per simulation month.and we have 36 months

df.SA.sum <- df.SA %>% group_by(param,quantile) %>% summarise(m = mean(GPP)) #reference is one value, and then a different value for all the quantiles
#of all the parameters, the GPP is the mean per quantile and parameter of the whole simulation, so 36 values.
df.SA.current <- df.SA.sum
df.SA.current.all <- bind_rows(list(do.call("rbind",df.SA.current %>% filter(param == "reference") %>% replicate(n = Nparam,simplify = FALSE)) %>%
                                      mutate(param = unique(df.SA.current %>% filter(param != "reference") %>% pull(param))),
                                    df.SA.current %>% filter(param != "reference"))) %>%
  pivot_wider(names_from = param, values_from = m) %>% arrange(quantile)
#first line repeats the value from the reference for 6 times (Nparam), and the second line gives them the name of the different parameters that
#are used for the SA. The third line adds this to the existing results. The pivot_wider is applied, everything is arranged by quantile.

sa.output <- df.SA.current.all #needed for the sensitivity analysis function.

#########################################################
#SA outputs for all pfts combined: take the sum of gpp###
#########################################################
df.SA2 <- readRDS("./outputs/df.seasonal.C.SA.RDS") 
df.SA2.sum <- df.SA %>% group_by(param,quantile) %>% summarise(m = sum(GPP)) 
df.SA2.current <- df.SA2.sum
df.SA2.current.all <- bind_rows(list(do.call("rbind",df.SA2.current %>% filter(param == "reference") %>% replicate(n = Nparam,simplify = FALSE)) %>%
                                      mutate(param = unique(df.SA2.current %>% filter(param != "reference") %>% pull(param))),
                                    df.SA2.current %>% filter(param != "reference"))) %>%
  pivot_wider(names_from = param, values_from = m) %>% arrange(quantile)

sa.output2 <- df.SA2.current.all 
#########################################################

outdir <- file.path("./outputs/") #needed for the sensitivity analysis function.


#do the actual sensitivity analysis:
SA.op <- sensitivity.analysis(trait.samples, sa.samples = sa.samples %>% filter(pft == 2), sa.output, outdir)
print(SA.op$variance.decomposition.output$partial.variances)

#all_gpp
# #change sa.output --> wrong!?
# SA.op <- sensitivity.analysis(trait.samples, sa.samples , rbind(sa.output2,sa.output2,sa.output2) %>% arrange(quantile), outdir)
# print(SA.op$variance.decomposition.output$partial.variances)
#change sa. samples
SA.op <- sensitivity.analysis(trait.samples, sa.samples2 , sa.output2 , outdir)
print(SA.op$variance.decomposition.output$partial.variances)


##' This function estimates the univariate responses of a model to a parameter for a set of traits, calculates the model sensitivity at the median,
##'  and performs a variance decomposition. This function results in a set of sensitivity plots (one per variable) and plot_variance_decomposition.
##'  
##' Spline estimate of univariate relationship between parameter value and model output
##' Creates a spline function using the splinefun function that estimates univariate response of parameter input to model output
##' 
##' 

#kurtosis: kurtosis describes the shape of a probability distribution and there are different ways of quantifying it for a theoretical distribution
           #and corresponding ways of estimating it from a sample from a population. Different measures of kurtosis may have different interpretations
#In sensitivity.analysis function:
    # In mathematics, a spline is a special function defined piecewise by polynomials. In interpolating problems, spline interpolation is often
    # preferred to polynomial interpolation because it yields similar results, even when using low degree polynomials, while avoiding Runge's 
    # phenomenon for higher degrees.
#spline.truncate: removes values lower than zero if (quantile(x, min.quantile) > 0)

#coef.vars = sqrt(var(set)) / median(set)
#the coefficient of variation (CV), also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution. It is often expressed as a percentage, and is defined as the ratio of the standard deviation {\displaystyle \ \sigma }\ \sigma  to the mean {\displaystyle \ \mu }\ \mu  (or its absolute value, {\displaystyle |\mu |}|\mu |).



#extra code from from Félicien to plot
df.SA.decomp <- data.frame(params = as.character(names(SA.op$variance.decomposition.output$coef.vars)),
                           CV = SA.op$variance.decomposition.output$coef.vars,
                           elast = SA.op$variance.decomposition.output$elasticities,
                           par.var = SA.op$variance.decomposition.output$partial.variances) %>% arrange(par.var) %>%
  mutate(params = factor(params,levels = params)) %>%
  pivot_longer(cols = c("CV","elast","par.var"),
               names_to = "type",
               values_to = "value")



# Plot
ggplot(data = df.SA.decomp) +
  geom_point(aes(x = value,y = params)) +
  geom_segment(aes(x = 0, y = params,
                   xend = value,yend = params)) +
  facet_wrap(~ type,scales = "free_x") +
  theme_bw()



