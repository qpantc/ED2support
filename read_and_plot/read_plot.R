setwd("./R")

source("./post.process.ED2.outputs.R")

load("../outputs/Results/BCI/analy/analysis.RData")

matplot(datum$szpft$gpp[,12,c(2,3,4,18)],type = "l")
