library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(simhelpers)
library(gghalves)

review <- "bender2018b"

rootdir <- "C:/DATA/WIV/Projects/GitHub/Metapreg/Data"

graphpath <- paste(rootdir,'/', review, '/Graphs', sep='')

#===========================================SIMULATION=====================================
sim <- read_dta(paste(rootdir,'/', review, '/simulatedresults.dta', sep=''))

lo <- 0
up <- 10

sim$inf <- sim$`Inf`
sim <- sim[!duplicated(cbind(sim$sim, sim$model, sim$stat)), ]

sim$CImethod <- ifelse(sim$CImethod=="na" & sim$package=="metafor", "z", sim$CImethod )

#===================================All=====================================

sim$id <- with(sim, paste(package, Env, Dist, param, Covariance, Weighting, Sigmethod, Prior, inf, sep='-'))

orange = "IND"  

for (run in 1:7) {
  #package
  if (run == 1) {
    package <- "meta"
    
  } else if (run == 2) {
    package <- "metan"
    
  } else if (run == 3) {
    package <- "metafor"
    
  } else if (run == 4) {
    package <- "metastan"
    
  } else if (run == 5) {
    package <- "meta"
    
  } else if (run == 6) {
    package <- "bayesmeta"
    
  } else if (run == 7) {
    package <- "metaplus"
  }
  
  #Environment
  if (run < 3) {
    Env <- "Stata" 
  } else {
    Env <- "R"  
  }
  
  simor <- sim[(sim$package==package & 
                  sim$Env==Env & ((sim$stat == 'median-orout' & sim$inf=="B") | 
                                    (sim$stat == 'mean-orout' & sim$inf=="F"))) |
                 
                 (sim$package=="metapreg" & 
                    ((sim$stat == 'median-orout' & sim$inf=="B" & sim$Prior == "ig(0.01, 0.01)") | 
                       (sim$stat == 'mean-orout' & sim$inf=="F"))  & 
                    sim$link == "logit" & 
                    sim$Covariance == orange  & sim$Dist == "BN" &
                    sim$Design == "comparative")
               ,]
  
  source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')
  
  min <- 0
  max <- 10
  rbias$left1 <- min - 0.3 - 0.3
  rbias$left2 <- min -1.6- 0.3
  rbias$left3 <- min -2.0- 0.3
  rbias$left4 <- min -3.8- 0.3
  rbias$left5 <- min -4.5- 0.3
  rbias$left6 <- min -5.8- 0.3
  rbias$left7 <- min -7.3- 0.3
  rbias$left8 <- min -8.9- 0.3
  
  rbias$right2 <- max - 2.25
  rbias$right1 <- max -2
  
  coverage$right3 <- max + 0.5
  coverage$right4 <- max + 3
  
  coverage$right5 <- max + 5.5
  rbias$right5 <- max + 7.15
  rbias$right6 <- max + 8
  
  mintick <- 2
  maxtick <- 10
  
  if (run == 1) {
    #Stata suite
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-smeta-perfomance-stats.R")
  } else if (run == 2) { 
    #metan
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-metan-perfomance-stats.R")
  } else if (run == 3) {
    #metafor
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-metafor-perfomance-stats.R")
  } else if (run == 4) {
    #metastan
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-metastan-perfomance-stats.R")
  } else if (run == 5) {
    #meta
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-R-meta-perfomance-stats.R")
  } else if (run == 6) {
    #bayesmeta
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-bayesmeta-perfomance-stats.R")
  } else if (run == 7) {
    #metaplus
    source("C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-metaplus-perfomance-stats.R")
  }
}
#=================================Selected Estimators =============================
sim$conditional <- with(sim,ifelse((package == "metapreg" & link == "logit" & 
                             ((Covariance == "IND"  & Design == "comparative" & 
                             Dist == "BN") | Dist %in% c('BB1', 'BB2', 'B2'))), 
                          1,0))

sim$conditional <- with(sim,ifelse(inference == "bayesian" & 
                                     package != "metapreg" & 
                                     Dist != "B", 1, conditional))

sim$conditional <- with(sim,ifelse((Sigmethod == "sj"| Sigmethod == "mp") & 
                                     Weighting == "SSW" & package == "meta" & 
                                     Env == "R", 1, conditional))

sim$conditional <- with(sim,ifelse(slope == "common", 0, conditional))

sim$conditional <- with(sim,ifelse((Dist == "QN" | Sigmethod == "pl") & 
                                     package == "metan"  , 1, conditional))

sim$conditional <- with(sim,ifelse((Sigmethod == "sj"| Sigmethod == "mp") & CImethod %in% c("t", "tkh") & 
                                     package == "meta" & Env == "Stata", 1, conditional))

sim$conditional <- with(sim,ifelse(package %in% c("metaplus") == TRUE, 1, conditional))

sim$conditional <- with(sim,ifelse(package == "metafor" & Dist == "BN" & 
                                     (Covariance == "IND" ), 1, conditional))


simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(package, Env, Dist, param, Sigmethod, Prior, inf, sep='-'))

#==========================Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')

#==================================================================

min <- 1.5
max <- 10
rbias$left1 <- min - 1
rbias$left2 <- min - 2.5
rbias$left3 <- min -3.5
rbias$left4 <- min -4.5
rbias$left5 <- min -5.5

rbias$right1 <- max + 1
rbias$right2 <- max + 2.15

coverage$right3 <- max + 3.75
coverage$right4 <- max + 5.5

mintick <- 2
maxtick <- 10

#==========================Plot Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-perfomance-stats.R')

windows()
g
#============================LINK + COVARIANCE + Misspecification=============================

sim$conditional <- with(sim,ifelse((package == "metapreg" & 
                                      Design %in% c("general", "comparative") &
                                      Dist %in% c("B1", "BN")), 
                                   1,0))

simor <- data.frame(sim[(sim$stat %in% c('median-orout', 'mean-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(Dist, Covariance, sep='-'))


# 
# sim$conditional <- with(sim,ifelse((package == "metapreg" & 
#                                       Dist %in% c("BB1", "BB2")), 
#                                    1,0))
# 
# simor <- data.frame(sim[(sim$stat %in% c('median-orout', 'mean-orout')==TRUE) & 
#                           (sim$conditional==1),])
# 
# simor$id <- with(simor, paste(Dist, Covariance, sep='-'))

#==========================Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')

#=============================================
max <- ceiling(max(simor$esthat))
min <- floor(min(simor$esthat))


trueor <- (min(sim$trueor))

rbias$right1 <- trueor + 2.5 

coverage$right3 <- trueor + 3.75
coverage$right4 <- trueor + 5.00

formatted_trueor <- format((min(sim$trueor)), digits=3)

minX <- min
maxX <- min(coverage$right4)

rbias$rel_bias <- format(rbias$rel_bias, digits=3)
rbias$rel_mse <- format(rbias$rel_mse, digits=1)
#==========================Plot Performance===========================================
orange = "IND"
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/Link-Cov-Mispecification.R')

windows()
B