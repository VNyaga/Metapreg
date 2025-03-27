library(haven)
library(ggplot2)

library(dplyr)
library(tidyr)

library(forcats)
library(simhelpers)
library(gghalves)

#=================diuretics ====================
review <- "diuretics"
lo <- 0
up <- 1
scale <- 0.5
scale <- 0.25
scale2 <- 1

#=====================================

rootdir <- "C:/DATA/WIV/Projects/GitHub/Metapreg/Data"

graphpath <- paste(rootdir,'/', review, '/Graphs', sep='')

sim <- read_dta(paste(rootdir,'/', review, '/simulatedresults.dta', sep=''))

sim$inf <- sim$`Inf`


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
  
  min <- lo
  max <- up
  
  
  rbias$left1 <- (trueor - 0.5 )
  rbias$left2 <- (trueor - 1.0)
  rbias$left3 <- (trueor -1.25)
  rbias$left4 <- (trueor -1.5)
  rbias$left5 <- (trueor -1.75)
  rbias$left6 <- (trueor -2.0)
  rbias$left7 <- (trueor -2.25)
  rbias$left8 <- ( trueor-2.5)
  
  rbias$right1 <- (trueor + 0.85 )
  rbias$right2 <- (trueor + 1.25 )
  coverage$right3 <- (trueor + 1.5 )
  coverage$right4 <- (trueor + 2.15)
  coverage$right5 <- (trueor + 2.75)
  
  rbias$right5 <- (trueor + 2.75)
  rbias$right6 <- (trueor + 3.25)
  
  
  #rbias$right5 <- max + 3
  
  mintick <- 0
  maxtick <- 1
  
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

#=================================Selected Estimators=============================
sim$conditional <- with(sim,ifelse((package == "metapreg" & link == "logit" & 
                                      ((Covariance == "IND"  & Design == "comparative" & 
                                          Dist == "BN" & (Prior == "ig(0.01, 0.01)" | Prior == "")) | 
                                         Dist %in% c('BB1', 'BB2'))), 
                                   1,0))



sim$conditional <- with(sim,ifelse(inference == "bayesian" & package != "metapreg" & Dist != "B", 1, conditional))
#sim$conditional <- with(sim,ifelse( Weighting == "SSW" & package == "meta" & Env == "R", 1, conditional))
sim$conditional <- with(sim,ifelse((Sigmethod == "sj" | Sigmethod == "mp") & 
                                     Weighting == "SSW" & package == "meta" & Env == "R", 1, conditional))
sim$conditional <- with(sim,ifelse(slope == "common", 0, conditional))
sim$conditional <- with(sim,ifelse((Dist == "QN" | Sigmethod == "mp" | Sigmethod == "pl") & package == "metan"  , 1, conditional))
sim$conditional <- with(sim,ifelse(Sigmethod == "sj" & CImethod == "tkh" & package == "meta" & Env == "Stata", 1, conditional))
sim$conditional <- with(sim,ifelse(package %in% c("metaplus") == TRUE, 1, conditional))
sim$conditional <- with(sim,ifelse(package == "metafor" & Dist == "BN" & (Covariance == "IND" ), 1, conditional))


simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])
#5900
simor$id <- with(simor, paste(package, Env, Dist, param, Sigmethod, Prior, inf, sep='-'))
#==========================Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')
#===================================================================================
rbias$left1 <- (trueor - 0.5 )
rbias$left2 <- (trueor - 1.0)
rbias$left3 <- (trueor -1.25)
rbias$left4 <- (trueor -1.5)
rbias$left5 <- (trueor -1.75)
rbias$left6 <- (trueor -2.0)
rbias$left7 <- (trueor -2.25)
rbias$left8 <- ( trueor-2.5)

rbias$right1 <- (trueor + 0.75 )
rbias$right2 <- (trueor + 1.25 )
coverage$right3 <- (trueor + 1.75 )
coverage$right4 <- (trueor + 2.25)

coverage$right5 <- (trueor + 2.75)
coverage$right6 <- (trueor + 3.25)

mintick <- 0
maxtick <- 1

#==========================Plot Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/plot-perfomance-stats.R')

windows()
g

#============================LINK + COVARIANCE + Misspecification=============================

sim$conditional <- with(sim,ifelse((package == "metapreg" & 
                                      Design %in% c("comparative", "general") 
                                    & Dist %in% c("B1", "BN") & 
                                      (Prior == "ig(0.01, 0.01)" | Prior == "" | Prior == "iw(2,3,i(2))")), 
                                   1,0))

simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(Dist, Covariance, sep='-'))

#==========================Performance===========================================
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')

#=============================================
max <- ceiling(max(simor$esthat))
min <- floor(min(simor$esthat))

rbias$left1 <- (trueor - 0.5 )
rbias$left2 <- (trueor - 1.0)
rbias$left3 <- (trueor -1.25)
rbias$left4 <- (trueor -1.5)
rbias$left5 <- (trueor -1.75)
rbias$left6 <- (trueor -2.0)
rbias$left7 <- (trueor -2.25)
rbias$left8 <- ( trueor-2.5)

rbias$right2 <- (trueor + 0.75 )
rbias$right1 <- (trueor + 1.25 )
coverage$right3 <- (trueor + 1.75 )
coverage$right4 <- (trueor + 2.25)

coverage$right5 <- (trueor + 2.75)
coverage$right6 <- (trueor + 3.25)

trueor <- (min(sim$trueor))

rbias$right1 <- trueor + 0.5 

coverage$right3 <- trueor + 1.25
coverage$right4 <- trueor + 2.25

formatted_trueor <- format((min(sim$trueor)), digits=3)

minX <- min
maxX <- min(coverage$right4)

rbias$rel_bias <- format(rbias$rel_bias, digits=3)
rbias$rel_mse <- format(rbias$rel_mse, digits=1)
#==========================Plot Performance===========================================
orange = "IND"
source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/Link-Cov-Mispecification.R')

windows()
F
