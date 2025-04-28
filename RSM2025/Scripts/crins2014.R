library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(simhelpers)
library(gghalves)


review <- "crins201"
lo <- 0
up <- 1
scale <- 0.5
scale <- 0.25
scale2 <- 1

#=====================================

rootdir <- "C:/DATA/WIV/Projects/GitHub/Metapreg/RSM2025/"

graphpath <- paste(rootdir,'/', review, '/Graphs', sep='')

sim <- read_dta(paste(rootdir,'/', review, '/simulatedresults.dta', sep=''))

sim$inf <- sim$`Inf`


sim$CImethod <- with(sim, ifelse(package=='metapreg' & CImethod == "" & inf=="B", 'hpd', CImethod))
sim$CImethod <- with(sim, ifelse(package=='metapreg' & CImethod == "" & inf=="F", 'z', CImethod))

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
  
  source(paste(rootdir, '/Scripts/perfomance-stats.R', sep=""))

  min <- lo
  max <- up
  
  
  rbias$left1 <- (trueor - 0.35 )
  rbias$left2 <- (trueor - 0.85)
  rbias$left3 <- (trueor -1.25)
  rbias$left4 <- (trueor -1.5)
  rbias$left5 <- (trueor -2.0)
  rbias$left6 <- (trueor -2.5)
  rbias$left7 <- (trueor -3.0)
  rbias$left8 <- ( trueor-3.5)
  
  rbias$right1 <- (trueor + 0.5 )
  rbias$right2 <- (trueor + 0.75 )
  coverage$right3 <- (trueor + 1.5 )
  coverage$right4 <- (trueor + 2.25)
  coverage$right5 <- (trueor + 3)
  
  rbias$right5 <- (trueor + 3.25)
  rbias$right6 <- (trueor + 4)
  
  mintick <- 0
  maxtick <- 1
  
  if (run == 1) {
    #Stata suite
    source(paste(rootdir, '/Scripts/plot-smeta-perfomance-stats.R', sep=""))
  } else if (run == 2) { 
    #metan
    source(paste(rootdir, '/Scripts/plot-metan-perfomance-stats.R', sep=""))
  } else if (run == 3) {
    #metafor
    source(paste(rootdir, '/Scripts/plot-metafor-perfomance-stats.R', sep=""))
  } else if (run == 4) {
    #metastan
    source(paste(rootdir, '/Scripts/plot-metastan-perfomance-stats.R', sep=""))
  } else if (run == 5) {
    #meta
    source(paste(rootdir, '/Scripts/plot-R-meta-perfomance-stats.R', sep=""))
  } else if (run == 6) {
    #bayesmeta
    source(paste(rootdir, '/Scripts/plot-bayesmeta-perfomance-stats.R', sep=""))
  } else if (run == 7) {
    #metaplus
    source(paste(rootdir, '/Scripts/plot-metaplus-perfomance-stats.R', sep=""))
  }
}

#=================================Selected Estimators =============================
sim$conditional <- with(sim,ifelse((package == "metapreg" & link == "logit" & 
                                      ((Covariance == orange  & Design == "comparative" & 
                                          Dist == "BN" &
                                          (Prior == "ig(0.01, 0.01)" | Prior == "")) | 
                                         Dist %in% c('BB1', 'BB2', 'B2'))), 
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
source(paste(rootdir, '/Scripts/perfomance-stats.R', sep=""))
#===================================================================================
# coverage$left1 <- (min - 0.25 - 0.3)*0.25
# coverage$left2 <- (min -0.75 - 0.3)*0.25
# coverage$left3 <- (min -1.0- 0.3)*.25
# coverage$left4 <- (min -1.5- 0.3)*.25
# coverage$left5 <- (min -2.5- 0.3)*.25
# coverage$left6 <- (min -3.0- 0.3)*.25
# coverage$left7 <- (min -4.0- 0.3)*.25
# coverage$left8 <- (min -5.0- 0.3)*.25
# 
# rbias$right2 <- (max - 0.5)
# rbias$right1 <- (max - 0.3)
# coverage$right3 <- (max - 0.1 )
# coverage$right4 <- (max + 0.1)
# 
# coverage$right5 <- (max + 0.3)
# coverage$right6 <- (max + 0.5)


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
maxtick <- 0.5
#==========================Plot Performance===========================================
source(paste(rootdir, '/Scripts/plot-perfomance-stats.R', sep=""))

windows()
g
#============================LINK + COVARIANCE + Misspecification=============================

sim$conditional <- with(sim,ifelse((package == "metapreg" & 
                                      Design %in% c("comparative", "general") & 
                                      Dist %in% c("B1", "BN") & 
                                      (Prior == "ig(0.01, 0.01)" | Prior == "" | Prior == "iw(2,3,i(2))")), 
                                   1,0))

simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(Dist, Covariance, sep='-'))

#==========================Performance===========================================
source(paste(rootdir, '/Scripts/perfomance-stats.R', sep=""))

#=============================================
max <- ceiling(max(simor$esthat))
min <- floor(min(simor$esthat))


trueor <- (min(sim$trueor))

rbias$right1 <- trueor + 0.5 

coverage$right3 <- trueor + 1.25
coverage$right4 <- trueor + 2.25


minX <- min
maxX <- min(coverage$right4)

rbias$rel_bias <- format(rbias$rel_bias, digits=3)
rbias$rel_mse <- format(rbias$rel_mse, digits=1)
#==========================Plot Performance===========================================
orange = "IND"
source(paste(rootdir, '/Scripts/Link-Cov-Mispecification.R', sep=""))

windows()
F

