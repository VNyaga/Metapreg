library(haven)
library(ggplot2)

library(dplyr)
library(tidyr)

library(forcats)
library(simhelpers)
library(gghalves)


review <- "nissen2010"
lo <- 0
up <- 2
scale <- 0.5
scale <- 0.25
scale2 <- 1

#=====================================

rootdir <- "C:/DATA/WIV/Projects/GitHub/Metapreg/RSM2025/Data"

graphpath <- paste(rootdir,'/', review, '/Graphs', sep='')

sim <- read_dta(paste(rootdir,'/', review, '/simulatedresults.dta', sep=''))

sim$inf <- sim$`Inf`

sim$CImethod <- ifelse(sim$CImethod=="na" & sim$package=="metafor", "z", sim$CImethod )

#===================================All=====================================

sim$id <- with(sim, paste(package, Env, Dist, param, Covariance, Weighting, Sigmethod, Prior, inf, sep='-'))

orange = "CS"  

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
                    ((sim$stat == 'median-orout' & sim$inf=="B") | 
                       (sim$stat == 'mean-orout' & sim$inf=="F"))  & 
                    sim$link == "logit" & 
                    sim$Covariance == orange  & sim$Dist == "BN" &
                    sim$Design == "comparative")
               ,]
  
  source(paste(rootdir, '/Scripts/perfomance-stats.R', sep=""))

  min <- lo
  max <- up
  
  rbias$left1 <- min - 0.25 
  rbias$left2 <- min -1
  rbias$left3 <- min -1.75
  rbias$left4 <- min -2.0
  rbias$left5 <- min -2.75
  rbias$left6 <- min -3.25
  rbias$left7 <- min -4.0
  rbias$left8 <- min -4.75
  
  rbias$right2 <- max - 0.75
  rbias$right1 <- max + 0.25 
  
  coverage$right3 <- max + 1.5
  coverage$right4 <- max + 2.75
  
  coverage$right5 <- max + 3.75
  rbias$right5 <- max + 5
  rbias$right6 <- max + 6
  
  #rbias$right5 <- max + 3
  
  mintick <- 0.5
  maxtick <- 2
  
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
#=================================conditional statistics=============================

sim$conditional <- with(sim,ifelse((package == "metapreg" & link == "logit" & 
                                      ((Covariance == "CS"  & Design == "comparative" & 
                                          Dist == "BN" &
                                          (Prior == "ig(0.01, 0.01)" | Prior == "")) | 
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

sim$conditional <- with(sim,ifelse(package == "metafor" & Dist == "BN" & (Covariance == "CS" ), 1, conditional))


simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(package, Env, Dist, param, Sigmethod, Prior, inf, sep='-'))

#==========================Performance===========================================
source(paste(rootdir, '/Scripts/perfomance-stats.R', sep=""))
#===================================================================================

min <- 0.5
max <- 2

rbias$left1 <- min - 0.25 
rbias$left2 <- min -1
rbias$left3 <- min -1.25
rbias$left4 <- min -1.75
rbias$left5 <- min -2


rbias$right1 <- max + 0.10
rbias$right2 <- max + 0.5 

coverage$right3 <- max + 1.25
coverage$right4 <- max + 2

mintick <- 0.5
maxtick <- 2


#==========================Plot Performance===========================================
source(paste(rootdir, '/Scripts/plot-perfomance-stats.R', sep=""))

windows ()
g

#============================LINK + COVARIANCE + Misspecification=============================

sim$conditional <- with(sim,ifelse((package == "metapreg" & 
                                      Design %in% c("comparative", "general") 
                                    & Dist %in% c("B1", "BN")), 
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

rbias$right1 <- trueor + 1.5 

coverage$right3 <- trueor + 2.75
coverage$right4 <- trueor + 4.00

formatted_trueor <- format((min(sim$trueor)), digits=3)

minX <- min
maxX <- min(coverage$right4)

rbias$rel_bias <- format(rbias$rel_bias, digits=3)
rbias$rel_mse <- format(rbias$rel_mse, digits=1)
#==========================Plot Performance===========================================
orange = "CS"

source(paste(rootdir, '/Scripts/Link-Cov-Mispecification.R', sep=""))

windows()
F