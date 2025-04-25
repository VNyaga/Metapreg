
sim <- data$sim[1]
mu0true <- data$mu0true[1]
ortrue <- data$ortrue[1]
tausqtrue <- data$tausqtrue[1]
sigmasqtrue <- data$sigmasqtrue[1]
nstudies <- data$nstudies[1]
studysize <- data$studysize[1]

repos <- "https://www.freestatistics.org/cran/"

packages <- c("bayesmeta", "metafor", "lme4", "BiasedUrn", "meta", 
              "metasens", "rema","RandMeta", "MetaStan", "metaSEM",
              "metaplus", "metaLik", "metaBMA")


#try(install.packages(packages, repos = repos, dependencies = TRUE), silent = FALSE)

event0 <- data$event0
event1 <- data$event1
total0 <- data$total0
total1 <- data$total1
slab <- data$studyid

esthat <- "."
esthatlo <- "." 
esthatup <- "." 

tau2hat <- "."
tau2hatlo <- "." 
tau2hatup <- "." 

sigma2hat <- "."
sigma2hatlo <- "." 
sigma2hatup <- "." 
IC <- "."


#======================================Metafor===============
suppressWarnings(library(metafor))
package <- "metafor"
inference <- "frequentist"
parameter <- "mean-orout"
tests <- c("z", "t")
run <- 0
models <- c("CM.EL", "CM.AL", "UM.FS", "UM.RS")
methods <- c("ML", "EE") #for umrs
#typecc <- c("none", "all", "only0", "if0all")
typecc <- c("only0")
#cor=TRUE, if UM.RS + ML

for (model in 1:length(models)) {
  if (models[model] == "UM.RS") {
    nmethods <- 2
  } else {
    nmethods <- 1
  }
  
  for (method in 1:nmethods) {
    if (models[model] == "UM.RS" & methods[method] == "ML"){
      nAGQ <- 1
      loops <- 2
    } else {
      loops <- 1
      nAGQ <- 7
    }
    
    for (cc in 1:length(typecc)) {
      for (loop in 1:loops) {
        
        if (loop == 2) corr <- TRUE else corr <- FALSE
        
          for (test in 1:length(tests)){
        
            try({
            #Options for independent covariance
            # corr <- FALSE
            # model <- 4
            # nAGQ <- 1
            # cc <- 1
            # method <- 1
            # test <- 1
            # 
            # start <- Sys.time()
              
            res1 <- rma.glmm(measure="OR", 
                             ai = event1 , 
                             bi = ( total1 - event1 ) , 
                             ci= event0 , 
                             di =( total0 - event0 ), 
                             model= models[model], 
                             method = methods[method], 
                             cor = corr, 
                             nAGQ = nAGQ,
                             to = typecc[cc],
                             coding = 0,
                             test=tests[test])
            # # Poisson
            # start <- Sys.time()
            # res1 <- rma.glmm(measure="IRR",
            #                  x1i = event1 ,
            #                  t1i =  total1 ,
            #                  x2i= event0 ,
            #                  t2i = total0 ,
            #                  model= models[model],
            #                  method = methods[method],
            #                  cor = corr,
            #                  nAGQ = nAGQ,
            #                  to = typecc[cc],
            #                  coding = 0,
            #                  test=tests[test])
            # 
            # summary(res1)
            # forest(res1, slab=slab,refline=1, transf = exp, header=TRUE, addpred=TRUE)
            # print( Sys.time() - start )
            
            citype <- tests[test]
            k <- res1$k
            tau2hat <- res1$sigma2
            
            sigma2hat <- res1$tau2
            sigma2hatlo <- res1$ci.lb.tau2
            sigma2hatup <- res1$ci.ub.tau2
            
            esthat <- exp(res1$b) 
            esthatlo <- exp(res1$ci.lb) 
            esthatup <- exp(res1$ci.ub)
            
            modeli <- models[model]
            methodi <- methods[method]
            
            line <-  paste(sim, ";",
                           paste(package, "-", inference, "-",  modeli,  "-", methodi,  "-", citype,  "-", "-", loop, "-", NA,  "-", NA, sep=""), ";",
                           parameter, ";",
                					esthat, ";", esthatlo, ";", esthatup, ";",
                					tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
                					sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
                					IC, ";",
                					mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
                					k, ";", nstudies, ";", studysize, 
                					sep="")
                					
    
            write(line, file=outfile, append = TRUE)
            run = run + 1
            print(run)
          }, silent=TRUE)
        }
      }
    }
  }
}
#==============================================ESCALC===============
# start0 <- Sys.time()
effs <- escalc(measure="OR",
               ai=event1,  
               n1i=total1,
               ci=event0, 
               n2i=total0,
               slab=slab, 
               to="only0")

# print( Sys.time() - start0)
# 
# start0 <- Sys.time()
# effsRR <- escalc(measure="RR",
#                ai=event1,
#                n1i=total1,
#                ci=event0,
#                n2i=total0,
#                slab=slab,
#                to="only0")
# 
# print( Sys.time() - start0)
# 
# start0 <- Sys.time()
# effsRD <- escalc(measure="RD",
#                  ai=event1,
#                  n1i=total1,
#                  ci=event0,
#                  n2i=total0,
#                  slab=slab,
#                  to="only0")
# print( Sys.time() - start0)
# 
# start0 <- Sys.time()
# effsASD <- escalc(measure="AS",
#                  ai=event1,
#                  n1i=total1,
#                  ci=event0,
#                  n2i=total0,
#                  slab=slab,
#                  to="only0")
# print( Sys.time() - start0)

index <- !is.na(effs$yi) & !is.na(effs$vi)

#nocc
K <- sum(!is.na(effs$yi) & !is.na(effs$vi))
effs <- effs[index,]
#=================================================bayesmeta===============
suppressWarnings(library(bayesmeta))
package <- "bayesmeta"
inference <- "bayesian"
model <- "RE"
citype <- c("shortest", "central")
run <- 0


  for (loop in 1:2) {
    if (loop == 1) {
      taupriordensity <- function(t){dhalfcauchy(t,scale=1)}
      
    } else {
      taupriordensity <- function(t){dhalfnormal(t,scale=1)}
    }
    
    tau2hat <- "."
    tau2hatlo <- "." 
    tau2hatup <- "." 
    
    sigma2hat <- "."
    sigma2hatlo <- "." 
    sigma2hatup <- "." 
    
    
      for (ci in 1:length(citype)) {
        try({
          
          # taupriordensity <- function(t){dhalfcauchy(t,scale=1)}
          # ci <- 2
          # start <- Sys.time()
          bayesma <- bayesmeta(effs,
                             mu.prior.mean=0, mu.prior.sd=sqrt(10),
                             tau.prior=taupriordensity, 
                             interval.type=citype[ci],
                             seed=1)
          # print(bayesma)
          # Forestplots
          # forestplot(bayesma, exponentiate=TRUE, xlog=TRUE, zero=1, digits=2)
          # forestplot(bayesma,  zero=0, digits=2)
          # print( Sys.time() - start)
          #model-based estimates
          #t(exp(bayesma$theta))[,c(4,7, 8)]
          
          #posterior weights
          #bayesma$weights
          #bayesma$weights.theta
          
        k <- bayesma$k

        sigma2hat <- bayesma$summary[2,1]
        sigma2hatlo <- bayesma$summary[5, 1]
        sigma2hatup <- bayesma$summary[6, 1]
        
        #median
        esthatlo <- exp(bayesma$summary[5, 2]) 
        esthatup <- exp(bayesma$summary[6, 2])
        
        tcci <- citype[ci]
        
        for (t in 1:2) {
          if (t == 1) {
            esthat <- exp(bayesma$summary[3,2]) #mean
            parameter <- "mean-orout"
            
          } else {
            esthat <- exp(bayesma$summary[2,2]) #median
            parameter <- "median-orout"
          }
          
          line <-  paste(sim, ";",
                         paste(package, "-", inference, "-",  model,  "-", tcci,  "-", NA,  "-", loop, "-", NA,  "-", NA, sep=""), ";",
                         parameter, ";",
                         esthat, ";", esthatlo, ";", esthatup, ";",
                         tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
                         sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
                         IC, ";",
                         mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
                         k, ";", nstudies, ";", studysize, 
                         sep="")
          
          
          write(line, file=outfile, append = TRUE)
          
        }
        run = run + 1
        print(run)
        }, silent = TRUE)
      }
  }

#====================================metaplus===============
suppressWarnings(library(metaplus))
package <- "metaplus"
inference <- "frequentist"
model <- "PL"
parameter <- "mean-orout"
run=0

typerandom <- c("normal", "t-dist", "mixture")


  for (dist in 1:length(typerandom)) {
    
    tau2hat <- "."
    tau2hatlo <- "." 
    tau2hatup <- "." 
    
    sigma2hat <- "."
    sigma2hatlo <- "." 
    sigma2hatup <- "." 
    
    try({
        # dist <- 1
        # start <- Sys.time()
        plus <- metaplus(yi, 
               sqrt(vi), 
               mods = NULL, 
               random = typerandom[dist], 
               plotci = FALSE, 
               justfit = FALSE, 
               data=effs,
               slab = slab)
      # summary(plus) 
      # plot(plus, transf=exp,  refline=1, header=TRUE)
      # plot(plus,  refline=0, header=TRUE)
      # print( Sys.time() - start)
      
      esthat <- exp(plus$results[1, 1])
      esthatlo <-exp(plus$results[1, 2])
      esthatup <-exp(plus$results[1, 3])
      
      sigma2hat <- plus$results[2, 1]
      sigma2hatlo <- plus$results[2, 2]
      sigma2hatup <- plus$results[2, 3]
      
      disti <- typerandom[dist]
      k <- K
      
      line <-  paste(sim, ";",
                     paste(package, "-", inference, "-",  model,  "-", disti,  "-", NA,  "-", NA, "-", NA, "-", NA,  "-", NA, sep=""), ";",
                     parameter, ";",
                     esthat, ";", esthatlo, ";", esthatup, ";",
                     tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
                     sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";",
                     IC, ";",
                     mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
                     k, ";", nstudies, ";", studysize, 
                     sep="")
      
      
      write(line, file=outfile, append = TRUE)
      run = run + 1
      print(run)
    }, silent= TRUE)
  }

#plot(plus)

#======================metastan===============
suppressWarnings(library(MetaStan))
params <- c("Smith", "Higgins")
intervals <- c("shortest", "central")
taupriordist <- c("half-normal", "half-cauchy", "uniform")
priordist <- c("normal", "cauchy", "uniform")
reffs <- c(TRUE, FALSE)
package <- "metastan"
inference <- "bayesian"
run = 0

dat_MetaStan <- create_MetaStan_dat(dat = data.frame(r1=event0, 
                                                     n1=total0, 
                                                     r2=event1, 
                                                     n2=total1),
                                    armVars = c(responders = "r", 
                                                sampleSize = "n"))

for (param in 1:2) {
    for (reff in 1:2) {
      
      if (reff ==  1) {
        model = "RE"
      } else {
        model = "FE"
      }
      
      if (reff == 1) {
        loops <- 3
      } else {
        loops <- 1
      }
      
      for (loop in 1:loops) {
        for (interval in 1:2) {
          
          try({
            tau2hat <- "."
            tau2hatlo <- "." 
            tau2hatup <- "." 
            
            sigma2hat <- "."
            sigma2hatlo <- "." 
            sigma2hatup <- "." 
            
             # loop <- 2
             # interval <- 1
             # reff <- 1
             # param <- 2
             # start <- Sys.time()
             # dat_MetaStan <- create_MetaStan_dat(dat = data.frame(r1=event0, 
             #                                                      n1=total0, 
             #                                                      r2=event1, 
             #                                                      n2=total1),
             #                                     armVars = c(responders = "r", 
             #                                                 sampleSize = "n"))
              stanfit  <- meta_stan(data = dat_MetaStan,
                                  likelihood = "binomial",
                                  mu_prior = c(0, 10),
                                  theta_prior = c(0, 10),
                                  tau_prior = 0.5,
                                  tau_prior_dist = taupriordist[loop],
                                  chains = 3,
                                  param = params[param],
                                  interval.type = intervals[interval],
                                  re = reffs[reff],
                                  seed=1)
              
              #Poisson-log & normal-identity, undocumented and no examples
              # library("shinystan")
              # 
              # stanfitshiny = as.shinystan(stanfit$fit)
              # launch_shinystan(stanfitshiny)

              # print(stanfit)  
              # forest_plot(stanfit, labels=slab)
              # print( Sys.time() - start)
              
            intervali <- intervals[interval]
            parami <- params[param]
            
            for (t in 1:2) {
              if (t==1) {
                parameter <- "mean-orout"
                esthat <- exp(stanfit$fit_sum["theta",1]) #mean
              } else {
                parameter <- "median-orout"
                esthat <- exp(stanfit$fit_sum["theta",6]) #median
              }
              
              if (reff == 1) {
                priori <- priordist[loop]
                
                sigma2hat <-  (stanfit$fit_sum["tau[1]",6])^2
                sigma2hatlo <- (stanfit$fit_sum["tau[1]",4])^2
                sigma2hatup <- (stanfit$fit_sum["tau[1]",8])^2
                
              } else {
                priori <- NA
                sigma2hat <- 0
                sigma2hatlo <- 0
                sigma2hatup <- 0
              }
              
              esthatlo <- exp(stanfit$fit_sum["theta",4])
              esthatup <- exp(stanfit$fit_sum["theta",8])
             
              line <-  paste(sim, ";",
                             paste(package, "-", inference, "-",  model,  "-", intervali,  "-", priori,  "-", parami, "-", loop,"-", NA,  "-", NA, sep=""), ";",
                             parameter, ";",
                             esthat, ";", esthatlo, ";", esthatup, ";",
                             tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
                             sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
                             IC, ";",
                             mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
                             nstudies, ";", nstudies, ";", studysize, 
                             sep="")
              
              
              write(line, file=outfile, append = TRUE)
            }
            run = run + 1
            print(run)
          }, silent= TRUE)
        
        }
    }
  }
}

#=================================================meta===============
suppressWarnings(library(meta))
#methods <- c("Inverse", "Peto", "SSW", "MH", "LRP", "GLMM")
methods <- c("SSW", "MH")

# While the Mantel-Haenszel and Peto method are defined under the common effect model, 
# random effects variants based on these methods are also implemented in metabin. 
# random effects meta-analysis using the Mantel-Haenszel or inverse variance method are typically very similar. 
#typecc <- c("", "only0", "if0all", "all")
typecc <- c("")
#taumethods <- c("REML", "PM", "DL", "ML", "HS", "SJ", "HE", "EB")
taumethods <- c("REML", "PM", "HS", "SJ", "HE", "EB")
#cimethods <- c("classic", "HK", "KR")
cimethods <- c("HK", "KR")
taucimethods <- c("")
#taucimethods <- c("J", "BJ", "QP", "PL", "")
#adhoc_hakn_cis <- c("", "se", "ci", "IQWiG6")
adhoc_hakn_cis <- c("")
package <- "meta"
inference <- "frequentist"
run = 0

for (method in 1:length(methods)) {
    for (taumethod in 1:length(taumethods)) {
      for (ci in 1:length(cimethods)){
        for (cimethod in 1:length(taucimethods)){
          for (loop in 1:length(adhoc_hakn_cis)) {
            try({
              tau2hat <- "."
              tau2hatlo <- "." 
              tau2hatup <- "." 
              
              sigma2hat <- "."
              sigma2hatlo <- "." 
              sigma2hatup <- "." 
              
              # method <- 2
              # taumethod <- 2
              # cimethod <- 1
              #  
              # start <- Sys.time()
              
              res <- metabin(event.c=event0, 
                            n.c=total0, 
                            event.e=event1, 
                            n.e=total1,
                            sm = "OR",  
                            method = methods[method],
                            method.tau = taumethods[taumethod],
                            method.random.ci = cimethods[cimethod],
                            studlab = slab)
              
              # summary(res)
              # forest(res, prediction = TRUE)
              # print( Sys.time() - start)
              
              k <- res$k

              for (m in 1:2) {
                
                if (m == 1) {
                  model <- "common" 
                  esthat <- exp(res$TE.common) 
                  esthatlo <- exp(res$lower.common) 
                  esthatup <- exp(res$upper.common)
                  
                  sigma2hat <- 0
                  sigma2hatlo <- 0
                  sigma2hatup <- 0
                  
                } else {
                  model <- "random" 
                  esthat <- exp(res$TE.random) 
                  esthatlo <- exp(res$lower.random) 
                  esthatup <- exp(res$upper.random)
                  
                  sigma2hat <- res$tau2
                  sigma2hatlo <- res$lower.tau2
                  sigma2hatup <- res$upper.tau2
                }
              
                for (t in 1:2) {
                  if (t == 1) {
                    parameter = "mean-orout" 
                  } else {
                    parameter = "median-orout"
                  }
                  
                  methodi <- methods[method]
                  cci <- typecc[cc]
                  taumethodi <- taumethods[taumethod]
                  cii <- cimethods[ci]
                  cimethodi <- taucimethods[cimethod]
                  hakni <- adhoc_hakn_cis[loop]
                  
                  line <-  paste(sim, ";",
                                 paste(package, "-", inference, "-",  model,  "-", methodi,  "-", cci,  "-", taumethodi, "-", cii, "-", cimethodi, "-", hakni, sep=""), ";",
                                 parameter, ";",
                                 esthat, ";", esthatlo, ";", esthatup, ";",
                                 tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
                                 sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
                                 IC, ";",
                                 mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
                                 k, ";", nstudies, ";", studysize, 
                                 sep="")
                  
                  
                  write(line, file=outfile, append = TRUE)
                }
              }
              
              run = run + 1
              print(run)
            }, silent= TRUE)
            
          }
        }
      }
    }
}

#=================================randmeta===============
#Only works for more than 5 studies
# typecc <- c("none")
# suppressWarnings(library(RandMeta))
# models <- c("DL", "wilcox", "wang", "median")
# package <- "randmeta"
# inference <- "frequentist"
# parameter <- "mean-orout"
# run = 0
# 
# 
#   for (model in 1:length(models)) {
#     try({
#       tau2hat <- "."
#       tau2hatlo <- "." 
#       tau2hatup <- "." 
#       
#       sigma2hat <- "."
#       sigma2hatlo <- "." 
#       sigma2hatup <- "." 
#       
#       fit=random.meta(effs$yi, effs$vi, type=models[model])
#       
#       esthat <- exp(fit$theta)
#       esthatlo <- exp(fit$ci95[1])
#       esthatup <- exp(fit$ci95[2])
#       
#       sigma2hat <- NA
#       sigma2hatlo <- NA
#       sigma2hatup <- NA
#       
#       tau2hat <- NA
#       tau2hatlo <- NA
#       tausigma2hatup <- NA
#       
#       modeli <- models[model]
#       k <- K
# 
#       line <-  paste(sim, ";",
#                      paste(package, "-", inference, "-",  modeli,  "-", NA,  "-", NA, "-", NA, "-", NA, "-", NA, sep=""), ";",
#                      parameter, ";",
#                      esthat, ";", esthatlo, ";", esthatup, ";",
#                      tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
#                      sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
#                      IC, ";",
#                      mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
#                      k, ";", nstudies, ";", studysize, 
#                      sep="")
#       
#       
#       write(line, file=outfile, append = TRUE)
#       
#       run = run + 1
#       print(run)
#     }, silent= TRUE)
#   }

#=================================rema===============
# library(rema)
# package <- "rema"
# inference <- "frequentist"
# parameter <- "mean-orout"
# 
# try({
#   fit <- rema(trt.events = event1,
#             trt.total = total1,
#             ctrl.events = event0,
#             ctrl.total = total0,
#             distr = FALSE)
#   
#   esthat <- fit$TE 
#   esthatlo <- fit$CI[1] 
#   esthatup <-fit$CI[2]
#   
#   sigma2hat <- NA
#   sigma2hatlo <- NA
#   sigma2hatup <- NA
#   
#   tau2hat <- NA
#   tau2hatlo <- NA
#   tausigma2hatup <- NA
#   
#   line <-  paste(sim, ";",
#                  paste(package, "-", inference, "-",  "exact",  "-", NA,  "-", NA,  "-", NA, "-", NA, "-", NA, "-", NA, sep=""), ";",
#                  parameter, ";",
#                  esthat, ";", esthatlo, ";", esthatup, ";",
#                  tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
#                  sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";",
#                  IC, ";",
#                  mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
#                  k, ";", studysize, 
#                  sep="")
#   
#   
#   write(line, file=outfile, append = TRUE)
# }, silent = TRUE)

#================================================metaBMA===============
# suppressWarnings(library(metaBMA))
# package <- "metabma"
# inference <- "bayesian"
# model <- "IV"
# method <- NULL
# optssummarize <-  c("integrate", "stan")
# optslogml <-  c("integrate", "stan")
# taudist <- c("cauchy", "norm", "invgamma", "invgamma")
# parm1 <- c(0, 0, 0.01, 0.001)
# parm2 <- c(1, 1, 0.01, 0.001)
# run <- 0
# 
# for (s in 1:2) {
#   for (l in 1:2) {
#     
#     logmli <- optslogml[l]
#     summi <- optssummarize[s]
#     model <- "iv"
#     
#     tau2hat <- "."
#     tau2hatlo <- "." 
#     tau2hatup <- "." 
#     
#     sigma2hat <- "."
#     sigma2hatlo <- "." 
#     sigma2hatup <- "." 
#     
#     try({
#       mf <- meta_fixed(effs$yi, 
#                        sqrt(effs$vi), 
#                        slab,
#                        d = prior("norm", c(mean = 0, sd = sqrt(10))),
#                        logml = optssummarize[s],
#                        summarize = optslogml[l],
#                        chains=3)
#       
#       esthatlo <- exp(mf$estimates[1, 3])
#       esthatup <- exp(mf$estimates[1,5])
#       sigma2hat <- 0
#       k <- mf$data$N
#       
#       for (t in 1:2) {
#         if (t == 1) {
#           esthat <- exp(mf$estimates[1, 1]) #mean 
#           parameter <- "mean-orout"
#         } else {
#           esthat <- exp(mf$estimates[1, 4]) #median
#           parameter <- "median-orout"
#         }
#         
#         line <-  paste(sim, ";",
#                        paste(package, "-", inference, "-",  model,  "-", logmli,  "-", summi,  "-",   NA,"-", NA,  "-", NA, sep=""), ";",
#                        parameter, ";",
#                        esthat, ";", esthatlo, ";", esthatup, ";",
#                        tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
#                        sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
#                        IC, ";",
#                        mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
#                        k, ";", nstudies, ";", studysize, 
#                        sep="")
#         
#         
#         write(line, file=outfile, append = TRUE)
#       }
#       run = run + 1
#       print(run)
#     }, silent = FALSE)  
#     
#     for (loop in 1:4) {
#       tau2hat <- "."
#       tau2hatlo <- "." 
#       tau2hatup <- "." 
#       
#       sigma2hat <- "."
#       sigma2hatlo <- "." 
#       sigma2hatup <- "." 
#       
#       try({
#         loop <- 1
#         s <-  2
#         l <- 2
#         me <- meta_random(y=effs$yi, 
#                           SE=sqrt(effs$vi), 
#                           d = prior("norm", c(mean = 0, sd = sqrt(10))),
#                           tau = prior(taudist[loop], 
#                                       c(parm1[loop], 
#                                         parm2[loop]), lower=0),
#                           logml = optssummarize[s],
#                           summarize = optslogml[l],
#                           chains=3)
#         
#         plot_forest(me)
#         
#         
#         esthatlo <- exp(me$estimates[1, 3])
#         esthatup <- exp(me$estimates[1,5])
#         
#         sigma2hat <- me$estimates[2, 4] #median
#         sigma2hatlo <- me$estimates[2, 3]
#         sigma2hatup <- me$estimates[2, 5]
#         model <- "DL"
#         k <- me$data$N
#         
#         for (t in 1:2) {
#           if (t == 1) {
#             esthat <- exp(me$estimates[1, 1]) #mean 
#             parameter <- "mean-orout"
#           } else {
#             esthat <- exp(me$estimates[1, 4]) #median
#             parameter <- "median-orout"
#           }
#           
#           line <-  paste(sim, ";",
#                          paste(package, "-", inference, "-",  model,  "-", logmli,  "-", summi,  "-",  loop, sep=""), ";",
#                          parameter, ";",
#                          esthat, ";", esthatlo, ";", esthatup, ";",
#                          tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
#                          sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
#                          IC, ";",
#                          mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
#                          k, ";", nstudies, ";", studysize, 
#                          sep="")
#           
#           
#           write(line, file=outfile, append = TRUE)
#         }
#         run = run + 1
#         print(run)
#       }, silent= TRUE)
#     }
#   }
# }

#======================================metasem===============
# suppressWarnings(library(metaSEM))
# intervals = c("LB", "z")
# constraints <- c(NULL, 0)
# parameter <- "mean-orout"
# package <- "metasem"
# inference <- "frequentist"
# run = 0
# 
# for (loop in 1:2) {
#   if (loop == 1) {
#     model <- "IV"  
#   } else {
#     model <- "DL"  
#   }
#   for (t in 1:2) {
#     tau2hat <- "."
#     tau2hatlo <- "." 
#     tau2hatup <- "." 
#     
#     sigma2hat <- "."
#     sigma2hatlo <- "." 
#     sigma2hatup <- "." 
#     
#     try({
#       # t <- 1
#       # loop <- 1
#       fitsem <- meta(y=matrix(effs$yi), 
#                      v=matrix(effs$vi),
#                      RE.constraints=constraints[t],
#                      intervals.type = intervals[loop])
#       
#       fitsem <- meta(y=matrix(effs$yi), 
#                      v=matrix(effs$vi),
#                      RE.constraints=NULL,
#                      intervals.type = "LB")
#       
#       summsem <- summary(fitsem)
#       
#       #metafor
#       #forest(rma(yi=matrix(effs$yi), vi=matrix(effs$vi), slab=slab, data=effs), refline=1, addfit=FALSE, header=TRUE, transf=exp, ylim=c(nstudies+3, -1))
#       
#       esthat <- exp(summsem$coefficients[1,1])
#       esthatlo <- exp(summsem$coefficients[1,3])
#       esthatup <- exp(summsem$coefficients[1,4])
#       
#       #abline(h=0)
#       #addpoly(x=esthat, ci.lb=esthatlo, ci.ub=esthatup,  mlab="Summary", rows=-1, transf=identity)
#       
#       if (t == 1) {
#         sigma2hat <- summsem$coefficients[2,1]
#         sigma2hatlo <- summsem$coefficients[2,3]
#         sigma2hatup <- summsem$coefficients[2,4]
#       } else {
#         sigma2hat <- 0
#         sigma2hatlo <- 0
#         sigma2hatup <- 0
#       }
#       
#       intervali <- intervals[loop]
#       contraint <- constraints[t]
#       k <- K
#       
#       line <-  paste(sim, ";",
#                      paste(package, "-", inference, "-",  model,  "-", intervali,  "-", contraint,  "-", NA, "-", NA, "-", NA,  "-", NA, sep=""), ";",
#                      parameter, ";",
#                      esthat, ";", esthatlo, ";", esthatup, ";",
#                      tau2hat, ";", tau2hatlo, ";", tau2hatup, ";",
#                      sigma2hat, ";",sigma2hatlo, ";",sigma2hatup, ";", 
#                      IC, ";",
#                      mu0true, ";", ortrue, ";", tausqtrue, ";", sigmasqtrue, ";",
#                      k, ";", nstudies, ";", studysize, 
#                      sep="")
#       
#       
#       write(line, file=outfile, append = TRUE)
#       run = run + 1
#       print(run)
#     }, silent= TRUE)
#   }
# }
