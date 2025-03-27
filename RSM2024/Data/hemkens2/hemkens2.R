library(haven)
library(ggplot2)

library(dplyr)
library(tidyr)

library(forcats)
library(simhelpers)
library(gghalves)

#=================hemkens2 ====================
review <- "hemkens2"
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

#===================================All=====================================

sim$id <- with(sim, paste(package, Env, Dist, param, Covariance, Weighting, Sigmethod, Prior, inf, sep='-'))

orange = "IND"  

Env <- "Stata"

package <- "meta"
package <- "metan"

Env <- "R"

package <- "metafor"
package <- "metastan"
package <- "meta"
package <- "bayesmeta"
package <- "metaplus"

#-------------------------------
simor <- sim[(sim$package==package & 
                sim$Env==Env & ((sim$stat == 'median-orout' & sim$inf=="B") | 
                                  (sim$stat == 'mean-orout' & sim$inf=="F"))) |
               
               (sim$package=="metapreg" & 
                  ((sim$stat == 'median-orout' & sim$inf=="B") | 
                     (sim$stat == 'mean-orout' & sim$inf=="F"))  & 
                  sim$link == "logit" & 
                  sim$Covariance == orange  & 
                  sim$Design == "comparative")
             ,]

source('C:/DATA/WIV/Projects/GitHub/Metapreg/Scripts/perfomance-stats.R')

min <- lo
max <- up
#=================================conditional statistics=============================
sim$conditional <- with(sim,ifelse((package == "metapreg" & link == "logit" & 
                                      (Covariance == "IND"  & 
                                         Design == "comparative") & 
                                      Dist == "BN"), 
                                   1,0))

sim$conditional <- with(sim,ifelse(inference == "bayesian" & package != "metapreg" & Dist != "B", 1, conditional))
sim$conditional <- with(sim,ifelse(Sigmethod == "SJ" & Weighting == "SSW" & package == "meta" & Env == "R", 1, conditional))
sim$conditional <- with(sim,ifelse(slope == "common", 0, conditional))
sim$conditional <- with(sim,ifelse((Sigmethod == "sj2s" | CImethod == "ivhet" | Sigmethod == "mp") & package == "metan"  , 1, conditional))
sim$conditional <- with(sim,ifelse(Sigmethod == "sj" & package == "meta" & Env == "Stata", 1, conditional))
sim$conditional <- with(sim,ifelse(package %in% c("metasem", "randmeta", "metaplus") == TRUE, 1, conditional))
sim$conditional <- with(sim,ifelse(package == "metafor" & Dist == "BN" & (Covariance == "IND" ), 1, conditional))


simor <- data.frame(sim[(sim$stat %in% c('mean-orout', 'median-orout')==TRUE) & 
                          (sim$conditional==1),])

simor$id <- with(simor, paste(package, Env, Dist, param, Sigmethod, Prior, inf, sep='-'))

# simor$esthat <- log(simor$esthat)
# simor$esthatlo <- log(simor$esthatlo)
# simor$esthatup <- log(simor$esthatup)
# simor$trueor <- log(simor$trueor)

coverage <- simor %>%
  mutate(params = model) %>%
  group_by(model,  id, stat, package, Dist, param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link) %>%
  group_modify(~ calc_coverage(.x, lower_bound = esthatlo, upper_bound = esthatup, true_param = trueor))  %>%
  as.data.frame()

rbias <- simor %>%
  mutate(params = model) %>%
  group_by(model,  id,stat, package, Dist, param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link) %>%
  group_modify(~ calc_relative(.x, estimate = esthat, true_param = trueor))  %>%
  as.data.frame()

coverage$coverage <- format(coverage$coverage, digits=2)
coverage$covtext <- with(coverage, ifelse(CImethod != '', paste(coverage, '(', CImethod, ')', sep=''), coverage))
rbias$rel_bias <- format(rbias$rel_bias, digits=2)
rbias$rel_mse <- format(rbias$rel_mse, digits=2)

# trueor <- log(min(sim$trueor))
# formatted_trueor <- format(log(min(sim$trueor)), digits=3)

trueor <- (min(sim$trueor))
formatted_trueor <- format((min(sim$trueor)), digits=2)

max <- ceiling(max(simor$esthat))
min <- floor(min(simor$esthat))

coverage$left1 <- (trueor - 5 )
coverage$left2 <- (trueor - 6)
coverage$left3 <- (trueor -7)
coverage$left4 <- (trueor -8)
coverage$left5 <- (trueor -9)
coverage$left6 <- (trueor -10)
coverage$left7 <- (trueor -11)
coverage$left8 <- (trueor -12)

rbias$right2 <- (trueor -1 )
rbias$right1 <- (trueor + 0.25 )
coverage$right3 <- (trueor + 1.75 )
coverage$right4 <- (trueor + 3.25)

coverage$right5 <- (trueor + 4.5)
coverage$right6 <- (trueor + 5.75)

# 
# coverage$left1 <- (min - 0.25 )
# coverage$left2 <- (min - 0.5)
# coverage$left3 <- (min -1.25)
# coverage$left4 <- (min -1.75)
# coverage$left5 <- (min -2.25)
# coverage$left6 <- (min -2.75)
# coverage$left7 <- (min -3.25)
# coverage$left8 <- (min -3.75)
# 
# rbias$right2 <- (max -1.25 )
# rbias$right1 <- (max - 0.75 )
# coverage$right3 <- (max - 0.25 )
# coverage$right4 <- (max + 0.25)
# 
# coverage$right5 <- (max + 0.75)
# coverage$right6 <- (max + 1.25)

minX <- min(coverage$left8)
maxX <- min(coverage$right4)

# sort dataframe by relbiase
order_median <- unique(rbias[rbias$Stat=="Median",] %>% arrange(desc(rel_bias)) %>% pull(id))
order_mean <- unique(rbias[rbias$Stat=="Mean",] %>% arrange(desc(rel_bias)) %>% pull(id))


# Plot
med <-  ggplot() +
  geom_half_violin(data = simor[simor$Stat=="Median",],
                   aes(x=factor(id, levels = order_median), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="gray", fill="gray") +
  
  geom_half_violin(data= simor[simor$package=="metapreg" & simor$Stat=="Median",], 
                   aes(x=factor(id, levels = order_median), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="orange", fill="orange") +
  
  
  geom_boxplot(data=simor[simor$Stat=="Median" & (simor$CImethod %in% c("eti", "HK")==TRUE),],  
               aes(x=factor(id, levels = order_median), y=esthat),
               width=0.1, color="black", alpha=0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left5, label=Env), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left7, label=package), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left6, label=inf), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left8, label=Dist), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left4, label=param), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),],
            aes(x=factor(id, levels = order_median), y=left3, label=Sigmethod), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=left1, label=Prior), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$Stat=="Median" & (rbias$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=right2, label=rel_bias), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$Stat=="Median" & (rbias$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=right1, label=rel_mse), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("eti", "HK")==TRUE),], 
            aes(x=factor(id, levels = order_median), y=right3, label=covtext), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Median" & (coverage$CImethod %in% c("hpd", 'KR')==TRUE),], 
            aes(x=factor(id, levels = order_median), y=right4, label=covtext), nudge_x = 0.2) +
  
  theme(legend.position="none") +
  coord_flip() +
  scale_y_continuous(position = 'right',
                     limits=c(-5, 20),
                     breaks=c(-5, 
                              min(coverage$left1), 
                              min(coverage$left3),
                              min(coverage$left4),
                              min(coverage$left5),
                              min(coverage$left6),
                              min(coverage$left7),
                              min(coverage$left8),
                              trueor,
                              min(rbias$right1),
                              min(rbias$right2),
                              min(coverage$right3),
                              20),
                     expand = c(0.04, 0.02),
                     labels = c("", "Prior", "Sig", "Param", "Env", "Inf",  "Package", 'Dist',  "True OR", "RMSE", "RBias", "Coverage", ""),
                     sec.axis = sec_axis(~.,breaks=c(trueor),
                                         labels=c(formatted_trueor))) +
  
  
  geom_hline(yintercept = trueor, linewidth=0.5, colour='black')   +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_line(color='black'),
        axis.text.x=element_text(size=12, face='bold'),
        panel.grid.major.y = element_line(color = 'gray'),
        panel.background = element_rect(fill = "white")) +
  xlab("") + ylab("") 

med

ggsave('Dist-Sim-Median.png',
       path = graphpath,
       width = 12,
       height = 8,
       units = 'in',
       dpi = 300)

ggsave('Dist-Sim-Median.pdf',
       path = graphpath,
       width = 12,
       height = 8,
       units = 'in',
       dpi = 300)

#-----------------------------------------
maxX <- min(coverage$right6)
# Plot
mean <-  ggplot() +
  
  geom_half_violin(data = simor[simor$Stat=="Mean",],
                   aes(x=factor(id, levels = order_mean), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="gray", fill="gray") +
  
  geom_half_violin(data = simor[simor$Stat=="Mean" & simor$package == "metapreg",],
                   aes(x=factor(id, levels = order_mean), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="orange", fill="orange") +
  
  geom_boxplot(data=simor[simor$Stat=="Mean",],  
               aes(x=factor(id, levels = order_mean), y=esthat),
               width=0.1, color="black", alpha=0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),], 
            aes(x=factor(id, levels = order_mean), y=left5, label=Env), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),], 
            aes(x=factor(id, levels = order_mean), y=left7, label=package), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),], 
            aes(x=factor(id, levels = order_mean), y=left6, label=inf), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),], 
            aes(x=factor(id, levels = order_mean), y=left8, label=Dist), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),],
            aes(x=factor(id, levels = order_mean), y=left4, label=param), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),],
            aes(x=factor(id, levels = order_mean), y=left3, label=Sigmethod), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & !duplicated(coverage$id),],
            aes(x=factor(id, levels = order_mean), y=left1, label=Prior), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$Stat=="Mean" & !duplicated(rbias$id),], 
            aes(x=factor(id, levels = order_mean), y=right2, label=rel_bias), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$Stat=="Mean" & !duplicated(rbias$id),],  
            aes(x=factor(id, levels = order_mean), y=right1, label=rel_mse), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & 
                            ((coverage$CImethod %in% c("eti", "HK",  'PL', 'LB', 't', 'ivhet')==TRUE) | coverage$CImethod ==""),], 
            aes(x=factor(id, levels = order_mean), y=right3, label=covtext), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & (coverage$CImethod %in% c("hpd", 'KR', 'z')==TRUE),], 
            aes(x=factor(id, levels = order_mean), y=right4, label=covtext), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & (coverage$CImethod %in% c("tkh")==TRUE),], 
            aes(x=factor(id, levels = order_mean), y=right5, label=covtext), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$Stat=="Mean" & (coverage$CImethod %in% c("tkht")==TRUE),], 
            aes(x=factor(id, levels = order_mean), y=right6, label=covtext), nudge_x = 0.2) +
  
  theme(legend.position="none") +
  coord_flip() +
  scale_y_continuous(position = 'right',
                     limits=c(minX, maxX),
                     breaks=c(minX, 
                              min(coverage$left1), 
                              min(coverage$left3),
                              min(coverage$left4),
                              min(coverage$left5),
                              min(coverage$left6),
                              min(coverage$left7),
                              min(coverage$left8),
                              trueor,
                              min(rbias$right1),
                              min(rbias$right2),
                              min(coverage$right3),
                              maxX),
                     expand = c(0.04, 0.01),
                     labels = c("", "Prior", "Sig", "Param", "Env", "Inf",  "Package", 'Dist',  "True LOR", "RMSE", "RBias", "Coverage", ""),
                     sec.axis = sec_axis(~.,breaks=c(trueor),
                                         labels=c(formatted_trueor))) +
  
  
  geom_hline(yintercept = trueor, linewidth=0.5, colour='black')   +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_line(color='black'),
        axis.text.x=element_text(size=12, face='bold'),
        panel.grid.major.y = element_line(color = 'gray'),
        panel.background = element_rect(fill = "white")) +
  xlab("") + ylab("") 

mean

ggsave('Dist-Sim-Mean.png',
       path = graphpath,
       width = 14,
       height = 8,
       dpi = 300)

ggsave('Dist-Sim-Mean.pdf',
       path = graphpath,
       width = 14,
       height = 8,
       dpi = 300)
