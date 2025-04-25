simor$truetotalvar = simor$truesigmasq + simor$truetausq
simor$null <- 1

# simor$esthat <- log(simor$esthat)
# simor$esthatlo <- log(simor$esthatlo)
# simor$esthatup <- log(simor$esthatup)
# simor$trueor <- log(simor$trueor)

power <- simor %>%
  mutate(params = model) %>%
  group_by(model, id, stat, Covariance, package, Dist,  Weighting, param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link, truetotalvar) %>%
  group_modify(~ calc_coverage(.x, lower_bound = esthatlo, upper_bound = esthatup, true_param = null))  %>%
  as.data.frame()

names(power)[18] <- 'power'

coverage <- simor %>%
  mutate(params = model) %>%
  group_by(model, id, stat, Covariance, package, Dist,  Weighting, Slope, param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link, truetotalvar) %>%
  group_modify(~ calc_coverage(.x, lower_bound = esthatlo, upper_bound = esthatup, true_param = trueor))  %>%
  as.data.frame()

lorcoverage <- simor %>%
  mutate(params = model) %>%
  group_by(model, id, stat, Covariance, package, Dist,  Weighting, Slope, param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link, truetotalvar) %>%
  group_modify(~ calc_coverage(.x, lower_bound = log(esthatlo), upper_bound = log(esthatup), true_param = log(trueor)))  %>%
  as.data.frame()

rbias <- simor %>%
  mutate(params = model) %>%
  group_by(model,  id,stat, Covariance, package, Dist,  Weighting, Slope,  param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link, truetotalvar) %>%
  group_modify(~ calc_relative(.x, estimate = esthat, true_param = trueor))  %>%
  as.data.frame()

abias <- simor[simor$esthat > 0,] %>%
  mutate(params = model) %>%
  group_by(model,  id,stat, Covariance, package, Dist,  Weighting, Slope,  param, CImethod, Sigmethod, Stat, Prior, inf, Measure, Env, Link, truetotalvar) %>%
  group_modify(~ calc_absolute(.x, estimate = log(esthat), true_param = log(trueor)))  %>%
  as.data.frame()

coverage$coverage <- format(coverage$coverage, digits=2)
coverage$width <- format(coverage$width, digits=2)
power$power <- format(power$power, digits=2)

coverage$covtext <- ifelse((coverage$CImethod != '' & coverage$K_coverage>0), 
                           paste(coverage$coverage, '(', coverage$width, ',', coverage$CImethod, ')', sep=''), 
                                          ifelse((coverage$CImethod == '' & coverage$K_coverage>0), 
                                                 paste(coverage$coverage, '(', coverage$width, ',', '-', ')', sep=''), 
                                                 ''))
abias$bias <- format(round(abias$bias, 2), nsmall=2)
abias$rmse <- format(round(abias$rmse, 2), nsmall=2)

rbias$rel_bias <- format(round(rbias$rel_bias, 2), nsmall=2)
rbias$rel_rmse <- format(round(rbias$rel_rmse, 2), nsmall=2)
rbias$rel_mse <- format(round(rbias$rel_mse, 2), nsmall=2)


max <- ceiling(max(simor$esthat))
min <- floor(min(simor$esthat))

trueor <- (min(sim$trueor))
formatted_trueor <- format((min(sim$trueor)), digits=3)

