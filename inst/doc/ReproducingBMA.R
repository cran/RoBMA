## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/ReproducingBMA/BMA_PowerPoseTest.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev      = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## ----include = FALSE----------------------------------------------------------
library(RoBMA)
# we pre-load the RoBMA models, the fitting time is around 2-5 minutes
fit_BMA_test   <- readRDS(file = "../models/ReproducingBMA/BMA_PowerPoseTest.RDS")
fit_BMA_est    <- readRDS(file = "../models/ReproducingBMA/BMA_PowerPoseEst.RDS")
fit_RoBMA_test <- readRDS(file = "../models/ReproducingBMA/PowerPoseTest.RDS")
fit_RoBMA_est  <- readRDS(file = "../models/ReproducingBMA/PowerPoseEst.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# 
# data("power_pose", package = "metaBMA")
# 
# fit_RoBMA_test <- RoBMA(d = power_pose$effectSize, se = power_pose$SE, study_names = power_pose$study,
#                         priors_effect  = prior(
#                           distribution = "cauchy",
#                           parameters = list(location = 0, scale = 1/sqrt(2)),
#                           truncation = list(0, Inf)),
#                         priors_heterogeneity = prior(
#                           distribution = "invgamma",
#                           parameters = list(shape = 1, scale = 0.15)),
#                         priors_bias = NULL,
#                         transformation = "cohens_d", seed = 1, parallel = TRUE)
# 
# fit_RoBMA_est  <- RoBMA(d = power_pose$effectSize, se = power_pose$SE, study_names = power_pose$study,
#                         priors_effect  = prior(
#                           distribution = "cauchy",
#                           parameters = list(location = 0, scale = 1/sqrt(2))),
#                         priors_heterogeneity = prior(
#                           distribution = "invgamma",
#                           parameters = list(shape = 1, scale = 0.15)),
#                         priors_bias = NULL,
#                         priors_effect_null = NULL,
#                         transformation = "cohens_d", seed = 2, parallel = TRUE)
# 
# saveRDS(fit_RoBMA_test, file = "../models/ReproducingBMA/PowerPoseTest.RDS")
# saveRDS(fit_RoBMA_est, file = "../models/ReproducingBMA/PowerPoseEst.RDS")

## -----------------------------------------------------------------------------
data("power_pose", package = "metaBMA")
power_pose[,c("study", "effectSize", "SE")]

## -----------------------------------------------------------------------------
fit_BMA_test$inclusion

round(fit_BMA_est$estimates,2)

## -----------------------------------------------------------------------------
summary(fit_RoBMA_test)

summary(fit_RoBMA_est)

## ----fig_mu_est, dpi = 300, fig.width = 4, fig.height = 4, out.width = "50%", fig.align = "center"----
plot(fit_RoBMA_est, parameter = "mu", prior = TRUE, xlim = c(-1, 1))

## ----fig_mu_test, dpi = 300, fig.width = 4, fig.height = 4, out.width = "50%", fig.align = "center"----
plot(fit_RoBMA_test, parameter = "mu", prior = TRUE, xlim = c(-.5, 1))

## ----fig_models, dpi = 300, fig.height = 4.5, fig.width = 7, out.width = '75%', fig.align = "center"----
plot_models(fit_RoBMA_est)

## ----fig_forest, dpi = 300, fig.height = 4.5, fig.width = 7, out.width = '75%', fig.align = "center"----
forest(fit_RoBMA_est)

