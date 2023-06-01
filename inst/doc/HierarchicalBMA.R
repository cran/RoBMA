## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/HierarchicalBMA/fit.RDS")
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
fit.0             <- readRDS(file = "../models/HierarchicalBMA/fit.0.RDS")
fit               <- readRDS(file = "../models/HierarchicalBMA/fit.RDS")
fit_BMA           <- readRDS(file = "../models/HierarchicalBMA/fit_BMA.RDS")
hierarchical_test <- readRDS(file = "../models/HierarchicalBMA/hierarchical_test.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
#  # R package version updating
#  library(RoBMA)
#  
#  data("dat.konstantopoulos2011", package = "metadat")
#  dat <- dat.konstantopoulos2011
#  
#  fit.0 <- RoBMA(d = dat$yi, v = dat$vi,
#               priors_effect_null        = NULL,
#               priors_heterogeneity_null = NULL,
#               priors_bias               = NULL,
#               parallel = TRUE, seed = 1)
#  
#  fit <- RoBMA(d = dat$yi, v = dat$vi, study_ids = dat$district,
#               priors_effect_null        = NULL,
#               priors_heterogeneity_null = NULL,
#               priors_bias               = NULL,
#               parallel = TRUE, seed = 1)
#  
#  fit_BMA <- RoBMA(d = dat$yi, v = dat$vi, study_ids = dat$district,
#                   priors_bias = NULL,
#                   parallel = TRUE, seed = 1)
#  
#  hierarchical_test <- RoBMA(d = dat$yi, v = dat$vi, study_ids = dat$district,
#                             priors_heterogeneity_null = NULL,
#                             priors_hierarchical_null = prior(distribution = "spike", parameters = list("location" = 0)),
#                             priors_bias = NULL,
#                             parallel = TRUE, seed = 1)
#  
#  saveRDS(fit.0,             file = "../models/HierarchicalBMA/fit.0.RDS",             compress = "xz")
#  saveRDS(fit,               file = "../models/HierarchicalBMA/fit.RDS",               compress = "xz")
#  saveRDS(fit_BMA,           file = "../models/HierarchicalBMA/fit_BMA.RDS",           compress = "xz")
#  saveRDS(hierarchical_test, file = "../models/HierarchicalBMA/hierarchical_test.RDS", compress = "xz")

## -----------------------------------------------------------------------------
data("dat.konstantopoulos2011", package = "metadat")
dat <- dat.konstantopoulos2011

head(dat)

## -----------------------------------------------------------------------------
fit_metafor.0 <- metafor::rma(yi = yi, vi = vi, data = dat)
fit_metafor.0

## -----------------------------------------------------------------------------
fit_metafor <- metafor::rma.mv(yi, vi, random = ~ school | district, data = dat)
fit_metafor

## -----------------------------------------------------------------------------
summary(fit.0, type = "individual")

## -----------------------------------------------------------------------------
summary(fit, type = "individual")

## ----fig_rho, dpi = 300, fig.width = 4, fig.height = 3, out.width = "50%", fig.align = "center"----
par(mar = c(2, 4, 0, 0))
plot(fit, parameter = "rho", prior = TRUE)

## -----------------------------------------------------------------------------
summary(fit_BMA)

## -----------------------------------------------------------------------------
summary(hierarchical_test)

