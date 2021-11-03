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
fit_BMA   <- readRDS(file = "../models/MedicineBMA/fit_BMA.RDS")
fit_BMAb  <- readRDS(file = "../models/MedicineBMA/fit_BMAb.RDS")
fit_RoBMA <- readRDS(file = "../models/MedicineBMA/fit_RoBMA.RDS")

## -----------------------------------------------------------------------------
library(RoBMA)

data("Poulsen2006", package = "RoBMA")
Poulsen2006

## -----------------------------------------------------------------------------
summary(fit_BMA, conditional = TRUE)

## ----fig_mu_BMA, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit_BMA, parameter = "mu", prior = TRUE)

## ----fig_mu_BMA_cond, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit_BMA, parameter = "mu", prior = TRUE, conditional = TRUE)

## ----fig_models, dpi = 300, fig.height = 4.5, fig.width = 7, out.width = '75%', fig.align = "center"----
plot_models(fit_BMA)

## ----fig_forest, dpi = 300, fig.height = 4.5, fig.width = 7, out.width = '75%', fig.align = "center"----
forest(fit_BMA, conditional = TRUE)

## -----------------------------------------------------------------------------
summary(fit_BMAb, conditional = TRUE)

## -----------------------------------------------------------------------------
summary(fit_RoBMA, conditional = TRUE)

## ----fig_mu_RoBMA_cond, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit_RoBMA, parameter = "mu", prior = TRUE, conditional = TRUE)

