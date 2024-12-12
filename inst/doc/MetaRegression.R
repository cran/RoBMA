## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/MetaRegression/fit_BMA.RDS")
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
fit_BMA   <- readRDS(file = "../models/MetaRegression/fit_BMA.RDS")
fit_RoBMA <- readRDS(file = "../models/MetaRegression/fit_RoBMA.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# library(RoBMA)
# 
# fit_BMA    <- NoBMA.reg(~ measure + age, data = Andrews2021, parallel = TRUE, seed = 1)
# fit_RoBMA  <- RoBMA.reg(~ measure + age, data = Andrews2021, parallel = TRUE, seed = 1, chains = 1)
# 
# saveRDS(fit_BMA,   file = "../models/MetaRegression/fit_BMA.RDS",   compress = "xz")
# saveRDS(fit_RoBMA, file = "../models/MetaRegression/fit_RoBMA.RDS", compress = "xz")

## -----------------------------------------------------------------------------
library(RoBMA)
data("Andrews2021", package = "RoBMA")
head(Andrews2021)

## -----------------------------------------------------------------------------
fit_rma <- metafor::rma(yi = r, sei = se, mods = ~ measure + age, data = Andrews2021)
fit_rma

## -----------------------------------------------------------------------------
emmeans::emmeans(metafor::emmprep(fit_rma), specs = "measure")

## -----------------------------------------------------------------------------
summary(fit_BMA, output_scale = "r")

## -----------------------------------------------------------------------------
marginal_summary(fit_BMA, output_scale = "r")

## ----fig_BMA, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
marginal_plot(fit_BMA, parameter = "measure", output_scale = "r", lwd = 2)

## -----------------------------------------------------------------------------
summary(fit_RoBMA, output_scale = "r")

## -----------------------------------------------------------------------------
marginal_summary(fit_RoBMA, output_scale = "r")

## ----fig_RoBMA, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
marginal_plot(fit_RoBMA, parameter = "measure", output_scale = "r", lwd = 2)

