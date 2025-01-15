## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/FastRoBMA/fit_RoBMA.bridge.RDS")
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
fit_RoBMA.bridge  <- readRDS(file = "../models/FastRoBMA/fit_RoBMA.bridge.RDS")
fit_RoBMA.ss      <- readRDS(file = "../models/FastRoBMA/fit_RoBMA.ss.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# library(RoBMA)
# 
# t1_bridge <- Sys.time()
# fit_RoBMA.bridge    <- RoBMA.reg(
#   ~ measure + age, data = Andrews2021,
#   algorithm = "bridge",
#   chains = 3, sample = 5000, burnin = 2000, adapt = 500, parallel = TRUE
# )
# t2_bridge <- Sys.time()
# 
# t1_ss <- Sys.time()
# fit_RoBMA.ss  <- RoBMA.reg(
#   ~ measure + age, data = Andrews2021,
#   algorithm = "ss",
#   chains = 6, sample = 10000, burnin = 2500, adapt = 2500, parallel = TRUE
# )
# t2_ss <- Sys.time()
# 
# # save memory space
# fit_RoBMA.bridge <- RoBMA:::.remove_model_posteriors(fit_RoBMA.bridge)
# fit_RoBMA.bridge <- RoBMA:::.remove_model_margliks(fit_RoBMA.bridge)
# fit_RoBMA.ss$model$fit <- NULL
# 
# saveRDS(fit_RoBMA.bridge, file = "../models/FastRoBMA/fit_RoBMA.bridge.RDS",    compress = "xz")
# saveRDS(fit_RoBMA.ss,     file = "../models/FastRoBMA/fit_RoBMA.ss.RDS",        compress = "xz")

## -----------------------------------------------------------------------------
library(RoBMA)
data("Andrews2021", package = "RoBMA")
head(Andrews2021)

## -----------------------------------------------------------------------------
summary(fit_RoBMA.bridge)
summary(fit_RoBMA.ss)

## ----fig_effect, dpi = 300, fig.width = 12, fig.height = 6, out.width = "100%", fig.align = "center"----
par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
plot(fit_RoBMA.bridge, parameter = "mu", prior = TRUE, main = "Bridge Sampling")
plot(fit_RoBMA.ss,     parameter = "mu", prior = TRUE, main = "Spike-and-Slab")

## -----------------------------------------------------------------------------
marginal_summary(fit_RoBMA.bridge)
marginal_summary(fit_RoBMA.ss)

## ----fig_marginal, dpi = 300, fig.width = 12, fig.height = 6, out.width = "100%", fig.align = "center"----
par(mfrow = c(1, 2))
marginal_plot(fit_RoBMA.bridge, parameter = "measure", conditional = TRUE, prior = TRUE, xlim = c(-1, 1), ylim = c(0, 7), main = "Bridge Sampling")
marginal_plot(fit_RoBMA.ss,     parameter = "measure", conditional = TRUE, prior = TRUE, xlim = c(-1, 1), ylim = c(0, 7), main = "Spike-and-Slab")

