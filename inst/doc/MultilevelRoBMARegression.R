## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/MultilevelRoBMARegression/zfit_Havrankova2025.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev       = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## ----include = FALSE----------------------------------------------------------
library(RoBMA)
zfit_reg <- readRDS(file = "../models/MultilevelRoBMARegression/zfit_Havrankova2025.RDS")
fit_reg  <- zfit_reg
class(fit_reg) <- class(fit_reg)[!class(fit_reg) %in% "zcurve_RoBMA"]

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# data("Havrankova2025", package = "RoBMA")
# 
# # Prior scaling
# fit_fe        <- metafor::rma(yi = y, sei = se, data = Havrankova2025, method = "FE")
# unti_scale    <- fit_fe$se * sqrt(sum(Havrankova2025$N))
# prior_scale   <- unti_scale * 0.5
# 
# df_reg  <- data.frame(
#   y               = Havrankova2025$y,
#   se              = Havrankova2025$se,
#   facing_customer = Havrankova2025$facing_customer,
#   study_id        = Havrankova2025$study_id
# )
# 
# fit_reg <- RoBMA.reg(
#   ~ facing_customer,
#   study_ids      = df_reg$study_id,
#   data           = df_reg,
#   rescale_priors = prior_scale,
#   prior_scale = "none", transformation = "none",
#   algorithm = "ss", sample = 20000, burnin = 10000, adapt = 10000,
#   thin = 5, parallel = TRUE, autofit = FALSE, seed = 1)
# )
# 
# zfit_reg    <- as_zcurve(fit_reg)
# saveRDS(zfit_reg, file = "../models/MultilevelRoBMARegression/zfit_Havrankova2025.RDS", compress = "xz")

## -----------------------------------------------------------------------------
library(RoBMA)
data("Havrankova2025", package = "RoBMA")

fit_fe        <- metafor::rma(yi = y, sei = se, data = Havrankova2025, method = "FE")
unti_scale    <- fit_fe$se * sqrt(sum(Havrankova2025$N))
prior_scale   <- unti_scale * 0.5

## ----eval = FALSE-------------------------------------------------------------
# df_reg  <- data.frame(
#   y               = Havrankova2025$y,
#   se              = Havrankova2025$se,
#   facing_customer = Havrankova2025$facing_customer,
#   study_id        = Havrankova2025$study_id
# )
# 
# fit_reg <- RoBMA.reg(
#   ~ facing_customer,
#   study_ids      = df_reg$study_id,
#   data           = df_reg,
#   rescale_priors = prior_scale,
#   prior_scale    = "none", transformation = "none",
#   algorithm = "ss", sample = 20000, burnin = 10000, adapt = 10000,
#   thin = 5, parallel = TRUE, autofit = FALSE, seed = 1)
# )

## -----------------------------------------------------------------------------
summary(fit_reg)

## -----------------------------------------------------------------------------
marginal_summary(fit_reg)

## ----eval = FALSE-------------------------------------------------------------
# zfit_reg <- as_zcurve(fit_reg)

## ----fig.width = 6, fig.height = 4--------------------------------------------
par(mar = c(4,4,0,0))
plot(zfit_reg, by.hist = 0.25, plot_extrapolation = FALSE, from = -4, to = 8)

