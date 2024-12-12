## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/MedicineBiBMA/fit.RDS")
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
fit <- readRDS(file = "../models/MedicineBiBMA/fit.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# 
# data("Poulsen2006", package = "RoBMA")
# 
# # p. 73 in https://www.cochranelibrary.com/cdsr/doi/10.1002/14651858.CD007094.pub5/epdf/full
# events_experimental        <- c(5, 2)
# events_control             <- c(0, 0)
# observations_experimental  <- c(35, 40)
# observations_control       <- c(39, 40)
# study_names <- c("Paul 2007", "Shadkam 2010")
# 
# # domain specific prior distributions: Acute Respiratory Infections
# fit <- BiBMA(
#   x1          = events_experimental,
#   x2          = events_control,
#   n1          = observations_experimental,
#   n2          = observations_control,
#   study_names = study_names,
#   priors_effect        = prior_informed("Acute Respiratory Infections", type = "logOR", parameter = "effect"),
#   priors_heterogeneity = prior_informed("Acute Respiratory Infections", type = "logOR", parameter = "heterogeneity"),
#   seed = 1
# )
# 
# saveRDS(fit,  file = "../models/MedicineBiBMA/fit.RDS")

## -----------------------------------------------------------------------------
library(RoBMA)

events_experimental        <- c(5, 2)
events_control             <- c(0, 0)
observations_experimental  <- c(35, 40)
observations_control       <- c(39, 40)
study_names <- c("Paul 2007", "Shadkam 2010")

## -----------------------------------------------------------------------------
summary(fit, conditional = TRUE, output_scale = "OR")

## ----fig_mu_BMA, dpi = 300, fig.width = 6, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit, parameter = "mu", prior = TRUE, conditional = TRUE)

