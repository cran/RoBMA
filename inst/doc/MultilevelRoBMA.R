## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/MultilevelRoBMA/fit_Johnides2025.RDS")
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
fit        <- readRDS(file = "../models/MultilevelRoBMA/fit_Johnides2025.RDS")
fit_simple <- readRDS(file = "../models/MultilevelRoBMA/fit_Johnides2025_single.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# data("Johnides2025", package = "RoBMA")
# 
# fit <- RoBMA(
#   d         = Johnides2025$d,
#   se        = Johnides2025$se,
#   study_ids = Johnides2025$study,
#   algorithm = "ss",
#   adapt     = 5000,
#   burnin    = 5000,
#   sample    = 10000,
#   parallel  = TRUE,
#   seed      = 1,
#   autofit   = FALSE
# )
# saveRDS(fit, file = "../models/MultilevelRoBMA/fit_Johnides2025.RDS", compress = "xz")
# 
# fit_simple <- RoBMA(
#   d         = Johnides2025$d,
#   se        = Johnides2025$se,
#   algorithm = "ss",
#   adapt     = 5000,
#   burnin    = 5000,
#   sample    = 10000,
#   parallel  = TRUE,
#   seed      = 1,
#   autofit   = FALSE
# )
# saveRDS(fit_simple, file = "../models/MultilevelRoBMA/fit_Johnides2025_single.RDS", compress = "xz")

## -----------------------------------------------------------------------------
library(RoBMA)
data("Johnides2025", package = "RoBMA")

## -----------------------------------------------------------------------------
head(Johnides2025)

## ----eval = FALSE-------------------------------------------------------------
# fit <- RoBMA(
#   d         = Johnides2025$d,
#   se        = Johnides2025$se,
#   study_ids = Johnides2025$study,
#   algorithm = "ss",
#   adapt     = 5000,
#   burnin    = 5000,
#   sample    = 10000,
#   parallel  = TRUE,
#   seed      = 1,
#   autofit   = FALSE
# )

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
summary_heterogeneity(fit)

## -----------------------------------------------------------------------------
summary(fit, type = "models")

## ----fig.width = 6, fig.height = 4--------------------------------------------
plot(fit, parameter = "weightfunction", rescale_x = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# fit_simple <- RoBMA(
#   d         = Johnides2025$d,
#   se        = Johnides2025$se,
#   algorithm = "ss",
#   adapt     = 5000,
#   burnin    = 5000,
#   sample    = 10000,
#   parallel  = TRUE,
#   seed      = 1,
#   autofit   = FALSE
# )

## -----------------------------------------------------------------------------
summary(fit_simple)

