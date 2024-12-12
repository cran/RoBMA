## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/Tutorial/fit_RoBMA_Lui2015.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev      = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
data("Lui2015", package = "RoBMA")
df <- Lui2015
# preload the fitted model
fit_RoBMA  <- readRDS(file = "../models/Tutorial/fit_RoBMA_Lui2015.RDS")
fit_RoBMA2 <- readRDS(file = "../models/Tutorial/fit_RoBMA_perinull_Lui2015.RDS")

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# data("Lui2015", package = "RoBMA")
# df <- Lui2015
# fit_RoBMA <- RoBMA(r = df$r, n = df$n, seed = 1, model = "PSMA", parallel = TRUE, save = "min")
# 
# fit_RoBMA2 <- RoBMA(r = df$r, n = df$n, seed = 2, parallel = TRUE, save = "min",
#                     priors_effect      = prior("normal", parameters = list(mean = 0.60, sd = 0.20), truncation = list(0, Inf)),
#                     priors_effect_null = prior("normal", parameters = list(mean = 0,    sd = 0.10)))
# 
# saveRDS(fit_RoBMA,  file = "../models/Tutorial/fit_RoBMA_Lui2015.RDS")
# saveRDS(fit_RoBMA2, file = "../models/Tutorial/fit_RoBMA_perinull_Lui2015.RDS")

## ----message = FALSE----------------------------------------------------------
library("metafor")
library("weightr")
library("RoBMA")

## -----------------------------------------------------------------------------
head(df)

## -----------------------------------------------------------------------------
df$r

## -----------------------------------------------------------------------------
dfz <- combine_data(r = df$r, n = df$n, study_names = df$study, transformation = "fishers_z")
head(dfz)

## -----------------------------------------------------------------------------
dfd <- combine_data(r = df$r, n = df$n, study_names = df$study, transformation = "cohens_d")
head(dfd)

## -----------------------------------------------------------------------------
fit_rma <- rma(yi = y, sei = se, data = dfz)
fit_rma

## -----------------------------------------------------------------------------
z2r(fit_rma$b)

## -----------------------------------------------------------------------------
fit_PET <- lm(y ~ se, weights = 1/se^2, data = dfz)
summary(fit_PET)

## -----------------------------------------------------------------------------
z2r(summary(fit_PET)$coefficients["(Intercept)", "Estimate"])

## -----------------------------------------------------------------------------
fit_PEESE <- lm(y ~ I(se^2), weights = 1/se^2, data = dfz)
summary(fit_PEESE)

## -----------------------------------------------------------------------------
fit_4PSM <- weightfunct(effect = dfd$y, v = dfd$se^2, steps = c(0.025, 0.05), table = TRUE)
fit_4PSM

## -----------------------------------------------------------------------------
fit_3PSM <- weightfunct(effect = dfd$y, v = dfd$se^2, steps = c(0.025), table = TRUE)
fit_3PSM

## -----------------------------------------------------------------------------
d2r(fit_3PSM$adj_est[2])

## -----------------------------------------------------------------------------
fit_rma_d <- rma(yi = y, sei = se, data = dfd)

## -----------------------------------------------------------------------------
fit_sel_d <- selmodel(fit_rma_d, type = "stepfun", steps = c(0.025))
fit_sel_d

## -----------------------------------------------------------------------------
summary(fit_RoBMA)

## -----------------------------------------------------------------------------
summary(fit_RoBMA, output_scale = "r")

## -----------------------------------------------------------------------------
interpret(fit_RoBMA, output_scale = "r")

## -----------------------------------------------------------------------------
summary(fit_RoBMA, type = "models", short_name = TRUE)

## -----------------------------------------------------------------------------
summary(fit_RoBMA, type = "diagnostics")

## ----dpi = 300, fig.width = 6, fig.height = 4, out.width = "50%", fig.align = "center"----
par(mar = c(4, 4, 1, 4))
plot(fit_RoBMA, prior = TRUE, output_scale = "r", )

## -----------------------------------------------------------------------------
summary(fit_RoBMA2, type = "models")

