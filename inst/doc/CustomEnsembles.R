## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
             any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
             !file.exists("../models/CustomEnsembles/Bem_update1.RDS") 
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = !is_check,
  dev      = "png"
)
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## -----------------------------------------------------------------------------
library(RoBMA)

data("Bem2011", package = "RoBMA")
Bem2011

## -----------------------------------------------------------------------------
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study,
             priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
             priors_effect_null        = prior("spike", parameters = list(location = 0)),
             priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
             priors_bias_null          = prior_none(),
             seed = 1)

## -----------------------------------------------------------------------------
summary(fit, type = "models")

## ----fig_mu_prior, dpi = 300, fig.width = 4, fig.height = 3, out.width = "50%", fig.align = "center"----
plot(prior("normal", parameters = list(mean = .15, sd = .10), truncation = list(lower = 0)))

## ----include = FALSE----------------------------------------------------------
# these fits are relatively fast, but we reduce the knitting time considerably
fit <- readRDS(file = "../models/CustomEnsembles/Bem_update1.RDS")

## -----------------------------------------------------------------------------
summary(fit, type = "models")

## ----include = FALSE----------------------------------------------------------
fit <- readRDS(file = "../models/CustomEnsembles/Bem_update2.RDS")

## -----------------------------------------------------------------------------
summary(fit, type = "models")

## -----------------------------------------------------------------------------
summary(fit)

## ----fig_mu_posterior, dpi = 300, fig.width = 4, fig.height = 3.5, out.width = "50%", fig.align = "center"----
plot(fit, parameter = "mu", prior = TRUE)

## ----fig_weightfunction_posterior, dpi = 300, fig.width = 5, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit, parameter = "weightfunction", prior = TRUE)

## ----fig_PETPEESE_posterior, dpi = 300, fig.width = 5, fig.height = 4, out.width = "75%", fig.align = "center"----
plot(fit, parameter = "PET-PEESE", prior = TRUE)

