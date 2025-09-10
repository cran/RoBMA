## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/zcurve/zcurve_RE_Hoppen2025.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev      = "png",
  fig.width = 7,
  fig.height = 5)
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## ----include = FALSE, eval = FALSE--------------------------------------------
# # R package version updating
# library(RoBMA)
# 
# # Ease of retrieval - Hoppen2025
# data("Weingarten2018", package = "RoBMA")
# Weingarten2018 <- Weingarten2018[Weingarten2018$standard_paradigm & Weingarten2018$proximal_dataset, ]
# fit_RE_Weingarten2018 <- NoBMA(r = Weingarten2018$r_xy, n = round(Weingarten2018$N), study_ids = Weingarten2018$paper_id,
#                            priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                            algorithm = "ss", sample = 10000, burnin = 10000, adapt = 10000,
#                            chains = 5, parallel = TRUE, seed = 1, save = "min")
# fit_RoBMA_Weingarten2018 <- RoBMA(r = Weingarten2018$r_xy, n = round(Weingarten2018$N), study_ids = Weingarten2018$paper_id,
#                               algorithm = "ss", sample = 10000, burnin = 10000, adapt = 10000,
#                               chains = 5, parallel = TRUE, seed = 1, save = "min")
# 
# zcurve_RE_Weingarten2018    <- as_zcurve(fit_RE_Weingarten2018)
# zcurve_RoBMA_Weingarten2018 <- as_zcurve(fit_RoBMA_Weingarten2018)
# 
# saveRDS(zcurve_RE_Weingarten2018,    file = "../models/zcurve/zcurve_RE_Weingarten2018.RDS", compress = "xz")
# saveRDS(zcurve_RoBMA_Weingarten2018, file = "../models/zcurve/zcurve_RoBMA_Weingarten2018.RDS", compress = "xz")
# 
# 
# # Social comparison example - Hoppen2025
# data("Hoppen2025", package = "RoBMA")
# fit_RE_Hoppen2025 <- NoBMA(d = Hoppen2025$d, se = sqrt(Hoppen2025$v),
#                            priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                            algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                            chains = 5, parallel = TRUE, seed = 1, save = "min")
# fit_RoBMA_Hoppen2025 <- RoBMA(d = Hoppen2025$d, se = sqrt(Hoppen2025$v),
#                               algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                               chains = 5, parallel = TRUE, seed = 1, save = "min")
# 
# zcurve_RE_Hoppen2025    <- as_zcurve(fit_RE_Hoppen2025)
# zcurve_RoBMA_Hoppen2025 <- as_zcurve(fit_RoBMA_Hoppen2025)
# 
# saveRDS(zcurve_RE_Hoppen2025,    file = "../models/zcurve/zcurve_RE_Hoppen2025.RDS", compress = "xz")
# saveRDS(zcurve_RoBMA_Hoppen2025, file = "../models/zcurve/zcurve_RoBMA_Hoppen2025.RDS", compress = "xz")
# 
# # ChatGPT example - Wang2025
# data("Wang2025", package = "RoBMA")
# Wang2025 <- Wang2025[Wang2025$Learning_effect == "Learning performance", ]
# fit_RE_Wang2025 <- NoBMA(d = Wang2025$g, se = Wang2025$se,
#                          priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                          algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                          chains = 5, parallel = TRUE, seed = 1, save = "min")
# fit_RoBMA_Wang2025 <- RoBMA(d = Wang2025$g, se = Wang2025$se,
#                             algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                             chains = 5, parallel = TRUE, seed = 1, save = "min")
# 
# zcurve_RE_Wang2025    <- as_zcurve(fit_RE_Wang2025)
# zcurve_RoBMA_Wang2025 <- as_zcurve(fit_RoBMA_Wang2025)
# 
# saveRDS(zcurve_RE_Wang2025,    file = "../models/zcurve/zcurve_RE_Wang2025.RDS", compress = "xz")
# saveRDS(zcurve_RoBMA_Wang2025, file = "../models/zcurve/zcurve_RoBMA_Wang2025.RDS", compress = "xz")
# 
# # Many Labs 2 example - ManyLabs16
# data("ManyLabs16", package = "RoBMA")
# fit_RE_ManyLabs16 <- NoBMA(d = ManyLabs16$y, se = ManyLabs16$se,
#                            priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                            algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                            chains = 5, parallel = TRUE, seed = 1, save = "min")
# fit_RoBMA_ManyLabs16 <- RoBMA(d = ManyLabs16$y, se = ManyLabs16$se,
#                               algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                               chains = 5, parallel = TRUE, seed = 1, save = "min")
# 
# zcurve_RE_ManyLabs16    <- as_zcurve(fit_RE_ManyLabs16)
# zcurve_RoBMA_ManyLabs16 <- as_zcurve(fit_RoBMA_ManyLabs16)
# 
# saveRDS(zcurve_RE_ManyLabs16,    file = "../models/zcurve/zcurve_RE_ManyLabs16.RDS", compress = "xz")
# saveRDS(zcurve_RoBMA_ManyLabs16, file = "../models/zcurve/zcurve_RoBMA_ManyLabs16.RDS", compress = "xz")

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

# preload the fitted models
zcurve_RE_Weingarten2018    <- readRDS(file = "../models/zcurve/zcurve_RE_Weingarten2018.RDS")
zcurve_RoBMA_Weingarten2018 <- readRDS(file = "../models/zcurve/zcurve_RoBMA_Weingarten2018.RDS")
zcurve_RE_Hoppen2025    <- readRDS(file = "../models/zcurve/zcurve_RE_Hoppen2025.RDS")
zcurve_RoBMA_Hoppen2025 <- readRDS(file = "../models/zcurve/zcurve_RoBMA_Hoppen2025.RDS")
zcurve_RE_Wang2025      <- readRDS(file = "../models/zcurve/zcurve_RE_Wang2025.RDS")
zcurve_RoBMA_Wang2025   <- readRDS(file = "../models/zcurve/zcurve_RoBMA_Wang2025.RDS")
zcurve_RE_ManyLabs16    <- readRDS(file = "../models/zcurve/zcurve_RE_ManyLabs16.RDS")
zcurve_RoBMA_ManyLabs16 <- readRDS(file = "../models/zcurve/zcurve_RoBMA_ManyLabs16.RDS")

fit_RE_Weingarten2018    <- zcurve_RE_Weingarten2018
fit_RoBMA_Weingarten2018 <- zcurve_RoBMA_Weingarten2018
fit_RE_Hoppen2025    <- zcurve_RE_Hoppen2025
fit_RoBMA_Hoppen2025 <- zcurve_RoBMA_Hoppen2025
fit_RE_Wang2025      <- zcurve_RE_Wang2025
fit_RoBMA_Wang2025   <- zcurve_RoBMA_Wang2025
fit_RE_ManyLabs16    <- zcurve_RE_ManyLabs16
fit_RoBMA_ManyLabs16 <- zcurve_RoBMA_ManyLabs16

class(fit_RE_Weingarten2018)    <- class(fit_RE_Weingarten2018)[!class(fit_RE_Weingarten2018)       %in% "zcurve_RoBMA"]
class(fit_RoBMA_Weingarten2018) <- class(fit_RoBMA_Weingarten2018)[!class(fit_RoBMA_Weingarten2018) %in% "zcurve_RoBMA"]
class(fit_RE_Hoppen2025)    <- class(fit_RE_Hoppen2025)[!class(fit_RE_Hoppen2025)       %in% "zcurve_RoBMA"]
class(fit_RoBMA_Hoppen2025) <- class(fit_RoBMA_Hoppen2025)[!class(fit_RoBMA_Hoppen2025) %in% "zcurve_RoBMA"]
class(fit_RE_Wang2025)      <- class(fit_RE_Wang2025)[!class(fit_RE_Wang2025)           %in% "zcurve_RoBMA"]
class(fit_RoBMA_Wang2025)   <- class(fit_RoBMA_Wang2025)[!class(fit_RoBMA_Wang2025)     %in% "zcurve_RoBMA"]
class(fit_RE_ManyLabs16)    <- class(fit_RE_ManyLabs16)[!class(fit_RE_ManyLabs16)       %in% "zcurve_RoBMA"]
class(fit_RoBMA_ManyLabs16) <- class(fit_RoBMA_ManyLabs16)[!class(fit_RoBMA_ManyLabs16) %in% "zcurve_RoBMA"]

## ----message = FALSE----------------------------------------------------------
library("RoBMA")

## -----------------------------------------------------------------------------
# Load the ease-of-retrieval dataset
data("Weingarten2018", package = "RoBMA")
# Filter to standard paradigm and proximal dataset as in the original analysis
Weingarten2018 <- Weingarten2018[Weingarten2018$standard_paradigm & Weingarten2018$proximal_dataset, ]
head(Weingarten2018)

## ----eval = FALSE-------------------------------------------------------------
# # Fit random-effects model (unadjusted for publication bias)
# fit_RE_Weingarten2018 <- NoBMA(r = Weingarten2018$r_xy, n = round(Weingarten2018$N),
#                                study_ids = Weingarten2018$paper_id,
#                                priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                                algorithm = "ss", sample = 10000, burnin = 10000, adapt = 10000,
#                                chains = 5, parallel = TRUE, seed = 1)
# 
# # Fit RoBMA model (adjusted for publication bias)
# fit_RoBMA_Weingarten2018 <- RoBMA(r = Weingarten2018$r_xy, n = round(Weingarten2018$N),
#                                   study_ids = Weingarten2018$paper_id,
#                                   algorithm = "ss", sample = 10000, burnin = 10000, adapt = 10000,
#                                   chains = 5, parallel = TRUE, seed = 1)

## -----------------------------------------------------------------------------
# Random-effects model results
summary(fit_RE_Weingarten2018, output_scale = "r")

# RoBMA model results  
summary(fit_RoBMA_Weingarten2018, output_scale = "r")

## ----eval = FALSE-------------------------------------------------------------
# # Create z-curve objects
# zcurve_RE_Weingarten2018    <- as_zcurve(fit_RE_Weingarten2018)
# zcurve_RoBMA_Weingarten2018 <- as_zcurve(fit_RoBMA_Weingarten2018)

## ----fig.cap="Ease-of-Retrieval Effect: Model Fit Assessment"-----------------
# Create histogram of observed z-statistics
hist(zcurve_RoBMA_Weingarten2018, from = -3, to = 6, by = 0.25)

# Add model-implied distributions
lines(zcurve_RE_Weingarten2018, from = -3, to = 6, col = "black", lty = 2, lwd = 2)
lines(zcurve_RoBMA_Weingarten2018, from = -3, to = 6, col = "blue", lty = 2, lwd = 2)

# Add legend
legend("topright", 
       legend = c("Random-Effects", "RoBMA"), 
       col = c("black", "blue"), 
       lty = 2, lwd = 2)

## ----fig.cap="Ease-of-Retrieval Effect: Extrapolation Analysis"---------------
# Plot extrapolation to pre-publication bias
plot(zcurve_RoBMA_Weingarten2018, from = -3, to = 6, by.hist = 0.25)

## -----------------------------------------------------------------------------
# Extract z-curve summary metrics
summary(zcurve_RoBMA_Weingarten2018)

## -----------------------------------------------------------------------------
# Load the social comparison dataset
data("Hoppen2025", package = "RoBMA")
head(Hoppen2025)

## ----eval = FALSE-------------------------------------------------------------
# # Fit random-effects model (unadjusted for publication bias)
# fit_RE_Hoppen2025 <- NoBMA(d = Hoppen2025$d, se = sqrt(Hoppen2025$v),
#                            priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                            algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                            chains = 5, parallel = TRUE, seed = 1)
# 
# # Fit RoBMA model (adjusted for publication bias)
# fit_RoBMA_Hoppen2025 <- RoBMA(d = Hoppen2025$d, se = sqrt(Hoppen2025$v),
#                               algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                               chains = 5, parallel = TRUE, seed = 1)

## -----------------------------------------------------------------------------
# Random-effects model results
summary(fit_RE_Hoppen2025)

## -----------------------------------------------------------------------------
# RoBMA model results  
summary(fit_RoBMA_Hoppen2025)

## ----eval = FALSE-------------------------------------------------------------
# # Create z-curve objects
# zcurve_RE_Hoppen2025    <- as_zcurve(fit_RE_Hoppen2025)
# zcurve_RoBMA_Hoppen2025 <- as_zcurve(fit_RoBMA_Hoppen2025)

## ----fig.cap="Social Comparison: Model Fit Assessment"------------------------
# Create histogram of observed z-statistics
hist(zcurve_RoBMA_Hoppen2025)

# Add model-implied distributions
lines(zcurve_RE_Hoppen2025, col = "black", lty = 2, lwd = 2)
lines(zcurve_RoBMA_Hoppen2025, col = "blue", lty = 2, lwd = 2)

# Add legend
legend("topright", 
       legend = c("Random-Effects", "RoBMA"), 
       col = c("black", "blue"), 
       lty = 2, lwd = 2)

## ----fig.cap="Social Comparison: Extrapolation Analysis"----------------------
# Plot extrapolation to pre-publication bias
plot(zcurve_RoBMA_Hoppen2025)

## -----------------------------------------------------------------------------
# Extract z-curve summary metrics
summary(zcurve_RoBMA_Hoppen2025)

## -----------------------------------------------------------------------------
# Load the ChatGPT dataset
data("Wang2025", package = "RoBMA")
# Select learning performance studies
Wang2025 <- Wang2025[Wang2025$Learning_effect == "Learning performance", ]
head(Wang2025)

## ----eval = FALSE-------------------------------------------------------------
# # Fit models
# fit_RE_Wang2025 <- NoBMA(d = Wang2025$g, se = Wang2025$se,
#                          priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                          algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                          chains = 5, parallel = TRUE, seed = 1)
# fit_RoBMA_Wang2025 <- RoBMA(d = Wang2025$g, se = Wang2025$se,
#                             algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                             chains = 5, parallel = TRUE, seed = 1)
# 
# # Create z-curve objects
# zcurve_RE_Wang2025    <- as_zcurve(fit_RE_Wang2025)
# zcurve_RoBMA_Wang2025 <- as_zcurve(fit_RoBMA_Wang2025)

## ----fig.cap="ChatGPT: Model Fit Assessment"----------------------------------
# Assess model fit
hist(zcurve_RoBMA_Wang2025, from = -2, to = 8)
lines(zcurve_RE_Wang2025, col = "black", lty = 2, lwd = 2, from = -2, to = 8)
lines(zcurve_RoBMA_Wang2025, col = "blue", lty = 2, lwd = 2, from = -2, to = 8)
legend("topright", 
       legend = c("Random-Effects", "RoBMA"), 
       col = c("black", "blue"), 
       lty = 2, lwd = 2)

## ----fig.cap="ChatGPT: Extrapolation Analysis"--------------------------------
# Examine extrapolation
plot(zcurve_RoBMA_Wang2025, from = -2, to = 8)

## -----------------------------------------------------------------------------
summary(fit_RE_Wang2025)
summary(fit_RoBMA_Wang2025)

## -----------------------------------------------------------------------------
# Extract z-curve summary metrics for ChatGPT data
summary(zcurve_RoBMA_Wang2025)

## -----------------------------------------------------------------------------
# Load the Many Labs 2 framing effect data
data("ManyLabs16", package = "RoBMA")
head(ManyLabs16)

## ----eval = FALSE-------------------------------------------------------------
# # Fit models and create z-curve objects
# fit_RE_ManyLabs16 <- NoBMA(d = ManyLabs16$y, se = ManyLabs16$se,
#                            priors_effect_null = NULL, priors_heterogeneity_null = NULL,
#                            algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                            chains = 5, parallel = TRUE, seed = 1)
# fit_RoBMA_ManyLabs16 <- RoBMA(d = ManyLabs16$y, se = ManyLabs16$se,
#                               algorithm = "ss", sample = 10000, burnin = 5000, adapt = 5000,
#                               chains = 5, parallel = TRUE, seed = 1)
# 
# zcurve_RE_ManyLabs16    <- as_zcurve(fit_RE_ManyLabs16)
# zcurve_RoBMA_ManyLabs16 <- as_zcurve(fit_RoBMA_ManyLabs16)

## ----fig.cap="Framing Effects: Model Fit Assessment"--------------------------
# Assess model fit - should show good agreement
hist(zcurve_RoBMA_ManyLabs16)
lines(zcurve_RE_ManyLabs16, col = "black", lty = 2, lwd = 2)
lines(zcurve_RoBMA_ManyLabs16, col = "blue", lty = 2, lwd = 2)
legend("topleft", 
       legend = c("Random-Effects", "RoBMA"), 
       col = c("black", "blue"), 
       lty = 2, lwd = 2)

## ----fig.cap="Framing Effects: Extrapolation Analysis"------------------------
# Extrapolation should show minimal change
plot(zcurve_RoBMA_ManyLabs16)

## ----eval = FALSE-------------------------------------------------------------
# summary(fit_RE_ManyLabs16)
# summary(fit_RoBMA_ManyLabs16)

## -----------------------------------------------------------------------------
# Extract z-curve summary metrics for Many Labs 2 data
summary(zcurve_RoBMA_ManyLabs16)

