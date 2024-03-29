% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RoBMA.R
\name{update.RoBMA}
\alias{update.RoBMA}
\title{Updates a fitted RoBMA object}
\usage{
\method{update}{RoBMA}(
  object,
  refit_failed = TRUE,
  extend_all = FALSE,
  prior_effect = NULL,
  prior_heterogeneity = NULL,
  prior_bias = NULL,
  prior_hierarchical = NULL,
  prior_weights = NULL,
  prior_effect_null = NULL,
  prior_heterogeneity_null = NULL,
  prior_bias_null = NULL,
  prior_hierarchical_null = NULL,
  study_names = NULL,
  chains = NULL,
  adapt = NULL,
  burnin = NULL,
  sample = NULL,
  thin = NULL,
  autofit = NULL,
  parallel = NULL,
  autofit_control = NULL,
  convergence_checks = NULL,
  save = "all",
  seed = NULL,
  silent = TRUE,
  ...
)
}
\arguments{
\item{object}{a fitted RoBMA object}

\item{refit_failed}{whether failed models should be refitted. Relevant only
if new priors or \code{prior_weights} are not supplied. Defaults to \code{TRUE}.}

\item{extend_all}{extend sampling in all fitted models based on \code{"sample_extend"}
argument in \code{\link[=set_autofit_control]{set_autofit_control()}} function. Defaults to \code{FALSE}.}

\item{prior_effect}{prior distribution for the effect size (\code{mu})
parameter that will be treated as belonging to the alternative hypothesis.
Defaults to \code{NULL}.}

\item{prior_heterogeneity}{prior distribution for the heterogeneity \code{tau}
parameter that will be treated as belonging to the alternative hypothesis.
Defaults to \code{NULL}.}

\item{prior_bias}{prior distribution for the publication bias adjustment
component that will be treated as belonging to the alternative hypothesis.
Defaults to \code{NULL}.}

\item{prior_hierarchical}{prior distribution for the correlation of random effects
(\code{rho}) parameter that will be treated as belonging to the alternative hypothesis. This setting allows
users to fit a hierarchical (three-level) meta-analysis when \code{study_ids} are supplied.
Note that this is an experimental feature and see News for more details. Defaults to a beta distribution
\code{prior(distribution = "beta", parameters = list(alpha = 1, beta = 1))}.}

\item{prior_weights}{either a single value specifying prior model weight
of a newly specified model using priors argument, or a vector of the
same length as already fitted models to update their prior weights.}

\item{prior_effect_null}{prior distribution for the effect size (\code{mu})
parameter that will be treated as belonging to the null hypothesis.
Defaults to \code{NULL}.}

\item{prior_heterogeneity_null}{prior distribution for the heterogeneity \code{tau}
parameter that will be treated as belonging to the null hypothesis.
Defaults to \code{NULL}.}

\item{prior_bias_null}{prior distribution for the publication bias adjustment
component that will be treated as belonging to the null hypothesis.
Defaults to \code{NULL}.}

\item{prior_hierarchical_null}{prior distribution for the correlation of random effects
(\code{rho}) parameter that will be treated as belonging to the null hypothesis. Defaults to \code{NULL}.}

\item{study_names}{an optional argument with the names of the studies}

\item{chains}{a number of chains of the MCMC algorithm.}

\item{adapt}{a number of adaptation iterations of the MCMC algorithm.
Defaults to \code{500}.}

\item{burnin}{a number of burnin iterations of the MCMC algorithm.
Defaults to \code{2000}.}

\item{sample}{a number of sampling iterations of the MCMC algorithm.
Defaults to \code{5000}.}

\item{thin}{a thinning of the chains of the MCMC algorithm. Defaults to
\code{1}.}

\item{autofit}{whether the model should be fitted until the convergence
criteria (specified in \code{autofit_control}) are satisfied. Defaults to
\code{TRUE}.}

\item{parallel}{whether the individual models should be fitted in parallel.
Defaults to \code{FALSE}. The implementation is not completely stable
and might cause a connection error.}

\item{autofit_control}{allows to pass autofit control settings with the
\code{\link[=set_autofit_control]{set_autofit_control()}} function. See \code{?set_autofit_control} for
options and default settings.}

\item{convergence_checks}{automatic convergence checks to assess the fitted
models, passed with \code{\link[=set_convergence_checks]{set_convergence_checks()}} function. See
\code{?set_convergence_checks} for options and default settings.}

\item{save}{whether all models posterior distributions should be kept
after obtaining a model-averaged result. Defaults to \code{"all"} which
does not remove anything. Set to \code{"min"} to significantly reduce
the size of final object, however, some model diagnostics and further
manipulation with the object will not be possible.}

\item{seed}{a seed to be set before model fitting, marginal likelihood
computation, and posterior mixing for reproducibility of results. Defaults
to \code{NULL} - no seed is set.}

\item{silent}{whether all print messages regarding the fitting process
should be suppressed. Defaults to \code{TRUE}. Note that \code{parallel = TRUE}
also suppresses all messages.}

\item{...}{additional arguments.}
}
\value{
\code{RoBMA} returns an object of class 'RoBMA'.
}
\description{
\code{update.RoBMA} can be used to
\enumerate{
\item{add an additional model to an existing \code{"RoBMA"} object by
specifying either a null or alternative prior for each parameter
and the prior odds of the model (\code{prior_weights}), see the
\code{vignette("CustomEnsembles")} vignette,}
\item{change the prior odds of fitted models by specifying a vector
\code{prior_weights} of the same length as the fitted models,}
\item{refitting models that failed to converge with updated settings
of control parameters,}
\item{or changing the convergence criteria and recalculating the ensemble
results by specifying new \code{control} argument and setting
\code{refit_failed == FALSE}.}
}
}
\details{
See \code{\link[=RoBMA]{RoBMA()}} for more details.
}
\examples{
\dontrun{
# using the example data from Bem 2011 and fitting the default (RoBMA-PSMA) model
fit <- RoBMA(d = Bem2011$d, se = Bem2011$se, study_names = Bem2011$study)

# the update function allows us to change the prior model weights of each model
fit1 <- update(fit, prior_weights = c(0, rep(1, 35)))

# add an additional model with different priors specification
# (see '?prior' for more information)
fit2 <- update(fit,
               priors_effect_null = prior("point", parameters = list(location = 0)),
               priors_heterogeneity = prior("normal",
                                  parameters = list(mean = 0, sd = 1),
                                  truncation = list(lower = 0, upper = Inf)),
               priors_bias = prior_weightfunction("one-sided",
                                    parameters = list(cuts = c(.05, .10, .20),
                                                      alpha = c(1, 1, 1, 1))))

# update the models with an increased number of sample iterations
fit3 <- update(fit, autofit_control = set_autofit_control(sample_extend = 1000), extend_all = TRUE)
}


}
\seealso{
\code{\link[=RoBMA]{RoBMA()}}, \code{\link[=summary.RoBMA]{summary.RoBMA()}}, \code{\link[=prior]{prior()}}, \code{\link[=check_setup]{check_setup()}}
}
