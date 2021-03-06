# the main functions
.fit_RoBMA_model       <- function(object, i){

  model              <- object[["models"]][[i]]
  priors             <- model[["priors"]]
  fit_control        <- object[["fit_control"]]
  autofit_control    <- object[["autofit_control"]]
  convergence_checks <- object[["convergence_checks"]]
  add_info           <- object[["add_info"]]

  errors   <- NULL
  warnings <- NULL

  # don't sample the complete null model
  if(!.is_model_constant(priors)){

    # generate the model syntax
    model_syntax <- .generate_model_syntax(priors, add_info[["effect_direction"]], add_info[["prior_scale"]], add_info[["effect_measure"]])

    # remove unnecessary objects from data to mitigate warnings
    fit_data     <- .fit_data(object[["data"]], priors, add_info[["effect_direction"]], add_info[["prior_scale"]])

    # fit the model
    fit <- BayesTools::JAGS_fit(
      model_syntax       = model_syntax,
      data               = fit_data,
      prior_list         = priors,
      chains             = fit_control[["chains"]],
      adapt              = fit_control[["adapt"]],
      burnin             = fit_control[["burnin"]],
      sample             = fit_control[["sample"]],
      thin               = fit_control[["thin"]],
      autofit            = fit_control[["autofit"]],
      autofit_control    = autofit_control,
      parallel           = fit_control[["parallel"]],
      cores              = fit_control[["cores"]],
      silent             = fit_control[["silent"]],
      seed               = fit_control[["seed"]],
      required_packages  = "RoBMA"
    )

    # assess the model fit and deal with errors
    if(inherits(fit, "error")){

      fit            <- list()
      converged      <- FALSE
      has_posterior  <- FALSE
      errors         <- c(errors, fit$message)
      # deal with failed models
      marglik        <- list()
      marglik$logml  <- NA
      class(marglik) <- "bridge"

    }else{

      has_posterior <- TRUE
      check_fit     <- BayesTools::JAGS_check_convergence(
        fit          = fit,
        prior_list   = priors,
        max_Rhat     = convergence_checks[["max_Rhat"]],
        min_ESS      = convergence_checks[["min_ESS"]],
        max_error    = convergence_checks[["max_error"]],
        max_SD_error = convergence_checks[["max_SD_error"]]
      )
      warnings    <- c(warnings, attr(fit, "warnings"), attr(check_fit, "errors"))
      if(convergence_checks[["remove_failed"]] && !check_fit){
        converged <- FALSE
      }else{
        converged <- TRUE
      }

    }

    # compute marginal likelihood
    if(length(fit) != 0){

      marglik <- BayesTools::JAGS_bridgesampling(
        fit              = fit,
        data             = fit_data,
        prior_list       = priors,
        log_posterior    = .marglik_function,
        maxiter          = 50000,
        silent           = fit_control[["silent"]],
        priors           = priors,
        effect_direction = add_info[["effect_direction"]],
        prior_scale      = add_info[["prior_scale"]],
        effect_measure   = add_info[["effect_measure"]]
      )

      # deal with failed marginal likelihoods
      if(inherits(marglik, "error")){

        errors         <- c(errors, marglik$message)
        converged      <- FALSE
        marglik        <- list()
        marglik$logml  <- NA
        class(marglik) <- "bridge"

      }else{

        # forward warnings if present
        warnings <- c(warnings, attr(marglik, "warnings"))

      }
    }


  }else{

    fit_data       <- .fit_data(object[["data"]], priors, add_info[["effect_direction"]], add_info[["prior_scale"]])
    converged      <- TRUE
    has_posterior  <- FALSE
    fit            <- list()
    class(fit)     <- "null_model"
    marglik        <- list()
    marglik$logml  <- sum(stats::dnorm(fit_data[["y"]], priors$mu$parameters[["location"]], fit_data[["se"]], log = TRUE))
    class(marglik) <- "bridge"

  }

  # add model summaries
  if(has_posterior){
    fit_summary   <- BayesTools::runjags_estimates_table(fit = fit, prior_list = priors, warnings = warnings)
    if(add_info[["prior_scale"]] != "y"){
      fit_summaries <- .runjags_summary_list(fit, priors, add_info[["prior_scale"]], warnings)
    }else{
      fit_summaries <- NULL
    }
  }else{
    fit_summary    <- BayesTools::runjags_estimates_empty_table()
    if(add_info[["prior_scale"]] != "y"){
      fit_summaries  <- list(
        "d"     = BayesTools::runjags_estimates_empty_table(),
        "r"     = BayesTools::runjags_estimates_empty_table(),
        "logOR" = BayesTools::runjags_estimates_empty_table(),
        "z"     = BayesTools::runjags_estimates_empty_table()
      )
    }else{
      fit_summaries <- NULL
    }
  }

  model <- c(
    model,
    fit           = list(fit),
    fit_summary   = list(fit_summary),
    fit_summaries = list(fit_summaries),
    marglik       = list(marglik),
    errors        = list(errors),
    warnings      = list(warnings),
    converged     = converged,
    has_posterior = has_posterior,
    output_scale  = add_info[["prior_scale"]],
    prior_scale   = add_info[["prior_scale"]]
  )

  return(model)
}

# tools
.fit_data              <- function(data, priors, effect_direction, prior_scale){

  # unlist the data.frame
  original_measure <- attr(data, "original_measure")
  effect_measure   <- attr(data, "effect_measure")

  fit_data <- list()
  # change the effect size direction (important for one-sided selection and PET/PEESE)
  if(effect_direction == "negative"){
    fit_data$y <- - data[,"y"]
  }else{
    fit_data$y <- data[,"y"]
  }
  fit_data$se <- data[,"se"]
  fit_data$K  <- length(data[["y"]])

  # add critical y-values
  if(!is.null(priors[["omega"]])){
    fit_data$crit_y  <- .get_cutoffs(fit_data[["y"]], fit_data[["se"]], priors[["omega"]], original_measure, effect_measure)
  }

  return(fit_data)
}
.generate_model_syntax <- function(priors, effect_direction, priors_scale, effect_measure){

  model_syntax <- "model{\n"

  ### prior transformations
  # the precise transformation for heterogeneity is not used due the inability to re-scale large variances
  # instead, approximate linear scaling is employed in the same way as in metaBMA package
  model_syntax <- paste0(model_syntax, .JAGS_scale(priors_scale, effect_measure, "mu",  "mu_transformed"))
  model_syntax <- paste0(model_syntax, .JAGS_scale(priors_scale, effect_measure, "tau", "tau_transformed"))
  if(!is.null(priors[["PET"]])){
    model_syntax <- paste0(model_syntax, paste0("PET_transformed = PET\n"))
  }else if(!is.null(priors[["PEESE"]])){
    # don't forget that the transformation is inverse for PEESE
    model_syntax <- paste0(model_syntax, .JAGS_scale(effect_measure, priors_scale, "PEESE", "PEESE_transformed"))
  }


  ### model
  model_syntax <- paste0(model_syntax, "for(i in 1:K){\n")

  # marginalized random effects and the effect size
  prec     <- "1 / ( pow(se[i],2) + pow(tau_transformed,2) )"
  reg_std  <- "pow( pow(se[i],2) + pow(tau_transformed,2), 1/2)"

  eff <- ifelse(effect_direction == "negative", "-1 * mu_transformed", "mu_transformed")
  # add PET/PEESE
  if(!is.null(priors[["PET"]])){
    eff <- paste0("(", eff, " + PET_transformed * se[i])")
  }else if(!is.null(priors[["PEESE"]])){
    eff <- paste0("(", eff, " + PEESE_transformed * pow(se[i],2))")
  }

  # the observed data
  if(is.null(priors[["omega"]])){
    model_syntax <- paste0(model_syntax, "  y[i] ~ dnorm(",     eff, ",", prec, " )\n")
  }else if(grepl("one.sided", priors[["omega"]]$distribution)){
    model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_1s(", eff, ",", prec, ", crit_y[i,], omega) \n")
  }else if(grepl("two.sided", priors[["omega"]]$distribution)){
    model_syntax <- paste0(model_syntax, "  y[i] ~ dwnorm_2s(", eff, ",", prec, ", crit_y[i,], omega) \n")
  }


  model_syntax <- paste0(model_syntax, "}\n")
  model_syntax <- paste0(model_syntax, "}")

  return(model_syntax)
}
.marglik_function      <- function(parameters, data, priors, effect_direction, prior_scale, effect_measure){

  # extract parameters
  mu  <- parameters[["mu"]]
  tau <- parameters[["tau"]]
  if(!is.null(priors[["PET"]])){
    PET    <- parameters[["PET"]]
  }else if(!is.null(priors[["PEESE"]])){
    PEESE  <- parameters[["PEESE"]]
  }else if(!is.null(priors[["omega"]])){
    omega  <- parameters[["omega"]]
  }

  ### re-scale parameters
  if(prior_scale != effect_measure){
    mu_transformed  <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(mu)
    )
    tau_transformed <- do.call(
      .get_scale(prior_scale, effect_measure),
      args = list(tau)
    )
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      # don't forget that the transformation is inverse for PEESE
      PEESE_transformed <- do.call(
        .get_scale(effect_measure, prior_scale),
        args = list(PEESE)
      )
    }
  }else{
    mu_transformed  <- mu
    tau_transformed <- tau
    if(!is.null(priors[["PET"]])){
      PET_transformed   <- PET
    }else if(!is.null(priors[["PEESE"]])){
      PEESE_transformed <- PEESE
    }
  }

  ### model
  # marginalized random effects and the effect size
  pop_sd  <- sqrt(data[["se"]]^2 + tau_transformed^2)

  # add PET/PEESE
  eff <- ifelse(effect_direction == "negative", -1, 1) * mu_transformed
  if(!is.null(priors[["PET"]])){
    eff <- eff + PET_transformed   * data[["se"]]
  }else if(!is.null(priors[["PEESE"]])){
    eff <- eff + PEESE_transformed * data[["se"]]^2
  }

  ### compute the marginal log_likelihood
  log_lik <- 0

  # the individual studies
  if(is.null(priors[["omega"]])){
    log_lik <- log_lik + sum(stats::dnorm(data[["y"]], mean = eff, sd = pop_sd, log = TRUE))
  }else if(priors[["omega"]]$distribution == "one.sided"){
    log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = data[["crit_y"]], type = "one.sided", log = TRUE))
  }else if(priors[["omega"]]$distribution == "two.sided"){
    log_lik <- log_lik + sum(.dwnorm_fast(data[["y"]], mean = eff, sd = pop_sd, omega = omega, crit_x = data[["crit_y"]], type = "two.sided", log = TRUE))
  }

  return(log_lik)
}

# additional tools
.fitting_priority       <- function(models){

  # model fitting difficulty using the following heuristic:
  # selection models > random effects | PET/PEESE > non-null models
  fitting_difficulty <- sapply(models, function(model){

    diffuculty <- 0

    if(is.prior.simple(model$priors[["mu"]])){
      diffuculty <- diffuculty + 1
    }
    if(is.prior.simple(model$priors[["tau"]])){
      diffuculty <- diffuculty + 3
    }
    if(is.null(model$priors[["PET"]])){
      diffuculty <- diffuculty + 1
    }else if(is.null(model$priors[["PEESE"]])){
      diffuculty <- diffuculty + 1
    }else if(is.null(model$priors[["omega"]])){
      diffuculty <- diffuculty + 5
    }
  })

  return(order(fitting_difficulty, decreasing = TRUE))
}
.runjags_summary_list   <- function(fit, priors, priors_scale, warnings){

  summary_list <- list()

  for(measure in c("d", "r", "z", "logOR")){

    # prepare transformations if necessary
    if(measure != priors_scale){
      transformations <- list()
      if("mu" %in% names(priors) && ((is.prior.point(priors[["mu"]]) && priors[["mu"]][["parameters"]][["location"]] != 0) || !is.prior.point(priors[["mu"]]))){
        transformations[["mu"]] <- list(
          "fun" = .transform_mu,
          "arg" = list(
            "from" = priors_scale,
            "to"   = measure
          )
        )
      }
      if("tau" %in% names(priors) && ((is.prior.point(priors[["mu"]]) && priors[["mu"]][["parameters"]][["location"]] != 0) || !is.prior.point(priors[["tau"]]))){
        transformations[["tau"]] <- list(
          "fun" = .scale,
          "arg" = list(
            "from" = priors_scale,
            "to"   = measure
          )
        )
      }
      if("PEESE" %in% names(priors)){
        # the transformation is inverse for PEESE
        transformations[["tau"]] <- list(
          "fun" = .scale,
          "arg" = list(
            "from" = measure,
            "to"   = priors_scale
          )
        )
      }
      if(length(transformations) == 0){
        transformations <- NULL
      }
    }else{
      transformations <- NULL
    }

    summary_list[[measure]] <- BayesTools::runjags_estimates_table(
      fit             = fit,
      prior_list      = priors,
      transformations = transformations,
      warnings        = warnings,
      footnotes       = .scale_note(priors_scale, measure),
    )
  }

  return(summary_list)
}


# JAGS tools for model building and marginal likelihood
.JAGS_transformation    <- function(from, to, from_par, to_par_name){

  return(paste0(to_par_name, " = ", .JAGS_transform(from, to, from_par), "\n"))

}
.JAGS_transform         <- function(from, to, par){

  if(from == to){
    return(par)
  }else{
    return(paste0(from, "2", to, "(", par,")"))
  }

}
.JAGS_transformation_se <- function(from, to, from_par, from_anc, to_par_name){

  return(paste0(to_par_name, " = ",.JAGS_transform_se(from, to, from_par, from_anc), "\n"))

}
.JAGS_transform_se      <- function(from, to, par, par_anc){

  if(from == to){
    return(par)
  }else{
    return(paste0("se_", from, "2", "se_", to, "(", par, if(!all(c(from, to) %in% c("logOR", "d"))) paste0(", ", par_anc),")"))
  }

}
.JAGS_scale             <- function(from, to, from_par, to_par_name){

  if(from == to){
    return(paste0(to_par_name, " = ", from_par, "\n"))
  }else{
    return(paste0(to_par_name, " = ", "scale_", from, "2", to, "(", from_par,")", "\n"))
  }

}
