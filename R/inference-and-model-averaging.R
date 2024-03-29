.balance_component_probability   <- function(object){

  converged <- .get_model_convergence(object)
  if(all(!converged))
    stop("All models included in the ensemble failed to converge.")

  # assess the main component type
  component_types <- cbind.data.frame(
    effect         = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "effect")),
    heterogeneity  = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "heterogeneity")),
    bias           = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "bias")),
    baseline       = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "baseline")),
    hierarchical   = sapply(object[["models"]], function(model) !.is_component_null(model[["priors"]], "hierarchical"))
  )
  # add regressions if neccessary
  if(!is.null(object$add_info[["predictors_test"]])){
    for(predictor_test in object$add_info[["predictors_test"]]){
      component_types[[predictor_test]] <- sapply(object[["models"]], function(model) predictor_test %in% model[["terms_test"]])
    }
  }
  # extract the prior odds set by user
  prior_weights  <- sapply(object[["models"]], function(model) model[["prior_weights_set"]])

  # check whether there is a comparable model for each non-converged models
  for(i in seq_along(object[["models"]])[!converged]){

    temp_same_structure <- sapply(seq_along(object$models)[converged], function(j) all(component_types[j,] == component_types[i,]))
    temp_same_structure <- seq_along(object$models)[converged][temp_same_structure]

    # if yes, transfer the prior odds
    if(length(temp_same_structure) >= 1){

      prior_weights[temp_same_structure] <- prior_weights[temp_same_structure] + prior_weights[i] / length(temp_same_structure)
      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge. However, there were other models with the same combination of model components and their prior probability was increased to account for the failed models.")

    }else{

      prior_weights[i] <- 0
      object$add_info[["warnings"]] <- c(object$add_info[["warnings"]], "Some of the models failed to converge and their prior probability couldn't be balanced over models with the same combination of model components since they don't exist.")

    }
  }

  for(i in seq_along(object[["models"]])){
    object[["models"]][[i]][["prior_weights"]] <- prior_weights[i]
  }

  return(object)
}
.restore_component_probability   <- function(object){

  # extract the prior odds set by user
  prior_weights  <- sapply(object[["models"]], function(model) model[["prior_weights_set"]])

  for(i in seq_along(object[["models"]])){
    object[["models"]][[i]][["prior_weights"]] <- prior_weights[i]
  }

  return(object)

}
.ensemble_inference    <- function(object){

  # use only converged models with prior weights > 0 for inference about parameters
  prior_weights <- sapply(object[["models"]], function(model) model[["prior_weights"]])
  models        <- object[["models"]][.get_model_convergence(object) & prior_weights > 0]

  # obtain the component type
  effect         <- sapply(models, function(model)!.is_component_null(model[["priors"]], "effect"))
  heterogeneity  <- sapply(models, function(model)!.is_component_null(model[["priors"]], "heterogeneity"))
  bias           <- sapply(models, function(model)!.is_component_null(model[["priors"]], "bias"))
  hierarchical   <- sapply(models, function(model)!.is_component_null(model[["priors"]], "hierarchical"))
  baseline       <- sapply(models, function(model)!.is_component_null(model[["priors"]], "baseline"))

  # obtain the parameter types
  weightfunctions <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.weightfunction)))
  PET             <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PET)))
  PEESE           <- sapply(models, function(model)any(sapply(model[["priors"]], is.prior.PEESE)))

  # define inference options: always effect and heterogeneity
  components      <- c("Effect", "Heterogeneity")
  parameters      <- c("mu", "tau")
  components_null <- list("Effect" = !effect, "Heterogeneity" = !heterogeneity)
  parameters_null <- list("mu"     = !effect, "tau"           = !heterogeneity)

  if(any(bias)){
    components      <- c(components,      "Bias")
    components_null <- c(components_null, "Bias" = list(!bias))
  }
  if(any(hierarchical)){
    components      <- c(components,      "Hierarchical")
    parameters      <- c(parameters,      "rho")
    components_null <- c(components_null, "Hierarchical" = list(!hierarchical))
    parameters_null <- c(parameters_null, "rho" = list(!hierarchical))
  }
  if(any(baseline)){
    components      <- c(components,      "Baseline")
    components_null <- c(components_null, "Baseline" = list(!baseline))
  }
  if(any(weightfunctions)){
    components      <- c(components,      "bias.selection-models")
    parameters      <- c(parameters,      "omega")
    components_null <- c(components_null, "bias.selection-models" = list(!(weightfunctions & bias)))
    parameters_null <- c(parameters_null, "omega" = list(!weightfunctions))
  }
  if(any(PET) & any(PEESE)){
    components      <- c(components,      "bias.PET-PEESE")
    parameters      <- c(parameters,      "PET", "PEESE")
    components_null <- c(components_null, "bias.PET-PEESE" = list(!((PET | PEESE) & bias)))
    parameters_null <- c(parameters_null, "PET" = list(!PET), "PEESE" = list(!PEESE))
  }else if(any(PET)){
    components      <- c(components,      "bias.PET")
    parameters      <- c(parameters,      "PET")
    components_null <- c(components_null, "bias.PET" = list(!(PET & bias)))
    parameters_null <- c(parameters_null, "PET" = list(!PET))
  }else if(any(PEESE)){
    components      <- c(components,      "bias.PEESE")
    parameters      <- c(parameters,      "PEESE")
    components_null <- c(components_null, "bias.PEESE" = list(!(PEESE & bias)))
    parameters_null <- c(parameters_null, "PEESE" = list(!PEESE))
  }

  # deal with meta-regression
  if(!is.null(object[["formula"]])){
    # use the intercept for the effect
    parameters[parameters == "mu"]                         <- "mu_intercept"
    names(parameters_null)[names(parameters_null) == "mu"] <- "mu_intercept"

    # add the terms
    model_predictors      <- lapply(models, function(model) model[["terms"]])
    model_predictors_test <- lapply(models, function(model) model[["terms_test"]])

    predictors      <- object$add_info[["predictors"]]
    predictors_test <- object$add_info[["predictors_test"]]

    # define inference options
    components_predictors      <- NULL
    parameters_predictors      <- "mu_intercept"
    parameters_marginal        <- "mu_intercept"
    components_predictors_null <- list()
    parameters_predictors_null <- list("mu_intercept" = !effect)
    parameters_marginal_null   <- list("mu_intercept" = !effect)

    components_predictors_distributions      <- NULL
    components_predictors_distributions_null <- list()


    # predictors
    for(i in seq_along(predictors_test)){
      components_predictors <- c(components_predictors, .BayesTools_parameter_name(predictors_test[i]))
      components_predictors_null[[.BayesTools_parameter_name(predictors_test[i])]] <-
        sapply(model_predictors_test, function(x) if(length(x) == 0) TRUE else !(predictors_test[i] %in% x))
    }

    for(i in seq_along(predictors)){
      parameters_predictors <- c(parameters_predictors, .BayesTools_parameter_name(predictors[i]))
      parameters_marginal   <- c(parameters_marginal,   .BayesTools_parameter_name(predictors[i]))
      parameters_predictors_null[[.BayesTools_parameter_name(predictors[i])]] <-
        sapply(model_predictors_test, function(x) !(predictors[i] %in% x))
      parameters_marginal_null[[.BayesTools_parameter_name(predictors[i])]] <-
        sapply(seq_along(models), function(j) !((predictors[i] %in% model_predictors_test[j]) || effect[j]))
    }


    # get models inference
    if(is.null(components_predictors)){
      inference_predictors <- NULL
    }else{
      inference_predictors <- BayesTools::ensemble_inference(
        model_list   = models,
        parameters   = components_predictors,
        is_null_list = components_predictors_null,
        conditional  = FALSE
      )
    }
    # deal with the possibility of only null models models
    if(all(sapply(parameters_predictors_null, all))){
      inference_predictors_conditional <- NULL
    }else{
      inference_predictors_conditional <- BayesTools::ensemble_inference(
        model_list   = models,
        parameters   = parameters_predictors[!sapply(parameters_predictors_null, all)],
        is_null_list = parameters_predictors_null[!sapply(parameters_predictors_null, all)],
        conditional  = TRUE
      )
    }

    # get model-averaged posteriors
    if(is.null(parameters_predictors)){
      posteriors_predictors <- NULL
    }else{
      posteriors_predictors <- BayesTools::mix_posteriors(
        model_list   = models,
        parameters   = parameters_predictors,
        is_null_list = parameters_predictors_null,
        seed         = object$add_info[["seed"]],
        conditional  = FALSE
      )
      posteriors_predictors <- BayesTools::transform_factor_samples(posteriors_predictors)
    }

    # deal with the possibility of only null models models
    if(all(sapply(parameters_predictors_null, all))){
      posteriors_predictors_conditional <- NULL
    }else{
      posteriors_predictors_conditional <- BayesTools::mix_posteriors(
        model_list   = models,
        parameters   = parameters_predictors[!sapply(parameters_predictors_null, all)],
        is_null_list = parameters_predictors_null[!sapply(parameters_predictors_null, all)],
        seed         = object$add_info[["seed"]],
        conditional  = TRUE
      )
      posteriors_predictors_conditional <- BayesTools::transform_factor_samples(posteriors_predictors_conditional)
    }

    # create marginal estimates and summary
    if(all(sapply(parameters_marginal_null, all))){
      inference_marginal <- NULL
    }else{
      inference_marginal <- BayesTools::marginal_inference(
        model_list          = models,
        marginal_parameters = parameters_marginal[!sapply(parameters_marginal_null, all)],
        parameters          = parameters_marginal,
        is_null_list        = parameters_marginal_null,
        formula             = object[["formula"]],
        seed                = object$add_info[["seed"]],
        silent              = TRUE
      )
    }
  }else{
    # create empty objects in case of no predictors
    inference_predictors              <- NULL
    inference_predictors_conditional  <- NULL
    posteriors_predictors             <- NULL
    posteriors_predictors_conditional <- NULL
    inference_marginal                <- NULL
  }

  ### get models inference
  inference <- BayesTools::ensemble_inference(
    model_list   = models,
    parameters   = components,
    is_null_list = components_null,
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    inference_conditional <- NULL
  }else{
    inference_conditional <- BayesTools::ensemble_inference(
      model_list   = models,
      parameters   = components[!sapply(components_null, all)],
      is_null_list = components_null[!sapply(components_null, all)],
      conditional  = TRUE
    )
  }


  ### get model-averaged posteriors
  posteriors <- BayesTools::mix_posteriors(
    model_list   = models,
    parameters   = parameters,
    is_null_list = parameters_null,
    seed         = object$add_info[["seed"]],
    conditional  = FALSE
  )
  # deal with the possibility of only null models models
  if(all(sapply(components_null, all))){
    posteriors_conditional <- NULL
  }else{
    posteriors_conditional <- BayesTools::mix_posteriors(
      model_list   = models,
      parameters   = parameters[!sapply(parameters_null, all)],
      is_null_list = parameters_null[!sapply(parameters_null, all)],
      seed         = object$add_info[["seed"]],
      conditional  = TRUE
    )
  }

  # rename mu_intercept back to mu
  if(any(names(posteriors) == "mu_intercept")){
    attr(posteriors[["mu_intercept"]], "parameter")        <- "mu"
    names(posteriors)[names(posteriors) == "mu_intercept"] <- "mu"
  }
  if(any(names(posteriors_conditional) == "mu_intercept")){
    attr(posteriors_conditional[["mu_intercept"]], "parameter")                    <- "mu"
    names(posteriors_conditional)[names(posteriors_conditional) == "mu_intercept"] <- "mu"
  }

  # return the results
  output <- list(
    inference                         = inference,
    inference_conditional             = inference_conditional,
    inference_predictors              = inference_predictors,
    inference_predictors_conditional  = inference_predictors_conditional,
    inference_marginal                = inference_marginal,
    posteriors                        = posteriors,
    posteriors_conditional            = posteriors_conditional,
    posteriors_predictors             = posteriors_predictors,
    posteriors_predictors_conditional = posteriors_predictors_conditional
  )
  return(output)
}

.compute_coeficients   <- function(RoBMA){
  if(!is.null(RoBMA[["posteriors_predictors"]])){
    coefs <- do.call(c, unname(lapply(RoBMA[["posteriors_predictors"]], function(posterior){
      if(inherits(posterior, "mixed_posteriors.factor")){
        out        <- apply(posterior, 2, mean)
        names(out) <- .output_parameter_names(names(out))
      }else{
        out        <- mean(posterior)
        names(out) <- .output_parameter_names(attr(posterior,"parameter"))
      }
      return(out)
    })))
  }else{
    coefs        <- c("mu" = mean(RoBMA$posteriors[["mu"]]))
  }
  return(c(
    coefs,
    "tau"    = mean(RoBMA$posteriors[["tau"]]),
    "rho"    = if(!is.null(RoBMA$posteriors[["rho"]]))   mean(RoBMA$posteriors[["rho"]]),
    if(!is.null(RoBMA$posteriors[["omega"]])) apply(RoBMA$posteriors[["omega"]], 2, mean),
    "PET"    = if(!is.null(RoBMA$posteriors[["PET"]]))   mean(RoBMA$posteriors[["PET"]]),
    "PEESE"  = if(!is.null(RoBMA$posteriors[["PEESE"]])) mean(RoBMA$posteriors[["PEESE"]])
  ))
}
