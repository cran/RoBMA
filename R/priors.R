#' @name prior
#' @inherit BayesTools::prior
#' @export
prior <- BayesTools::prior

#' @name prior_none
#' @inherit BayesTools::prior_none
#' @export
prior_none <- BayesTools::prior_none

#' @name prior_factor
#' @inherit BayesTools::prior_factor
#' @export
prior_factor <- BayesTools::prior_factor

#' @name prior_PET
#' @inherit BayesTools::prior_PET
#' @export
prior_PET  <- BayesTools::prior_PET

#' @name prior_PEESE
#' @inherit BayesTools::prior_PEESE
#' @export
prior_PEESE <- BayesTools::prior_PEESE

#' @name prior_weightfunction
#' @inherit BayesTools::prior_weightfunction
#' @export
prior_weightfunction <- BayesTools::prior_weightfunction

#' @name prior_informed
#' @inherit BayesTools::prior_informed
#' @details Further details can be found in \insertCite{erp2017estimates;textual}{RoBMA},
#' \insertCite{gronau2017bayesian;textual}{RoBMA}, and
#' \insertCite{bartos2021bayesian;textual}{RoBMA}.
#' @references
#' \insertAllCited{}
#' @export
prior_informed <- BayesTools::prior_informed

#' @title Orthornomal contrast matrix
#'
#' @description Return a matrix of orthornomal contrasts.
#' Code is based on \code{stanova::contr.bayes} and corresponding to description
#' by \insertCite{rouder2012default;textual}{BayesTools}
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.orthonormal <- BayesTools::contr.orthonormal

#' @title Mean difference contrast matrix
#'
#' @description Return a matrix of mean difference contrasts.
#' This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#' distributions on difference between the gran mean and factor level are identical independent
#' of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#' the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#' distribution on difference between each factor level and the grand mean -- this is approximately
#' twice the scale of \code{contr.orthonormal}.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.meandif <- BayesTools::contr.meandif

#' @title Independent contrast matrix
#'
#' @description Return a matrix of independent contrasts -- a level for each term.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
#' @export
contr.independent <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n >= 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(x = 1, nrow = n, ncol = n)

  return(cont)
}
