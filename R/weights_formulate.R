#' Posterior distribution parameters for mixture priors of Beta distribution
#'
#' \code{weights_formulate} calculates posterior distribution paramters when a Binomial Likelihood is given and mixture of two Beta distribution
#' is taken as prior for the parameter p in the Binomial Distribution
#'
#' @param w Weight on the first Beta distribution of the mixture Prior
#' @param x Total Number of events (In Dose Escalation Oncology Trials, this may be defined as
#' number of people who have experienced Dose Limiting Toxicities through administration of current Dose Level)
#' @param n Trial size (In Dose Escalation Oncology Trials, this may be defined as
#' total number of people who have been administered current Dose Level (missing responses will be excluded). Necessarily n will be greater than or equal to x
#' @param a1 alpha parameter ( > 0) for 1st Beta distribution, must be input properly when \eqn{w = 0 or 1}
#' @param b1 beta parameter ( > 0) for 1st Beta distribution, must be input properly properly when \eqn{w = 0 or 1}
#' @param a2 alpha parameter ( > 0) for 2nd Beta distribution, will not be used if \eqn{w = 0 or 1}
#' @param b2 beta parameter ( > 0) for 2nd Beta distribution, will not be used if \eqn{w = 0 or 1}
#'
#' @details
#' While using the function for Oncology Dose Escalation Trials, w is assumed to be the weight taken on the 1st Prior Distribution of the mixture,
#' named as Informative Prior, that is, the Beta Distribution aring from the past studies. Hence, (1 - w)
#' represents the weight on the 2nd part of the Beta distribution. When w takes the value 0 or 1(implies exstence of a single Prior Beta distribution),
#' \code{weights_formulate} by default takes a1 and b1 as model parameters and ignores a2 and b2.
#' @details
#' Let, \code{x|p ~ Binom(n,p) hence, f(x|p) = choose(n, x) * p^x (1 - p)^(n-x)} and \code{ p ~ g(p) = w * dbeta(p, a1, b1) + (1 - w) * dbeta(p, a2, b2)} is the prior for p
#' @details
#' Then the unconditional distribution of x is given by,
#' @details
#' \eqn{P(X = x) = f(x) = choose(n,x) * [w * (beta(a1 + x, b1 + n - x)/beta(a1, b1)) + (1 - w) * (beta(a2 + x, b2 + n - x)/beta(a2, b2)) ]}
#' @details
#' The prior distribution for p becomes,
#' @details
#' \eqn{g(p|x) = (f(x | p) * g(p) / f(x)) = w_1 * dbeta(p, a1 + x, b1 + n - x) + (1 - w_1) * dbeta(p, a2 + x, b2 + n - x)}
#' @details
#' where, \eqn{w1 = (w * beta(a1 + x , b1 + n - x) / beta(a1, b1))/((w * beta(a1 + x , b1 + n - x) / beta(a1, b1)) + ((1 - w) * beta(a2 + x , b2 + n - x) / beta(a2, b2))) }
#' @details
#' Please remember that, dbeta(p, a, b) refers to the pdf of Beta distribution of p with parameters a and b
#' , while beta(a,b) gives us the value of beta function with parameters a and b
#' @return \code{weights}  Weight on the 1st part of the Posterior Mixture Beta Distribution. When there is only one Beta distribution,
#' this value will always return 0 (or 1) , if we pass the value of w as 0 (or 1)
#' @return \code{param_inform}  Parameters (alpha, beta) for the 1st Beta distribution. When there is only one prior,
#' these values are returned.
#' @return \code{param_noninform}  Parameters (alpha, beta) for the 2nd Beta distribution. When there is only one prior,
#' these values are returned as NULL
#' @seealso
#' \code{\link[base]{Special}} to know about beta(a,b) function
#' @seealso
#' \code{\link[stats]{Beta}} to know about beta distribution function
#' @seealso
#' \code{\link[stats]{Binomial}} to know about Binomial Distribution function
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Beta-binomial_distribution} to know about Beta Binomial distribution, the unconditional distribution of x
#' @seealso
#' \url{https://www.cancer.gov/publications/dictionaries/cancer-terms/def/dose-limiting} to know about Dose Limiting Toxicity (DLT)
#'
#' @export
#'
#' @examples weights_formulate(w = 1, x = 1, n = 3, a1 = 1, b1= 1, a2 = 1, b2 = 1) #Will show an warning but return values
#' @examples weights_formulate(w = 0.1, x = 1, n = 3, a1 = 1, b1= 1, a2 = 1, b2 = 1)
weights_formulate = function(w = NULL, x, n, a1 = NULL, b1= NULL, a2 = NULL, b2 = NULL)
{
  #Checking the value of weight

  if (isTRUE(w > 1))
  {
    stop("Weight on informative prior can be at most 1")
  }
  else if (isTRUE(w < 0))
  {
    stop("Weight on a prior can not be negative")
  }

  #Checking the eligibility of the parameters

  if (isTRUE(sum(c(a1, b1, a2, b2) <= 0) > 0))
  {
    stop("Beta parameters must be non-negative")
  }

  #Checking the number of events happened is less than total number of trials

  if (isTRUE(n < 1))
  {
    stop("The trial size must be at least 1")
  }

  if(isTRUE(x > n))
  {
    stop("Number of successes for the event (i.e. experiencing DLT 's) must be lower than total number of trials (i.e. patients treated)")
  }

  #Checking feasibility condition of prior parameters
  a1_null = is.null(a1)
  b1_null = is.null(b1)
  a2_null = is.null(a2)
  b2_null = is.null(b2)
  total_null = a1_null + b1_null + a2_null + b2_null

  if (isTRUE(total_null == 4))
  {
    stop("Please input a1, a2, b1, b2 properly. ")
  }

  if(w %in% c(0, 1))
  {
    if (isTRUE(total_null == 2))
    {
      if(isTRUE((a2_null + b2_null) == 1))
      {
        stop("Please input either both a1 and b1, or both a2 and b2, (ai,bi) is the pair of parameters. For Uniform Distribution, either put a1 = 1, b1 = 1, or put, a2 = 1 and b2 = 1")
      }
      else if (isTRUE((a2_null + b2_null) == 0))
      {
        a1 = a2
        b1 = b2
        warning("You should put the parameter values for a1 and b1 instead of a2 and b2")
      }
    }
    else if (total_null %in% c(1,3))
    {
      stop("Please input a1, b1, a2, b2 properly, (ai,bi) is the pair of parameters. For Uniform Distribution, either put a1 = 1, b1 = 1, or put, a2 = 1 and b2 = 1")
    }
    else if (isTRUE(total_null == 0))
    {
      warning("Check inputs for prior parameters, taking a1 and b1 as original parameters")
    }
  }
  else
  {
    if (isTRUE(total_null > 0) )
    {
      stop("Please input model parameters  for both priors properly")
    }
  }

  #Calculating parameters based on extreme weights

  if (w %in% c(0, 1))
  {
    param_inform = c(a1 + x, b1 + n - x)
    param_noninform = c(NULL, NULL)
  }

  #tackling general case

  else
  {
    param_inform = c(a1 + x, n - x + b1)
    param_noninform = c(a2 + x, n - x + b2)
    w = (w  * (beta(a1 + x, b1 + n - x) / beta(a1, b1))) /
      (w * (beta(a1 + x, b1 + n - x) / beta(a1, b1)) +
         (1 - w) * (beta(a2 + x, b2 + n - x) / beta(a2, b2)))
  }

  ls = list()
  ls$weight = w
  ls$param_inform = param_inform
  ls$param_noninform = param_noninform
  return(ls)
}
