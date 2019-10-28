#' Posterior distribution parameters for mixture priors
#'
#' \code{weights_formulate} calculates posterior distribution paramters for a Binomial likelihood and mixture Beta prior
#'
#' @param w Weight on the informative prior (Distribution coming from a belief)
#' @param x Total Number of events
#' @param n Trial size
#' @param a1 alpha parameter for informative distribution, must be input properly when \eqn{w = 0 or 1}
#' @param b1 beta parameter for informative distribution, must be input properly properly when \eqn{w = 0 or 1}
#' @param a2 alpha parameter for non informative distribution, will not be used if \eqn{w = 0 or 1}
#' @param b2 beta parameter for non informative distribution, will not be used if \eqn{w = 0 or 1}
#'
#' @details
#' When w takes the value 0 or 1, which implies existence of a single Prior Beta distribution, \code{weights_formulate} by default takes a1 and b1 as model parameters.
#' @return \code{weights}  Weights on posterior distribution components. The first value refers to the weight for informative part of the posterior distribution.
#' @return \code{param_inform}  Parameters (alpha, beta) for the informative distribution
#' @return \code{param_noninform}  Parameters (alpha, beta) for the non informative distribution
#' @export
#'
#' @examples weights_formulate(w = 1, x = 1, n = 3, a1 = 1, b1= 1, a2 = 1, b2 = 1)
weights_formulate = function(w = 1, x, n, a1 = 1, b1= 1, a2 = NULL, b2 = NULL)
{
  #Checking the value of weight

  if (w > 1)
  {
    stop("Weight on informative prior can be at most 1")
  }
  else if (w < 0)
  {
    stop("Weight on a prior can not be negative")
  }

  #Checking the eligibility of the parameters

  if (sum(c(a1, b1, a2, b2) <= 0) > 0)
  {
    stop("Beta parameters must be non-negative")
  }

  #Checking the number of events happened is less than total number of trials

  if (n < 1)
  {
    stop("The trial size must be at least 1")
  }

  if(x > n)
  {
    stop("Number of successes for the event (i.e. experiencing DLT 's) must be lower than total number of trials (i.e. patients treated)")
  }

  #Checking feasibility condition of prior parameters
  a1_null = is.null(a1)
  b1_null = is.null(b1)
  a2_null = is.null(a2)
  b2_null = is.null(b2)
  total_null = a1_null + b1_null + a2_null + b2_null

  if (total_null == 4)
  {
    stop("Please input a1, a2, b1, b2 properly. ")
  }

  if(w %in% c(0, 1))
  {
    if (total_null == 2)
    {
      if((a2_null + b2_null) == 1)
      {
        stop("Please input either both a1 and b1, or both a2 and b2, (ai,bi) is the pair of parameters. For Uniform Distribution, either put a1 = 1, b1 = 1, or put, a2 = 1 and b2 = 1")
      }
      else if ((a2_null + b2_null) == 0)
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
    else if (total_null == 0)
    {
      warning("Check inputs for prior parameters, taking a1 and b1 as original parameters")
    }
  }
  else
  {
    if (total_null > 0 )
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
