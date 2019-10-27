#' Posterior distribution parameters for mixture priors
#'
#' \code{weights_formulate} calculates posterior distribution paramters for a Binomial likelihood and mixture Beta prior
#'
#' @param w Weight on the informative prior (Distribution coming from a belief)
#' @param x Total Number of events
#' @param n Trial size
#' @param a1 alpha parameter for informative distribution
#' @param b1 beta parameter for informative distribution
#' @param a2 alpha parameter for non informative distribution
#' @param b2 beta parameter for non informative distribution
#'
#' @return \code{weights}  Weights on posterior distribution components. The first value refers to the weight for informative part of the posterior distribution.
#' @return \code{param_inform}  Parameters (alpha, beta) for the informative distribution
#' @return \code{param_noninform}  Parameters (alpha, beta) for the non informative distribution
#' @export
#'
#' @examples weights_formulate(w = 1, x = 1, n = 3, a1 = 1, b1= 1, a2 = 1, b2 = 1)
weights_formulate = function(w = 1, x, n, a1 = 1, b1= 1, a2 = 1, b2 = 1)
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

  if (a1 <= 0 || b1 <= 0 || a2 <= 0 || b2 <= 0)
  {
    stop("Beta parameters can not be negative")
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

  #Calculating parameters based on extreme weights

  if (w == 1)
  {
    param_inform = c(a1 + x, b1 + n - x)
    param_noninform = c(0, 0)
  }

  else if (w == 0)
  {
    param_inform = c(0,0)
    param_noninform = c(a2 + x, b2 + n - x)
  }

  #tackling general case

  else
  {
    param_inform = c(a1 + x, n - x + b1)
    param_noninform = c(a2 + x, n - x + b2)
    w = (w *( beta(a1 + x, b1 + n - x)/beta(a1, b1)))/
      (w *( beta(a1 + x, b1 + n - x)/beta(a1, b1)) +
         (1 - w) *( beta(a2 + x, b2 + n - x)/beta(a2, b2)))
  }

  ls = list()
  ls$weights = c(w, 1 - w)
  ls$param_inform = param_inform
  ls$param_noninform = param_noninform
  return(ls)
}
