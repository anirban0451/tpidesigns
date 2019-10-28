#' Calculation of Unit Probability Mass
#'
#' \code{UPM} calculates Unit Probability Mass for an interval (a, b) when the Underlying distribution is beta or mixture of two beta distributions.
#' @importFrom stats pbeta
#' @param w Weight on the first component of mixture distribution, i.e, the informative Prior
#' @param a,b Range Parameters
#' @param a1,b1 alpha and beta parameters for informative prior component, must be input when w equals 1
#' @param a2,b2 alpha and beta parameters for noninformative prior component, must be input when w equals 0
#'
#' @details
#' UPM(a,b) = \eqn{(F(b) - F(a))/(b - a)}, defined for an interval (a,b), when X~F().
#' In this function, F() is assumed to be Cumulative Beta distribution function or mixture of two cumulative Beta distribution functions.
#' If F() consists of a single Beta distribution, and not a mixture, then we must either input \eqn{w = 1} and a1, b1 , or \eqn{w = 0} and a2,b2
#' @return Unit Probability Mass value or the UPM value
#' @export
#'
#' @examples
#' UPM(w = 1, a = 0.3, b = 0.4, a1 = 2, b1 = 5)
#' UPM(w = 0, a = 0.3, b = 0.4, a2 = 2, b2 = 5)
#' UPM(w = 0.3, a = 0.3, b = 0.4, a1 = 3, b1 = 6, a2 = 2, b2 = 5)
UPM <- function(w, a = 0, b = 1, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  #Checking if the weight value is at most 1 or at least 0
  if(w < 0 || w >1)
  {
    stop("w is weight taken on first prior (informative), which can lie between 0 and 1")
  }

  #Checking the feasibility of the domain
  if(a < 0)
  {
    a = 0
    warning("Domain of Beta distribution is (0,1), changing a to 0")
  }
  if(b > 1)
  {
    b = 1
    warning("Domain of Beta distribution is (0,1), changing a to 1")
  }
  if(a >= b) stop("a must be less than b and both should lie within (0,1")
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
  #Checking the over toxicity of the dose
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
      stop("Please input model parameters properly")
    }
  }
  #calculating the value

  if (w %in% c(0,1))
  {
    val = (pbeta(b, shape1 = a1, shape2 = b1) - pbeta(a, shape1 = a1, shape2 = b1)) / (b - a)
  }

  else
  {
    val =  w * ((pbeta(b, shape1 = a1, shape2 = b1) - pbeta(a, shape1 = a1, shape2 = b1)) / (b - a)) +
      (1 - w) * ((pbeta(b, shape1 = a2, shape2 = b2) - pbeta(a, shape1 = a2, shape2 = b2)) / (b - a))
  }

  return(val)
}
