#' Calculation of Unit Probability Mass
#'
#' \code{UPM} calculates Unit Probability Mass for an interval (a, b) when the Underlying distribution is beta or mixture of two beta distributions.
#' @import stats
#' @param a,b Range Parameters between which UPM is needed to be calculated.
#' @inheritParams weights_formulate
#'
#' @details
#' Unit Probability MASS or UPM(a,b) = \eqn{(F(b) - F(a))/(b - a)}, defined for an interval (a,b), when X~F().
#' In this function, F() is assumed to be Cumulative Beta distribution function or mixture of two cumulative Beta distribution functions.
#' @details
#' Hence, \eqn{F(x) =  w * pbeta(x, a1, b1) + (1 - w) * pbeta(x, a2, b2)}, pbeta is cumulative Beta distribution.
#' @details
#' If F() consists of a single Beta distribution, and not a mixture, then the convention here assumed is
#' to input \eqn{w = 1} and a1, b1 , or \eqn{w = 0} and a2,b2
#' @return Unit Probability Mass value or the UPM value for the interval (a, b)
#' @seealso
#' \code{\link{weights_formulate}}, \code{\link[stats]{Beta}}
#' @export
#'
#' @examples UPM(w = 1, a = 0.3, b = 0.4, a1 = 2, b1 = 5)
#' @examples UPM(w = 0, a = 0.3, b = 0.4, a2 = 2, b2 = 5)
#' @examples UPM(w = 0.3, a = 0.3, b = 0.4, a1 = 3, b1 = 6, a2 = 2, b2 = 5)
#' @examples UPM(w = 1, a = 0.3, b = 0.4, a1 = 2, b1 = 5, a2 = 7, b2 = 8) #will give warning
UPM <- function(w, a = 0, b = 1, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  #Checking if the weight value is at most 1 or at least 0
  if(isTRUE(w < 0 || w > 1))
  {
    stop("w is weight taken on first prior (informative), which can lie between 0 and 1")
  }

  #Checking the feasibility of the domain
  if(isTRUE(a < 0))
  {
    a = 0
    warning("Domain of Beta distribution is (0,1), changing a to 0")
  }
  if(isTRUE(b > 1))
  {
    b = 1
    warning("Domain of Beta distribution is (0,1), changing a to 1")
  }
  if(isTRUE(a >= b)) stop("a must be less than b and both should lie within (0,1")
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
  #Checking the over toxicity of the dose
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
    if (isTRUE(total_null > 0))
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



#' Graphical plot of Unit Probability MASS
#'
#' \code{upmplot} Produces a graphical plot of Unit Probability Mass for a given set of parameters.
#' @import ggplot2
#' @inheritParams weights_formulate
#' @param pt Target toxicity proportion to achieve in current Dose Level (Less Toxicity means under- dosing, where as more toxicity means over - dosing)
#' @param e1 Amount of variation that can be allowed to the left of the pt value to conclude that target toxicity has been achieved.
#' Default value is 0.05
#' @param e2 Amount of variation that can be allowed to the right of the pt value to conclude that target toxicity has been achieved.
#' Default value is 0.05
#' @param design The Design that is implemented in the trials. This arguement includes values "mtpi" and "mmtpi"
#'
#' @return A graph that includes Probability Distributions of the Dose Limiting Toxocity Rate and value of Unit Probability Mass at corresponding intervals.
#' @inherit UPM details
#' @section Decision Making Based on UPM values:
#' For modified Toxicity Probability Interval (mTPI) Design, the toxicity range (0,1) is divided into
#' three ranges, (1) Under-Dosing Interval [0, pt - e1), (2) Target-Toxicity Interval [pt - e1, pt - e2], (3) Over-Dosing Interval (pt + e2, 1].
#' UPM is calculated for the the above intervals and Decision is taking accordingly,\cr if the UPM is maximum for interval (1),
#' then the strength of the current Dosage is escalated,\cr if its maximum for Interval (2), then more patients are administered with
#' current dose,\cr if the UPM is maximum in interval (3), then strength of the current Dose is de-escalated.\cr For Modified Toxicity Interval Design-2 (mTPI -2, encoded as "mmtpi")
#' the intervals (1) and (3) are again divided into another sub- intervals and same steps are followed.\cr But, before that, we must ensure that the Dose is not severely toxic
#' and hence it is advised to run the \code{\link{decisiontpi}} function to know about the severity of current Dose.
#' @seealso
#' \code{\link{UPM}}, \code{\link{weights_formulate}}
#' @export
#'
#' @examples require(ggplot2)
#' @examples n = 13 #must be a value >= 3
#' @examples x = sample.int(n, 1)
#' @examples upmplot(x = 5, n = 7, pt = 0.3, design = "mmtpi", w = 0.1, a1 = 1, a2 = 1, b1 = 4, b2 = 6)
upmplot <- function(x , n , pt, e1 = 0.05, e2 = 0.05, design = c("mtpi", "mmtpi"), w, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  if(isTRUE(pt > 1 || pt < 0))
  {
    stop("Target toxicity Probability should take values between 0 and 1")
  }
  if(isTRUE(pt - e1 < 0 || pt + e2 > 1))
  {
    stop ("e1 and e2, two thresholds should be small compared to the target probability pt")
  }
  if (isTRUE(w > 1))
  {
    stop("Weight on informative prior can be at most 1")
  }
  else if (isTRUE(w < 0))
  {
    stop("Weight on a prior can not be negative")
  }

  #Checking the eligibility of the parameters

  if (isTRUE(any(c(a1, b1, a2, b2) <= 0) == TRUE))
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

  if (design == "mtpi")
  {
    interval = c(0, pt - e1, pt + e2, 1)
    length_interval = length(interval)
  }
  else if (design %in% c("mtpi", "mmtpi"))
  {
    breaks_lower = floor((pt - e1) / 0.1)
    breaks_upper = floor((1 - pt - e2)  / 0.1)
    interval = c(0, pt - e1 - 0.1 * (breaks_lower : 0) , pt + e2 + 0.1 * (0 : breaks_upper) , 1)
    length_interval = length(interval)
  }
  else
  {
    stop("Please input one input among the designs: mtpi, mmtpi")
  }
  params = weights_formulate(w = w, x = x, n = n, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
  w = params$weight
  a1 = params$param_inform[1]
  b1 = params$param_inform[2]
  a2 = params$param_noninform[1]
  b2 = params$param_noninform[2]
  upm_array= rep(0, length_interval - 1)
  for(i in 1: (length_interval - 1))
  {
    upm_array[i] = UPM(w = w, a = interval [i], b = interval [i + 1], a1 = a1, b1 = b1, a2 = a2, b2 = b2)
  }
  if (w %in% c(0,1))
  {
    plotupm = ggplot(data.frame(x=seq(0.01,1,0.01)), aes(x)) +
      stat_function(fun=function(x) dbeta(x, shape1 = a1, shape2 = b1))
  }
  else
  {
    plotupm = ggplot(data.frame(x=seq(0.01,1,0.01)), aes(x)) +
      stat_function(fun=function(x) w * dbeta(x, shape1 = a1, shape2 = b1) + (1 - w) * dbeta(x, shape1 = a2, shape2 = b2))
  }
  plotupm_addY <- plotupm + geom_vline(xintercept = interval, linetype="dashed", color = "steelblue", size = 0.7)

  segment_x <- interval[-length_interval]
  segment_xend <- interval[-1]
  segment_y <- segment_yend <- upm_array
  segment_data <- data.frame(x = segment_x, y = segment_y, x_end = segment_xend,y_end = segment_yend)
  plotupm_addsegments <- plotupm_addY +
    geom_segment(data = segment_data, mapping = aes(x = x, xend = segment_data$x_end, y = segment_data$y, yend = segment_data$y_end))

  plotupm_addfootnote <- plotupm_addsegments +
    labs(title = " Plotting of UPM values and Posterior DLT distribution", x = "DLT Rate", y = "Unit Probability Mass (UPM)",
         caption = "The Dashed lines represent the intervals and the Horizontal lines represent the UPM")
  return(plotupm_addfootnote)
}
