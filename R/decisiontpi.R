#' Dosing Decisions based on different tpi designs
#'
#' \code{decisiontpi} gives whether to escalate, stay or de-escalate to a level of Dose given total number of patients and total number of people experiencing toxicities in a cohort
#' @param eta threshold value to check if the Dose is severely toxic
#' @param design Design parameter, tells us which design to use. Options are \code{"tpi", "mtpi", "mmtpi"}
#' @inheritParams upmplot
#' @return
#' \code{"E"} We should increase the current level of Dose.\cr
#' \code{"S"} We should stay at the current level of Dose and treat more patients.\cr
#' \code{"D"} We should decrease the current level of Dose.\cr
#' \code{"DU"} We should decrease the current level of Dose and we should not go beyond this Dose level for treating more patients
#' @details
#' \code{decisiontpi} checks if the DLT rate within the sample at current Dose level is severely toxic. At first, it calculates the posterior distribution
#' of the DLT(or, Dose Limiting Toxicity) Rate based on the number of DLT 's in the sample. If the probability is too much (> eta), then the
#' function returns the value "DU", which implies that the Current Dose level is unacceptably toxic and can't ever be used for further administration.\cr
#' If the dose is not severely toxic, then it gives us the decision rule based on the design we provide.\cr
#' If we work with 'tpi' design, then the range of DLT rate, ie, [0,1], is broken up into three intervals, Under-Dosing [0, pt - e1), Target-Toxicity [pt - e1, pt + e2] and Over-Dosing (pt + e2, 1]\cr
#' Then, probability for these intervals are calculated and the interval for which the probability is high, leads to the Decision- making.\cr
#' For mTPI and mTPI-2 (coded as \code{"mmtpi"} in the package), the decision making after ensuring the non-severe toxicity of the current dose level
#' is based on the Unit Probability Mass among the Intervals.
#' @seealso \code{\link{UPM}} for definition on Unit Probability Mass, \code{\link{upmplot}} for the Decision Making Criterion based on UPM
#' @export
#'
#' @examples n = 13 #must be a value >= 3
#' @examples x = sample.int(n, 1)
#' @examples decisiontpi(x = x, n = n, design = "mmtpi", pt = 0.4, e1 = 0.06, e2 = 0.04, eta = 0.95, w = 0.4, a1 = 4, b1 = 3, a2 = 1, b2 = 1)
#'
decisiontpi <- function(pt, e1 = 0.05, e2 = 0.05, x, n, eta, design = c("tpi", "mtpi", "mmtpi"), w, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  #Checking feasibility conditions for pt, e1 and e2
  if(isTRUE(pt > 1 || pt < 0))
  {
    stop("Target toxicity Probability should take values between 0 and 1")
  }
  if(isTRUE(pt - e1 < 0 || pt + e2 > 1))
  {
    stop ("e1 and e2, two thresholds should be small compared to the target probability pt")
  }

  #checking feasibility condition for x and n
  if(isTRUE(n < 0))
  {
    stop ("n, total number of patients in a cohort must be positive")
  }
  if (isTRUE(x > n))
  {
    stop("Number of patients experiencing DLT 's must be less than or equal to total number of patients treated in that cohort.")
  }
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
  #Checking over toxicity of Doses
  threshold = 0
  if (w %in% c(0,1))
  {
    threshold = pbeta(pt, a1 + x, b1 + n - x, lower.tail = FALSE)
    if(isTRUE(threshold >= eta)) {return("DU")}
    else
    {
      a1 = a1 + x
      b1 = b1 + n - x
    }
  }
  else
  {
    params = weights_formulate(w = w, x = x, n = n, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
    w = params$weight
    a1 = params$param_inform[1]
    b1 = params$param_inform[2]
    a2 = params$param_noninform[1]
    b2 = params$param_noninform[2]
    threshold = w * pbeta(pt, a1 + x, b1 + n - x, lower.tail = FALSE) +
                         (1 - w) * pbeta(pt, a2 + x, b2 + n - x, lower.tail = FALSE)
    if(isTRUE(threshold >= eta)){return("DU")}
  }


  #breaking up the domain in compatible ranges according to Professor Yuan Ji 's paper
  if (design %in% c("tpi", "mtpi"))
  {
    interval = c(0, pt - e1, pt + e2 , 1)
    length_interval = length(interval)
    if(isTRUE(design == "mtpi"))
    {
      #Developing mTPI based Dose calculations
      upm_array = rep(0, length_interval - 1)
      for (i  in 1:(length_interval - 1))
      {
        upm_array[i] = UPM(w = w, a = interval[i], b = interval[i + 1], a1 = a1, b1 = b1, a2 = a2, b2 = b2)
      }
      max_upm = max(upm_array)
      location = which.max(upm_array)
    }
    else
    {
      #Developing TPI based Dose Calculations
      prob_array = rep(0, length_interval - 1)
      for (i in 1:(length_interval - 1))
      {
        prob_array[i] = (interval[i + 1] - interval[i]) * UPM(w = w, a = interval[i], b = interval[i + 1], a1 = a1, b1 = b1, a2 = a2, b2 = b2)
      }
      location = which.max(prob_array)
    }
    if(location == 1)
    {
      return("E")
    }
    else if (location == 2)
    {
      return("S")
    }
    else
    {
      return("D")
    }
  }
  else
  {
    breaks_lower = floor((pt - e1) / 0.1)
    breaks_upper = floor((1 - pt - e2)  / 0.1)
    interval = c(0, pt - e1 - 0.1 * (breaks_lower : 0) , pt + e2 + 0.1 * (0 : breaks_upper) , 1)
    length_interval = length(interval)
    upm_array = rep(0, length_interval - 1)

    #Dose Finding Decision
    for (i  in 1:(length_interval - 1))
    {
      upm_array[i] = UPM(w = w, a = interval[i], b = interval[i + 1], a1 = a1, a2 = a2, b1 = b1, b2 = b2)
    }
    max_upm = max(upm_array)
    location = which.max(upm_array)
    if(isTRUE(interval[location] < (pt - e1)))
    {
      return("E")
    }
    else if (isTRUE(interval[location] == (pt - e1)))
    {
      return("S")
    }
    else
    {
      return("D")
    }
  }
}
