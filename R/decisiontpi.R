#' Dosing Decisions based on different tpi designs
#'
#' \code{decisiontpi} gives whether to escalate, stay or de-escalate to a level of Dose given total number of patients and total number of people experiencing toxicities in a cohort
#' @param pt Target toxicity proportion to achieve in the sample (lower toxicity implies underdosing)
#' @param e1 Allowable deviation towards the left of target toxicity \code{pt}, usually lies between 0 to 0.05
#' @param e2 Allowable deviation towards the right of target toxicity \code{pt}, usually lies between 0 to 0.05
#' @param x Number of patients experiencing DLT (Dose Limiting Toxicity) 's in the cohort
#' @param n Total number of patients in the cohort
#' @param eta threshold value to check if the Dose is severely toxic
#' @param design Design parameter, tells us which design to use. Options are \code{"tpi", "mtpi", "mmtpi"}
#' @param w weight on the informative prior
#' @param a1,b1 alpha and beta parameters for informative Beta Prior component
#' @param a2,b2 alpha and Beta parameters for noninformative Beta Prior component
#' @return
#' \code{"E"} We should increase the current level of Dose.\cr
#' \code{"S"} We should stay at the current level of Dose and treat more patients.\cr
#' \code{"D"} We should decrease the current level of Dose.\cr
#' \code{"DU"} We should decrease the current level of Dose and we should not go beyond this Dose level for treating more patients\cr
#' @export
#'
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
    params = weights_formulate(w = w, x = x, n = n, a1 = a1, a2 = a2, b1 = b1, b2 = b2)
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
    upm_array = rep(0, length_interval - 1)
    for (i  in 1:(length_interval - 1))
    {
      upm_array[i] = UPM(w = w, a1 = a1, a2 = a2, b1 = b1, b2 = b2, a = interval[i], b = interval[i + 1])
    }
    max_upm = max(upm_array)
    location = which.max(upm_array)
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
      upm_array[i] = UPM(w = w, a1 = a1, a2 = a2, b1 = b1, b2 = b2, a = interval[i], b = interval[i + 1])
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
