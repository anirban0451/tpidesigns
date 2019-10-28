#' Title
#'
#' @param pt
#' @param e1
#' @param e2
#' @param x
#' @param n
#' @param eta
#' @param design
#' @param w
#' @param a1
#' @param b1
#' @param a2
#' @param b2
#'
#' @return
#' @export
#'
#' @examples
decisiontpi <- function(pt, e1 = 0.05, e2 = 0.05, x, n, eta, design = c("tpi", "mtpi", "mmtpi"), w, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  #Checking feasibility conditions for pt, e1 and e2
  if(pt > 1 || pt < 0)
  {
    stop("Target toxicity Probability should take values between 0 and 1")
  }
  if(pt - e1 < 0 || pt + e2 > 1)
  {
    stop ("e1 and e2, two thresholds should be small compared to the target probability pt")
  }

  #checking feasibility condition for x and n
  if(n < 0)
  {
    stop ("n, total number of patients in a cohort must be positive")
  }
  if (x > n)
  {
    stop("Number of patients experiencing DLT 's must be less than or equal to total number of patients treated in that cohort.")
  }

  #Checking over toxicity of Doses
  threshold = 0
  if (w %in% c(0,1))
  {
    threshold = pbeta(pt, a1 + x, b1 + n -x, lower.tail = FALSE)
    if (threshold >= eta)return("DU")
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
    if(threshold >= eta)
    {
      return("DU")
    }
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
    if(interval[location] < (pt - e1))
    {
      return("E")
    }
    else if (interval[location] == (pt - e1))
    {
      return("S")
    }
    else
    {
      return("D")
    }
  }
}
