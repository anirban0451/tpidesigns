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

  if (w %in% c(0,1) && isTRUE(pbeta(pt, a1 + x, b1 + n -x, lower.tail = FALSE) > eta))
  {
    return("DU")
  }
  else
  {
    params = weights_formulate(w = w, x = x, n = n, a1 = a1, a2 = a2, b1 = b1, b2 = b2)
    w = params$weight
    a1 = params$param_inform[1]
    b1 = params$param_inform[2]
    a2 = params$param_noninform[1]
    b2 = params$param_noninform[2]
    if(isTRUE(w * pbeta(pt, a1 + x, b1 + n - x, lower.tail = FALSE) +
              (1 - w) * pbeta(pt, a2 + x, b2 + n - x, lower.tail = FALSE)> eta))
    {
      return("DU")
    }
  }

  #breaking up the domain in compatible ranges according to Professor Yuan Ji 's paper
  if (design %in% c("tpi", "mtpi"))
  {
    interval = c(0, pt - e1, pt + e2)
  }
  else
  {
    breaks_lower = floor((pt - e1) / 0.1)
  }

}
