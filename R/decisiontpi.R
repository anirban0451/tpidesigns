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

    #Calculation of Decision for w = 1

    if(pbeta(pt, a1 + x, b1 + n - x, lower.tail = FALSE) > eta)
    {
      return("DU")
    }
  }

}
