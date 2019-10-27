decisiontpi <- function(pt, e1 = 0.05, e2 = 0.05, x, n, design = c("tpi", "mtpi", "mmtpi"), w, a1, b1, a2, b2)
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

}
