#' Title
#'
#' @param nmax A
#' @param design B
#' @param pt C
#' @param e1 D
#' @param e2 E
#' @param eta F
#' @param w G
#' @param a1 H
#' @param b1 T
#' @param a2 G
#' @param b2 G
#'
#' @return g
#' @export
#'
#'
tpitable = function(nmax, design, pt, e1 = 0.05, e2 = 0.05, eta, w, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
{
  if(isTRUE(nmax < 3))
  {
    if(isTRUE(nmax <= 0))
    {
      stop("Number of patients must be positive and at least 3")
    }
    else
    {
      stop("Number of patients must be atleast 3")
    }
  }
  table = matrix(0, nrow = 4, ncol = nmax - 2)
  for(i in 3:nmax)
  {
    decisions = array(dim = 1)
    for(dlt in 0:i)
    {
      decisions[dlt + 1] = decisiontpi(pt, e1, e2, x = dlt, n = i, eta, design, w, a1, b1, a2, b2)
      if(isTRUE(decisions[dlt + 1] == "DU"))
      {
        table[4 , i - 2] = dlt
        break
      }
    }

    table[1 , i - 2] = ifelse(any(decisions == "E"), max(which(decisions == "E")) - 1, NA)

    table[2 , i - 2] = ifelse(any(decisions == "S"), max(which(decisions == "S")) - 1, NA)

    table[3 , i - 2] = ifelse(any(decisions == "D"), max(which(decisions == "D")) - 1, NA)
  }
  rownames_table = c("Maximum no. of DLT's allowed to Escalate from the current Dosage",
                     "Maximum no. of DLT's allowed to Stay at the current Dosage",
                     "Maximum no. of DLT's allowed to De-escalate to a lower Dosage",
                     "Minimum number of DLT's to declare current dosage as Unacceptably toxic")
  data.frame(table, row.names = rownames_table) -> tableframe
  names(tableframe) <- paste("n = ", 3:nmax)
  return(tableframe)
}
