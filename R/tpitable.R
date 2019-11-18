#' Tabular display of Dose Escalation Decisions
#'
#' \code{tpitable} gives a table of Dose Escalation Decisions when a certain number of DLT 's  are observed
#' in different sample size
#' @param nmax Maximum number of patients to be treated in current level of Dose, must be at least 3
#' @inheritParams decisiontpi
#'
#' @return A table containing threshold number of DLT 's for Dose Escalation decisions in different cohort size
#' upto the maximum number of patients allowed for treatment
#' @details
#' In Oncology Trials, often the maximum number of people allowed in the treatment is fixed beforehand, and the Dose Escalation
#' Decision starts from including three patients upto maximum number of patients allowed in the cohort. Hence, the motive of this table
#' is to create a frame which will be useful when a Clinician working on these type of Trials experiences a certain number of DLT 's in
#' the sample, he can straightway look at the table and take Dose Escalation Decisions.It is important to note that,
#' the number of patients in a cohort must be atleast 3, addressed by Continual Reassessment Method
#' @seealso
#' \url{https://www.cancer.gov/publications/dictionaries/cancer-terms/def/dose-limiting} for the Definition of Dose Limiting Toxicity\cr
#' \code{\link{decisiontpi}} to know about Dose Escalation Decision Strategy
#' @export
#'
#' @examples
#' tpitable(nmax = 5, design = "mtpi", pt = 0.3, eta = 0.75, w = 1, a1 = 2, b1 = 3)
#' tpitable(nmax = 13, design = "mmtpi", pt = 0.4, e1 = 0.06, e2 = 0.04, eta = 0.95, w = 0.4, a1 = 4, b1 = 3, a2 = 1, b2 = 1)
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
