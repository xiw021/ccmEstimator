#' @export

summary.ccmEstimation <- function(object,...)
{
  point.estimate = c(object$ATE1,object$ATE2,object$ACME1,object$ACME2)
  lower.ci = unname(as.list(object$confidenceIntervals[1,1:4]))
  upper.ci = unname(as.list(object$confidenceIntervals[2,1:4]))
  tab1 = cbind(point.estimate,lower.ci,upper.ci)
  rownames(tab1) <- c("ATE1","ATE2","ACME1*","ACME2*")

  point.estimate = c(object$estimand1,object$estimand2)
  lower.ci = unname(as.list(object$confidenceIntervals[1,5:6]))
  upper.ci = unname(as.list(object$confidenceIntervals[2,5:6]))
  tab2 = cbind(point.estimate,lower.ci,upper.ci)
  rownames(tab2) <- c("ACME2/ACME1","(ACME2/ATE2)/(ACME1/ATE1)")
  res <- list(call = object$call,
              tab1 = tab1,
              tab2 = tab2,
              sigLevel = object$sigLevel)
  class(res) <- "summary.ccmEstimation"
  res
}


#'@export

print.summary.ccmEstimation <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$tab1)
  cat("\n")
  cat('The significance level is set to ')
  cat(x$sigLevel)
  cat("\n")
  cat("* Biased estimation. See Bansak (2020) for details.\n")
  cat("\n")
  print(x$tab2)
}
