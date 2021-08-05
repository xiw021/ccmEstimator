#' Summary of ccmEstimation
#'
#' @param object
#' @param ...
#'
#' @return Summary tables
#' @export
#' @examples
#' \dontrun{
#' df = checkArguments('dapprp','trt1','trt2','immorp',testdata)
#' summary(df)
#' }

summary.ccmEstimation <- function(object, ...)
{
  point.estimation = c(object$ATE1,object$ATE2,object$ACME1,object$ACME2)
  lower.ci = unname(as.list(object$confidenceIntervals[1,1:4]))
  upper.ci = unname(as.list(object$confidenceIntervals[2,1:4]))
  tab1 = cbind(point.estimation,lower.ci,upper.ci)
  rownames(tab1) <- c("ATE1","ATE2","ACME1*","ACME2*")

  point.estimation = c(object$estimand1,object$estimand2)
  lower.ci = unname(as.list(object$confidenceIntervals[1,5:6]))
  upper.ci = unname(as.list(object$confidenceIntervals[2,5:6]))
  tab2 = cbind(point.estimation,lower.ci,upper.ci)
  rownames(tab2) <- c("ACME2/ACME1","(ACME2/ATE2)/(ACME1/ATE1)")
  res <- list(call = object$call,
              tab1 = tab1,
              tab2 = tab2)
  class(res) <- "summary.ccmEstimation"
  res
}


#' Print summary tables
#'
#' @param x
#' @param ...
#'
#' @return Print the summary tables
#' @export
#'
#' @examples
#' #' \dontrun{
#' df = checkArguments('dapprp','trt1','trt2','immorp',testdata)
#' summary(df)
#' }
print.summary.ccmEstimation <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$tab1)
  cat("\n")
  cat("* Biased estimation. See Bansak(2019) for details.\n")
  cat("\n")
  print(x$tab2)
}
