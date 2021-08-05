
#' Check Interaction between mediator and treatments
#'
#' @param para.df A data.frame of arguments generated using checkArguments().
#' @param sigLevel The significance level, 0.05 by default.
#'
#' @return Logical. 1 if the interaction between treatments and mediator exists.
#' @importFrom stats lm
#' @export
#'
#' @examples
#' \dontrun{
#' df = checkArguments('dapprp','trt1','trt2','immorp',testdata)
#' interaction = checkInteration(df)
#' }
#'
checkInteraction <- function(para.df,sigLevel = 0.05){
  interactionMod <- lm(Y ~ T1 + T2 + M + T1*M + T2*M, data=para.df)
  pval.gamma1 <- summary(interactionMod)$coef["T1:M","Pr(>|t|)"]
  pval.gamma2 <- summary(interactionMod)$coef["T2:M","Pr(>|t|)"]
  if ((pval.gamma1 < 0.05)|(pval.gamma1 < 0.05)){
   interaction = 1
  }else{
    interaction = 0
  }
  return(interaction)
}
