#' A function that checks the interaction
#'
#' @param Y
#' @param T1
#' @param T2
#' @param M
#' @param data
#' @param sigLevel
#'
#' @return a logical value. 1 if interaction between mediator and treatment exists.
#' @export
#'
#' @examples
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
