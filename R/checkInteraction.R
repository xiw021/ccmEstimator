#' @title Subsidiary Function for Comparative Causal Mediation Analysis
#' @description Subsidiary function to assess interactions between treatments and mediator
#' @usage checkInteraction(para.df,sigLevel = 0.05)
#' @param para.df a data frame containing the final data to be analyzed, generated using \code{checkData()}.
#' @param sigLevel significance level used to conduct hypothesis tests.
#' @return logical. 1 if there is evidence of interaction between treatment(s) and mediator.
#' @note This function is called internally and should not be used directly.
#' @author Kirk Bansak and Xiaohan Wu
#' @references Bansak, K. (2020). Comparative causal mediation and relaxing the assumption of no mediator-outcome confounding: An application to international law and audience costs. Political Analysis, 28(2), 222-243.
#' @importFrom stats lm
#' @export
#' @examples
#' data(ICAapp)
#' final.dat <- checkData(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp)
#' checkInteraction(final.dat,sigLevel = 0.05)
#'
checkInteraction <- function(para.df,sigLevel = 0.05){
  interactionMod <- lm(Y ~ T1 + T2 + M + T1*M + T2*M, data=para.df)
  pval.gamma1 <- summary(interactionMod)$coef["T1:M","Pr(>|t|)"]
  pval.gamma2 <- summary(interactionMod)$coef["T2:M","Pr(>|t|)"]
  if ((pval.gamma1 < sigLevel)|(pval.gamma1 < sigLevel)){
    interaction = 1
  }else{
    interaction = 0
  }
  return(interaction)
}
