#' @title Comparative Causal Mediation Analysis
#' @description Function to perform comparative causal mediation analysis to compare the mediation effects of different treatments via a common mediator.
#' Function requires two separate treaments (as well as a control condition) and one mediator.
#' @usage getCCM(Y,T1,T2,M,data = NULL,
#'   noInteraction = TRUE,sigLevel = 0.05,
#'   boots = 1000)
#' @param Y numeric outcome variable. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param T1 binary indicator for first treatment. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param T2 binary indicator for second treatment. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param M numeric mediator variable. Should be a vector if a data frame is not provided through the \code{data} argument, or the ("character") name of the variable in the data frame if provided.
#' @param data an optional data frame containing the variables to be used in analysis
#' @param noInteraction logical. If \code{TRUE} (the default), the assumption of no interaction between the treatments and mediator is employed in the analysis.
#' @param sigLevel significance level to use in construction of confidence intervals. Default is 0.05 (i.e. 95 percent confidence intervals).
#' @param boots number of bootstrap resamples taken for construction of confidence intervals.
#' @details Function will automatically assess the comparative causal mediation analysis scope conditions
#' (i.e. for each comparative causal mediation estimand, a numerator and denominator that are both statistically significantly estimated and of the same sign).
#' Results will be returned for each comparative causal mediation estimand only if scope conditions are met for it.
#' See "Scope Conditions" section in Bansak (2020) for more information.
#' Results will also be returned for the ATE and ACME for each treatment.
#'
#' If \code{noInteraction = TRUE} (the default setting), function will automatically assess the possibility of interactions between treatments and mediator and return a warning in case evidence of such interactions are found.
#' @importFrom stats var
#' @return A ccmEstimation object, which contains the estimates and confidence intervals for the two comparative causal mediation analysis estimands, as well as the ATE and ACME for each treatment.
#' User should input the ccmEstimation object into the \code{summary()} function to view the estimation results.
#'
#' Note also that the comparative causal mediation analysis results and interpretation of the results
#' will be printed in the console.
#' @author Kirk Bansak and Xiaohan Wu
#' @references Bansak, K. (2020). Comparative causal mediation and relaxing the assumption of no mediator-outcome confounding: An application to international law and audience costs. Political Analysis, 28(2), 222-243.
#' @export
#'
#' @examples
#' #Example from application in Bansak (2020)
#' data(ICAapp)
#' set.seed(321, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' ccm.results <-
#'    getCCM(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp,
#'    noInteraction = TRUE, sigLevel = 0.05, boots = 10000)
#' summary(ccm.results)

getCCM <- function(Y,T1,T2,M,data = NULL,
                   noInteraction = TRUE,sigLevel = 0.05,
                   boots = 1000){
  para.df = checkData(Y,T1,T2,M,data)

  if (noInteraction == TRUE){
    interaction.check = checkInteraction(para.df,sigLevel)
    if (interaction.check == 1){
      warning('Interactions between treatments and mediator may exist. Please consider setting noInteration as FALSE.')
    }
  }

  decide.output <- decideOutput(para.df,noInteraction,sigLevel,boots)
  flag1 <- decide.output[[1]]
  flag2 <- decide.output[[2]]
  df.ci <- decide.output[[3]]

  if (noInteraction == TRUE){
    mod1 <- lm(M ~ T1 + T2, data=para.df)
    mod2 <- lm(Y ~ T1 + T2 + M, data=para.df)
    fmod <- lm(Y ~ T1 + T2, data=para.df)

    model.list = list(mod1,mod2,fmod)

    alpha1 <- mod1$coef["T1"]
    alpha2 <- mod1$coef["T2"]
    beta <- mod2$coef["M"]
    ATE1 <- fmod$coef["T1"]
    ATE2 <- fmod$coef["T2"]
    ACME1 <- alpha1*beta
    ACME2 <- alpha2*beta

    if (flag1 == 0){
      estimand1 <- ACME2/ACME1
    }else{
      estimand1 <- NA
    }

    if (flag2 == 0){
      PropMed1 <- alpha1*beta/ATE1
      PropMed2 <- alpha2*beta/ATE2
      estimand2 <- PropMed2/PropMed1
    }else{
      estimand2 <- NA
    }

    out <- list(models = model.list,
                ATE1 = unname(ATE1),
                ATE2 = unname(ATE2),
                ACME1 = unname(ACME1),
                ACME2 = unname(ACME2),
                estimand1 = unname(estimand1),
                estimand2 = unname(estimand2),
                confidenceIntervals = df.ci)

  }else{

    mod1 <- lm(M ~ T1 + T2, data=para.df)
    mod2 <- lm(Y ~ T1 + T2 + M + T1*M + T2*M, data=para.df)
    fmod <- lm(Y ~ T1 + T2, data=para.df)

    model.list = list(mod1,mod2,fmod)

    alpha1 <- mod1$coef["T1"]
    alpha2 <- mod1$coef["T2"]
    beta <- mod2$coef["M"]
    gamma1 <- mod2$coef["T1:M"]
    gamma2 <- mod2$coef["T2:M"]
    ATE1 <- tau1 <- fmod$coef["T1"]
    ATE2 <- tau2 <- fmod$coef["T2"]
    omega1 <- beta + gamma1
    omega2 <- beta + gamma2
    ACMET1 <- alpha1*(beta+gamma1)
    ACMET2 <- alpha2*(beta+gamma2)

    if (flag1 == 0){
      estimand1 <- ACMET2/ACMET1
    }else{
      estimand1 <- NA
    }

    if (flag2 == 0){
      PropMed1 <- alpha1*(beta+gamma1)/ATE1
      PropMed2 <- alpha2*(beta+gamma2)/ATE2
      estimand2 <- PropMed2/PropMed1
    }else{
      estimand2 <- NA
    }

    out <- list(models = model.list,
                ATE1 = unname(ATE1),
                ATE2 = unname(ATE2),
                ACME1 = unname(ACMET1),
                ACME2 = unname(ACMET2),
                estimand1 = unname(estimand1),
                estimand2 = unname(estimand2),
                confidenceIntervals = df.ci)

  }

  if (!is.na(estimand1)){
    if (estimand1 > 1){
      cat(paste0('The mediation effect for treatment 2 is ' , round(estimand1,3), ' times larger than the mediation effect for treatment 1 \n   (with ',(1-sigLevel)*100,'% CI: [', round(out$confidenceIntervals['2.5%','estimand1.ci'],3),',',round(out$confidenceIntervals['97.5%','estimand1.ci'],3),'])\n'))
    }else{
      cat(paste0('The mediation effect for treatment 1 is ' , round(1/estimand1,3), ' times larger than the mediation effect for treatment 2 \n   (with ',(1-sigLevel)*100,'% CI: [', round(1/out$confidenceIntervals['97.5%','estimand1.ci'],3),',',round(1/out$confidenceIntervals['2.5%','estimand1.ci'],3),'])\n'))
    }
  }

  if (!is.na(estimand2)){
    if (estimand2 > 1){
      cat(paste0('The proportion mediated for treatment 2 is ' , round(estimand2,3), ' times larger than that for treatment 1 \n   (with ',(1-sigLevel)*100,'% CI: [', round(out$confidenceIntervals['2.5%','estimand2.ci'],3),',',round(out$confidenceIntervals['97.5%','estimand2.ci'],3),'])\n'))
    }else{
      cat(paste0('The proportion mediated for treatment 1 is ' , round(1/estimand2,3), ' times larger than that for treatment 2 \n   (with ',(1-sigLevel)*100,'% CI: [', round(1/out$confidenceIntervals['97.5%','estimand2.ci'],3),',',round(1/out$confidenceIntervals['2.5%','estimand2.ci'],3),'])\n'))
    }
  }

  cat('Please use summary() to view the results in more detail.')
  out$call <- match.call()
  class(out) <- "ccmEstimation"
  out

}
