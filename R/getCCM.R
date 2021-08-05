#' Get CCM estimands
#'
#' @param Y A numeric vector of outcome variable, or the name of the outcome variable (as character/string) in data.frame if data provided.
#' @param T1 A binary vector of the first treatment. 1 if receiving the first treatment. Or the name of the first treatment variable (as character/string) in data.frame if data provided.
#' @param T2 A binary vector of the second treatment. 1 if receiving the second treatment. Or the name of the second treatment variable (as character/string) in data.frame if data provided.
#' @param M A vector of an intermediary variable that is affected by T1 and T2 and that affects Y. Or the name of the intermediary variable (as character/string) in data.frame if data provided.
#' @param data  An optional data frame that contains the variables (Y, T1, T2, M)  in the model.
#' @param noInteraction Logical, TRUE by default.  If TRUE, there is no expected interaction between the treatments and mediator.
#' @param sigLevel The significance level, 0.05 by default.
#' @param boots Number of bootstrap iterations, must be larger than 1000, 1000 by default.
#' @param finiteSample Logical, FALSE by default. If False, non finite sample adjustment.
#' @importFrom stats var
#' @return A ccmEstimation object.
#' @export
#'
#' @examples
#' \dontrun{
#' ccm.model = getCCM('dapprp','trt1','trt2','immorp',testdata)
#' }

getCCM <- function(Y,T1,T2,M,data = NULL,
                   noInteraction = TRUE,sigLevel = 0.05,
                   boots = 1000, finiteSample = FALSE){
  para.df = checkArguments(Y,T1,T2,M,data)

  if (noInteraction == TRUE){
    interaction.check = checkInteraction(para.df,sigLevel)
    if (interaction.check == 1){
      warning('Interactions between treatment and mediator exist.  Please consider setting noInteration as FALSE.')
    }
  }

  decide.output <- decideOutput(para.df,noInteraction,boots,sigLevel)
  flag1 <- decide.output[[1]]
  flag2 <- decide.output[[2]]
  df.ci <- decide.output[[3]]

  if (noInteraction == TRUE){
    mod1 <- lm(M ~ T1 + T2, data=para.df)
    mod2 <- lm(Y ~ T1 + T2 + M, data=para.df)
    fmod <- lm(Y ~ T1 + T2, data=para.df)
    mmod <- lm(cbind(M,Y) ~ T1 + T2, data=para.df)

    model.list = list(mod1,mod2,fmod,mmod)

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
      cat(paste0('The mediation effect for treatment 2 is ' , round(estimand1,3), ' times larger than the mediation effect for treatment 1 (with ',(1-sigLevel)*100,'% CI: [', round(out$confidenceIntervals['2.5%','estimand1.ci'],3),',',round(out$confidenceIntervals['97.5%','estimand1.ci'],3),'])\n'))
    }else{
      cat(paste0('The mediation effect for treatment 1 is' , round(1/estimand1,3), ' times larger than the mediation effect for treatment 2 (with ',(1-sigLevel)*100,'% CI: [', round(1/out$confidenceIntervals['97.5%','estimand1.ci'],3),',',round(1/out$confidenceIntervals['2.5%','estimand1.ci'],3),'])\n'))
    }
  }

  if (!is.na(estimand2)){
    if (estimand2 > 1){
      cat(paste0('The proportion mediated for treatment 2 is ' , round(estimand2,3), ' times larger than that for treatment 1 (with ',(1-sigLevel)*100,'% CI: [', round(out$confidenceIntervals['2.5%','estimand2.ci'],3),',',round(out$confidenceIntervals['97.5%','estimand2.ci'],3),'])\n'))
    }else{
      cat(paste0('The proportion mediated for treatment 1 is ' , round(1/estimand2,3), ' times larger than that for treatment 2 (with ',(1-sigLevel)*100,'% CI: [', round(1/out$confidenceIntervals['97.5%','estimand2.ci'],3),',',round(1/out$confidenceIntervals['2.5%','estimand2.ci'],3),'])\n'))
    }
  }
  cat('Please use summary() to view the estimation results.')
  out$call <- match.call()
  class(out) <- "ccmEstimation"
  out
}

