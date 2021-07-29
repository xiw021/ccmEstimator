
#' Decide output and calculate confidence intervals
#'
#' @param para.df
#' @param noInteration
#' @param boots
#' @param sigLevel
#'
#' @return
#' @export
#'
#' @examples
decideOutput <- function(para.df,
                         noInteraction = TRUE,boots = 1000,sigLevel = 0.05){
  if (noInteraction == TRUE){
    print("Causal mediation estimates assuming no interactions between treatments and mediator.")
    ACME1hat.vec <- c()
    ACME2hat.vec <- c()
    ATE1hat.vec <- c()
    ATE2hat.vec <- c()
    PM1.vec <- c()
    PM2.vec <- c()
    ACMEratiohat.vec <- c()
    PMratio.vec <- c()
    n <- length(para.df$Y)
    for (j in 1:boots){
      if (j == 1){
        ssdat <- para.df
      }
      if (j > 1){
        k <- sample(seq(1:n),n,replace = T)
        ssdat <- para.df[k,]
      }
      mod1 <- lm(M ~ T1 + T2, data=ssdat)
      mod2 <- lm(Y ~ T1 + T2 + M, data=ssdat)
      fmod <- lm(Y ~ T1 + T2, data=ssdat)
      mmod <- lm(cbind(M,Y) ~ T1 + T2, data=ssdat)
      alpha1hat <- mod1$coef["T1"]
      alpha2hat <- mod1$coef["T2"]
      betahat <- mod2$coef["M"]
      ATE1hat <- ATE1hat.vec[j] <- fmod$coef["T1"]
      ATE2hat <- ATE2hat.vec[j] <- fmod$coef["T2"]
      ACME1hat.vec[j] <- alpha1hat*(betahat)
      ACME2hat.vec[j] <- alpha2hat*(betahat)
      ACMEratiohat.vec[j] <- (alpha2hat*betahat)/(alpha1hat*betahat)
      PM1.vec[j] <- (alpha1hat*(betahat))/ATE1hat
      PM2.vec[j] <- (alpha2hat*(betahat))/ATE2hat
      PMratio.vec[j] <- ((alpha2hat*(betahat))/ATE2hat)/((alpha1hat*(betahat))/ATE1hat)
    }

    # scope conditions CIs do not contain 0.
    scope.condition1 <- quantile(ACME1hat.vec,probs = sigLevel/2)*quantile(ACME1hat.vec,probs = 1 - sigLevel/2)
    scope.condition2 <- quantile(ACME2hat.vec,probs = sigLevel/2)*quantile(ACME2hat.vec,probs = 1 - sigLevel/2)
    scope.condition3 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM1.vec,probs = 1 - sigLevel/2)
    scope.condition4 <- quantile(PM2.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = 1 - sigLevel/2)
    if ((scope.condition1 <= 0)|(scope.condition2 <= 0)){
      flag1 <- 1
      flag2 <- 1
      warning("Scope conditions for estimand 1 and estimand 2 are not satisfied. Estimand1 and estimand2 are not available.")
    }else{
      flag1 <- 0
      if ((scope.condition3 <0)|(scope.condition4 <0)){
        flag2 <- 1
        warning("Scope conditions for  estimand 2 are not satisfied.")
      }else{
        flag2 <- 0
      }
    }

    scope.condition5 <- quantile(ACME1hat.vec,probs = sigLevel/2)*quantile(ACME2hat.vec,probs = sigLevel/2)
    scope.condition6 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = sigLevel/2)
    if (scope.condition5 <= 0){
      flag1 <- 1
      flag2 <- 1
      warning("ACME1 and ACME2 have different sign. Estimand1 and estimand2 are not available.")
    }else{
      if (scope.condition6 <= 0){
        flag2 <- 1
        warning("Proportions mediated have different sign. Estimand2 is not available.")
      }
    }

    ATE1.ci <- c(quantile(ATE1hat.vec,probs = sigLevel/2),quantile(ATE1hat.vec,probs = 1 - sigLevel/2))
    ATE2.ci <- c(quantile(ATE2hat.vec,probs = sigLevel/2),quantile(ATE2hat.vec,probs = 1 - sigLevel/2))
    ACME1.ci <- c(quantile(ACME1hat.vec,probs = sigLevel/2),quantile(ACME1hat.vec,probs = 1 - sigLevel/2))
    ACME2.ci <- c(quantile(ACME2hat.vec,probs = sigLevel/2),quantile(ACME2hat.vec,probs = 1 - sigLevel/2))
    if (flag1 == 0){
      estimand1.ci <- c(quantile(ACMEratiohat.vec,probs = sigLevel/2),quantile(ACMEratiohat.vec,probs = 1 - sigLevel/2))
    }else{
      estimand1.ci <- c(NA,NA)
    }

    if (flag2 == 0){
      estimand2.ci <- c(quantile(PMratio.vec,probs = sigLevel/2),quantile(PMratio.vec,probs = 1 - sigLevel/2))
    }else{
      estimand2.ci <- c(NA,NA)
    }
    #estimand1.ci <- c(quantile(ACMEratiohat.vec,probs = sigLevel/2),quantile(ACMEratiohat.vec,probs = 1 - sigLevel/2))
    #estimand2.ci <- c(quantile(PMratio.vec,probs = sigLevel/2),quantile(PMratio.vec,probs = 1 - sigLevel/2))
    df_ci <- data.frame(ATE1.ci,ATE2.ci,ACME1.ci,ACME2.ci,estimand1.ci,estimand2.ci)
    out <- list(flag1,flag2,df_ci)
  }else{
    print("Causal mediation estimates allowing for interactions between treatments and mediator.")
    ACMET1hat.vec <- c()
    ACMET2hat.vec <- c()
    ACMETratiohat.vec <- c()
    PM1.vec <- c()
    PM2.vec <- c()
    PMratio.vec <- c()

    alpha1.vec <- c()
    alpha2.vec <- c()
    omega1.vec <- c()
    omega2.vec <- c()
    tau1.vec <- c()
    tau2.vec <- c()

    n <- length(para.df$Y)

    for (j in 1:boots){

      if (j == 1){
        ssdat <- para.df
      }

      if (j > 1){
        k <- sample(seq(1:n),n,replace = T)
        ssdat <- para.df[k,]
      }

      mod1 <- lm(M ~ T1 + T2, data=ssdat)
      mod2 <- lm(Y ~ T1 + T2 + M + T1*M + T2*M, data=ssdat)
      fmod <- lm(Y ~ T1 + T2, data=ssdat)
      mmod <- lm(cbind(M,Y) ~ T1 + T2, data=ssdat)

      alpha1 <- alpha1.vec[j] <- mod1$coef["T1"]
      alpha2 <- alpha2.vec[j] <- mod1$coef["T2"]
      beta <- mod2$coef["M"]
      gamma1 <- mod2$coef["T1:M"]
      gamma2 <- mod2$coef["T2:M"]
      omega1.vec[j] <- beta + gamma1
      omega2.vec[j] <- beta + gamma2
      ATE1 <- tau1.vec[j] <- fmod$coef["T1"]
      ATE2 <- tau2.vec[j] <- fmod$coef["T2"]

      ACMET1hat.vec[j] <- alpha1*(beta+gamma1)
      ACMET2hat.vec[j] <- alpha2*(beta+gamma2)
      PM1.vec[j] <- alpha1*(beta+gamma1)/ATE1
      PM2.vec[j] <- alpha2*(beta+gamma2)/ATE2

      ACMETratiohat.vec[j] <- (alpha2*(beta+gamma2))/(alpha1*(beta+gamma1))
      PMratio.vec[j] <- (alpha2*(beta+gamma2)/ATE2)/(alpha1*(beta+gamma1)/ATE1)

    }

    scope.condition1 <- quantile(ACMET1hat.vec,probs = sigLevel/2)*quantile(ACMET1hat.vec,probs = 1 - sigLevel/2)
    scope.condition2 <- quantile(ACMET2hat.vec,probs = sigLevel/2)*quantile(ACMET2hat.vec,probs = 1 - sigLevel/2)
    scope.condition3 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM1.vec,probs = 1 - sigLevel/2)
    scope.condition4 <- quantile(PM2.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = 1 - sigLevel/2)
    if ((scope.condition1 <= 0)|(scope.condition2 <= 0)){
      flag1 <- 1
      flag2 <- 1
      warning("Scope conditions for estimand 1 and estimand 2 are not satisfied.Estimand1 and estimand2 are not available.")
    }else{
      flag1 <- 0
      if ((scope.condition3 <= 0)|(scope.condition4 <= 0)){
        flag2 <- 1
        warning("Scope conditions for  estimand 2 are not satisfied. Estimand2 is not available.")
      }else{
        flag2 <- 0
      }
    }

    scope.condition5 <- quantile(ACMET1hat.vec,probs = sigLevel/2)*quantile(ACMET2hat.vec,probs = sigLevel/2)
    scope.condition6 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = sigLevel/2)
    if (scope.condition5 <= 0){
      flag1 <- 1
      flag2 <- 1
      warning("ACMET1 and ACMET2 have different sign. Estimand1 and estimand2 are not available.")
    }else{
      if (scope.condition6 <= 0){
        flag2 <- 1
        warning("Proportions mediated have different sign. Estimand2 is not available.")
      }
    }

    ATE1.ci <- c(quantile(tau1.vec,probs = sigLevel/2),quantile(tau1.vec,probs = 1 - sigLevel/2))
    ATE2.ci <- c(quantile(tau2.vec,probs = sigLevel/2),quantile(tau2.vec,probs = 1 - sigLevel/2))
    ACMET1.ci <- c(quantile(ACMET1hat.vec,probs = sigLevel/2),quantile(ACMET1hat.vec,probs = 1 - sigLevel/2))
    ACMET2.ci <- c(quantile(ACMET2hat.vec,probs = sigLevel/2),quantile(ACMET2hat.vec,probs = 1 - sigLevel/2))
    if (flag1 == 0){
      estimand1.ci <- c(quantile(ACMETratiohat.vec,probs = sigLevel/2),quantile(ACMETratiohat.vec,probs = 1 - sigLevel/2))

    }else{
      estimand1.ci <- c(NA,NA)
    }

    if (flag2 == 0){
      estimand2.ci <- c(quantile(PMratio.vec,probs = sigLevel/2),quantile(PMratio.vec,probs = 1 - sigLevel/2))
    }else{
      estimand2.ci <- c(NA,NA)
    }

    df_ci <- data.frame(ATE1.ci,ATE2.ci,ACMET1.ci,ACMET2.ci,estimand1.ci,estimand2.ci)
    out <- list(flag1,flag2,df_ci)
  }
  return(out)
}
