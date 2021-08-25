#' @title Subsidiary Function for Comparative Causal Mediation Analysis
#' @description Subsidiary function to determine appropriate output and calculate confidence intervals
#' @usage decideOutput(para.df,noInteraction = TRUE,
#'   sigLevel = 0.05,boots = 1000)
#' @param para.df a data frame containing the final data to be analyzed, generated using \code{checkData()}.
#' @param noInteraction logical. If \code{TRUE} (the default), the assumption of no interaction between the treatments and mediator is employed in the analysis.
#' @param sigLevel significance level to use in construction of confidence intervals. Default is 0.05 (i.e. 95 percent confidence intervals).
#' @param boots number of bootstrap resamples taken for construction of confidence intervals.
#' @importFrom stats quantile
#' @return A list containing confidence intervals along with flags indicating the appropriate output to return to getCCM().
#' @note This function is called internally and thus should not be used directly.
#' @author Kirk Bansak and Xiaohan Wu
#' @references Bansak, K. (2020). Comparative causal mediation and relaxing the assumption of no mediator-outcome confounding: An application to international law and audience costs. Political Analysis, 28(2), 222-243.

decideOutput <- function(para.df,noInteraction = TRUE,
                         sigLevel = 0.05,boots = 1000){

  if (noInteraction == TRUE){
    cat("Estimation assuming no interactions between treatments and mediator.\n")
    ACME1.vec <- rep(NA,boots)
    ACME2.vec <- rep(NA,boots)
    ATE1.vec <- rep(NA,boots)
    ATE2.vec <- rep(NA,boots)
    PM1.vec <- rep(NA,boots)
    PM2.vec <- rep(NA,boots)
    ACMEratio.vec <- rep(NA,boots)
    PMratio.vec <- rep(NA,boots)
    n <- length(para.df$Y)
    for (j in 1:boots){

      k <- sample(seq(1:n),n,replace = T)
      ssdat <- para.df[k,]

      mod1 <- lm(M ~ T1 + T2, data=ssdat)
      mod2 <- lm(Y ~ T1 + T2 + M, data=ssdat)
      fmod <- lm(Y ~ T1 + T2, data=ssdat)

      alpha1 <- mod1$coef["T1"]
      alpha2 <- mod1$coef["T2"]
      beta <- mod2$coef["M"]
      ATE1 <- ATE1.vec[j] <- fmod$coef["T1"]
      ATE2 <- ATE2.vec[j] <- fmod$coef["T2"]
      ACME1.vec[j] <- alpha1*(beta)
      ACME2.vec[j] <- alpha2*(beta)
      ACMEratio.vec[j] <- (alpha2*beta)/(alpha1*beta)
      PM1.vec[j] <- (alpha1*(beta))/ATE1
      PM2.vec[j] <- (alpha2*(beta))/ATE2
      PMratio.vec[j] <- ((alpha2*(beta))/ATE2)/((alpha1*(beta))/ATE1)

      rm(k,ssdat,mod1,mod2,fmod,alpha1,alpha2,beta,ATE1,ATE2)

    }

    # scope conditions CIs do not contain 0.
    scope.condition1 <- quantile(ACME1.vec,probs = sigLevel/2)*quantile(ACME1.vec,probs = 1 - sigLevel/2)
    scope.condition2 <- quantile(ACME2.vec,probs = sigLevel/2)*quantile(ACME2.vec,probs = 1 - sigLevel/2)
    scope.condition3 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM1.vec,probs = 1 - sigLevel/2)
    scope.condition4 <- quantile(PM2.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = 1 - sigLevel/2)
    scope.condition5 <- quantile(ACME1.vec,probs = sigLevel/2)*quantile(ACME2.vec,probs = sigLevel/2)
    scope.condition6 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = sigLevel/2)
    if ((scope.condition1 <= 0)|(scope.condition2 <= 0)){
      flag1 <- 1
      flag2 <- 1
      warning("Scope conditions for estimand 1 and estimand 2 are not satisfied. Estimand 1 and estimand 2 cannot be reliably estimated.")
    }else{
      if (scope.condition5 <= 0){
        flag1 <- 1
        flag2 <- 1
        warning("ACME1 and ACME2 have different signs. Estimand 1 and estimand 2 cannot be reliably estimated.")
      }else{
        flag1 <- 0
        if ((scope.condition3 < 0)|(scope.condition4 < 0)){
          flag2 <- 1
          warning("Scope conditions for estimand 2 are not satisfied. Estimand 2 cannot be reliably estimated.")
        }else{
          if (scope.condition6 <= 0){
            flag2 <- 1
            warning("Proportions mediated have different signs. Estimand 2 cannot be reliably estimated.")
          }else{
            flag2 <- 0
          }
        }
      }
    }


    ATE1.ci <- c(quantile(ATE1.vec,probs = sigLevel/2),quantile(ATE1.vec,probs = 1 - sigLevel/2))
    ATE2.ci <- c(quantile(ATE2.vec,probs = sigLevel/2),quantile(ATE2.vec,probs = 1 - sigLevel/2))
    ACME1.ci <- c(quantile(ACME1.vec,probs = sigLevel/2),quantile(ACME1.vec,probs = 1 - sigLevel/2))
    ACME2.ci <- c(quantile(ACME2.vec,probs = sigLevel/2),quantile(ACME2.vec,probs = 1 - sigLevel/2))

    if (flag1 == 0){
      estimand1.ci <- c(quantile(ACMEratio.vec,probs = sigLevel/2),quantile(ACMEratio.vec,probs = 1 - sigLevel/2))
    }else{
      estimand1.ci <- c(NA,NA)
    }

    if (flag2 == 0){
      estimand2.ci <- c(quantile(PMratio.vec,probs = sigLevel/2),quantile(PMratio.vec,probs = 1 - sigLevel/2))
    }else{
      estimand2.ci <- c(NA,NA)
    }

    df_ci <- data.frame(ATE1.ci,ATE2.ci,ACME1.ci,ACME2.ci,estimand1.ci,estimand2.ci)
    out <- list(flag1,flag2,df_ci)

  }else{

    print("Estimation allowing for interactions between treatments and mediator.")
    ACMET1.vec <- rep(NA,boots)
    ACMET2.vec <- rep(NA,boots)
    ACMETratio.vec <- rep(NA,boots)
    PM1.vec <- rep(NA,boots)
    PM2.vec <- rep(NA,boots)
    PMratio.vec <- rep(NA,boots)

    alpha1.vec <- rep(NA,boots)
    alpha2.vec <- rep(NA,boots)
    omega1.vec <- rep(NA,boots)
    omega2.vec <- rep(NA,boots)
    tau1.vec <- rep(NA,boots)
    tau2.vec <- rep(NA,boots)

    n <- length(para.df$Y)

    for (j in 1:boots){

      k <- sample(seq(1:n),n,replace = T)
      ssdat <- para.df[k,]

      mod1 <- lm(M ~ T1 + T2, data=ssdat)
      mod2 <- lm(Y ~ T1 + T2 + M + T1*M + T2*M, data=ssdat)
      fmod <- lm(Y ~ T1 + T2, data=ssdat)

      alpha1 <- alpha1.vec[j] <- mod1$coef["T1"]
      alpha2 <- alpha2.vec[j] <- mod1$coef["T2"]
      beta <- mod2$coef["M"]
      gamma1 <- mod2$coef["T1:M"]
      gamma2 <- mod2$coef["T2:M"]
      omega1.vec[j] <- beta + gamma1
      omega2.vec[j] <- beta + gamma2
      ATE1 <- tau1.vec[j] <- fmod$coef["T1"]
      ATE2 <- tau2.vec[j] <- fmod$coef["T2"]

      ACMET1.vec[j] <- alpha1*(beta+gamma1)
      ACMET2.vec[j] <- alpha2*(beta+gamma2)
      PM1.vec[j] <- alpha1*(beta+gamma1)/ATE1
      PM2.vec[j] <- alpha2*(beta+gamma2)/ATE2

      ACMETratio.vec[j] <- (alpha2*(beta+gamma2))/(alpha1*(beta+gamma1))
      PMratio.vec[j] <- (alpha2*(beta+gamma2)/ATE2)/(alpha1*(beta+gamma1)/ATE1)

      rm(k,ssdat,mod1,mod2,fmod,alpha1,alpha2,beta,gamma1,gamma2,ATE1,ATE2)

    }

    scope.condition1 <- quantile(ACMET1.vec,probs = sigLevel/2)*quantile(ACMET1.vec,probs = 1 - sigLevel/2)
    scope.condition2 <- quantile(ACMET2.vec,probs = sigLevel/2)*quantile(ACMET2.vec,probs = 1 - sigLevel/2)
    scope.condition3 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM1.vec,probs = 1 - sigLevel/2)
    scope.condition4 <- quantile(PM2.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = 1 - sigLevel/2)
    scope.condition5 <- quantile(ACMET1.vec,probs = sigLevel/2)*quantile(ACMET2.vec,probs = sigLevel/2)
    scope.condition6 <- quantile(PM1.vec,probs = sigLevel/2)*quantile(PM2.vec,probs = sigLevel/2)
    if ((scope.condition1 <= 0)|(scope.condition2 <= 0)){
      flag1 <- 1
      flag2 <- 1
      warning("Scope conditions for estimand 1 and estimand 2 are not satisfied. Estimand 1 and estimand 2 cannot be reliably estimated.")
    }else{
      if (scope.condition5 <= 0){
        flag1 <- 1
        flag2 <- 1
        warning("ACME1 and ACME2 have different signs. Estimand 1 and estimand 2 cannot be reliably estimated.")
      }else{
        flag1 <- 0
        if ((scope.condition3 < 0)|(scope.condition4 < 0)){
          flag2 <- 1
          warning("Scope conditions for estimand 2 are not satisfied. Estimand 2 cannot be reliably estimated.")
        }else{
          if (scope.condition6 <= 0){
            flag2 <- 1
            warning("Proportions mediated have different signs. Estimand 2 cannot be reliably estimated.")
          }else{
            flag2 <- 0
          }
        }
      }
    }

    ATE1.ci <- c(quantile(tau1.vec,probs = sigLevel/2),quantile(tau1.vec,probs = 1 - sigLevel/2))
    ATE2.ci <- c(quantile(tau2.vec,probs = sigLevel/2),quantile(tau2.vec,probs = 1 - sigLevel/2))
    ACMET1.ci <- c(quantile(ACMET1.vec,probs = sigLevel/2),quantile(ACMET1.vec,probs = 1 - sigLevel/2))
    ACMET2.ci <- c(quantile(ACMET2.vec,probs = sigLevel/2),quantile(ACMET2.vec,probs = 1 - sigLevel/2))

    if (flag1 == 0){
      estimand1.ci <- c(quantile(ACMETratio.vec,probs = sigLevel/2),quantile(ACMETratio.vec,probs = 1 - sigLevel/2))
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
