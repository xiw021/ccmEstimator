context("test getCCM function")

test_that("whether the getCCM function gives correct output using ICAapp data",{
  set.seed(321, kind = "Mersenne-Twister", normal.kind = "Inversion")
  data(ICAapp)
  out = getCCM('dapprp','trt1','trt2','immorp',ICAapp)
  expect_equal(round(out[["ATE1"]],3),0.195)
  expect_equal(round(out[["ATE2"]],3),0.320)
  expect_equal(round(out[["ACME1"]],3),0.113)
  expect_equal(round(out[["ACME2"]],3),0.177)
  expect_equal(round(out[["estimand1"]],3),1.563)
  expect_equal(round(out[["estimand2"]],3),0.952)
  expect_equal(round(out[["confidenceIntervals"]],3),structure(list(ATE1.ci = c(0.144, 0.251),
                                                                    ATE2.ci = c(0.266,0.379),
                                                                    ACME1.ci = c(0.076, 0.15),
                                                                    ACME2.ci = c(0.138, 0.219),
                                                                    estimand1.ci = c(1.183, 2.16),
                                                                    estimand2.ci = c(0.755, 1.236)),
                                                               row.names = c("2.5%", "97.5%"), class = "data.frame"))
  set.seed(321, kind = "Mersenne-Twister", normal.kind = "Inversion")
  out = getCCM('nvotep','trt1','trt2','immorp',ICAapp,noInteraction = FALSE)
  expect_equal(round(out[["ATE1"]],3),0.182)
  expect_equal(round(out[["ATE2"]],3),0.281)
  expect_equal(round(out[["ACME1"]],3),0.096)
  expect_equal(round(out[["ACME2"]],3),0.176)
  expect_equal(round(out[["estimand1"]],3),1.829)
  expect_equal(round(out[["estimand2"]],3),1.184)
  expect_equal(round(out[["confidenceIntervals"]],3),structure(list(ATE1.ci = c(0.126, 0.235), ATE2.ci = c(0.225,0.337),
                                                                    ACMET1.ci = c(0.062, 0.131),
                                                                    ACMET2.ci = c(0.135, 0.219),
                                                                    estimand1.ci = c(1.287, 2.581),
                                                                    estimand2.ci = c(0.883, 1.566)),
                                                               row.names = c("2.5%", "97.5%"), class = "data.frame"))


  expect_warning(getCCM('nvotep','trt1','trt2','immorp',ICAapp), "Interactions between treatments and mediator may exist. Please consider setting noInteration as FALSE.", fixed = TRUE)

})

test_that("Whether the getCCM function give correct output under different scope conditions",{
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  T1 <- c(rep(0,50),rep(1,50),rep(0,50))
  T2 <- c(rep(0,50),rep(0,50),rep(1,50))
  dat <- data.frame(T1,T2)
  u1 <- rnorm(n=150,mean=0,sd = 2)
  u3 <- rnorm(n=150,mean=0,sd = 1)
  dat$M <- 1 + 0.1*dat$T1 + 6*dat$T2 + u1
  dat$Y <- 3 + 0.22*dat$T1 + 2*dat$T2 + 0.6*dat$M + u3
  expect_warning(getCCM('Y','T1','T2','M',data = dat),'Scope conditions for estimand 1 and estimand 2 are not satisfied. Estimand 1 and estimand 2 cannot be reliably estimated.',fixed = TRUE)
  expect_warning(getCCM('Y','T1','T2','M',data = dat,noInteraction = FALSE),'Scope conditions for estimand 1 and estimand 2 are not satisfied. Estimand 1 and estimand 2 cannot be reliably estimated.',fixed = TRUE)

  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  T1 <- c(rep(0,50),rep(1,50),rep(0,50))
  T2 <- c(rep(0,50),rep(0,50),rep(1,50))
  dat <- data.frame(T1,T2)
  u1 <- rnorm(n=150,mean=0,sd = 1)
  u3 <- rnorm(n=150,mean=0,sd = 2)
  dat$M <- 1 - 2*dat$T1 + 6*dat$T2 + u1
  dat$Y <- 3 + 1*dat$T1 + 2*dat$T2 + 0.6*dat$M + u3
  expect_warning(getCCM('Y','T1','T2','M',data = dat),'ACME1 and ACME2 have different signs. Estimand 1 and estimand 2 cannot be reliably estimated.',fixed = TRUE)
  expect_warning(getCCM('Y','T1','T2','M',data = dat,noInteraction = FALSE),'ACME1 and ACME2 have different signs. Estimand 1 and estimand 2 cannot be reliably estimated.',fixed = TRUE)


  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  T1 <- c(rep(0,50),rep(1,50),rep(0,50))
  T2 <- c(rep(0,50),rep(0,50),rep(1,50))
  dat <- data.frame(T1,T2)
  u1 <- rnorm(n=150,mean=0,sd = 1)
  u3 <- rnorm(n=150,mean=0,sd = 2)
  dat$M <- 1 + 2*dat$T1 + 6*dat$T2 + u1
  dat$Y <- 3 - 0.2*dat$T1 + 2*dat$T2 + 0.6*dat$M + u3
  expect_warning(getCCM('Y','T1','T2','M',data = dat),'Scope conditions for estimand 2 are not satisfied. Estimand 2 cannot be reliably estimated.',fixed = TRUE)
  expect_warning(getCCM('Y','T1','T2','M',data = dat,noInteraction = FALSE),'Scope conditions for estimand 2 are not satisfied. Estimand 2 cannot be reliably estimated.',fixed = TRUE)


  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
  T1 <- c(rep(0,50),rep(1,50),rep(0,50))
  T2 <- c(rep(0,50),rep(0,50),rep(1,50))
  dat <- data.frame(T1,T2)
  u1 <- rnorm(n=150,mean=0,sd = 1)
  u3 <- rnorm(n=150,mean=0,sd = 2)
  dat$M <- 1 + 2*dat$T1 + 6*dat$T2 + u1
  dat$Y <- 3 - 5*dat$T1 + 2*dat$T2 + 0.6*dat$M + u3
  expect_warning(getCCM('Y','T1','T2','M',data = dat),'Proportions mediated have different signs. Estimand 2 cannot be reliably estimated.',fixed = TRUE)
  expect_warning(getCCM('Y','T1','T2','M',data = dat,noInteraction = FALSE),'Proportions mediated have different signs. Estimand 2 cannot be reliably estimated.',fixed = TRUE)
})









