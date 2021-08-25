context("test getCCM function")

test_that("whether the getCCM function gives correct output using ICAapp data",{
  set.seed(1)
  data(ICAapp)
  out = getCCM('dapprp','trt1','trt2','immorp',ICAapp)
  expect_equal(round(out[["ATE1"]],4),0.1949)
  expect_equal(round(out[["ATE2"]],4),0.3202)
  expect_equal(round(out[["ACME1"]],4),0.1132)
  expect_equal(round(out[["ACME2"]],4),0.1770)
  expect_equal(round(out[["estimand1"]],4),1.5629)
  expect_equal(round(out[["estimand2"]],4),0.9516)
  expect_equal(round(out[["confidenceIntervals"]],4),structure(list(ATE1.ci = c(0.1361, 0.2478),
                                                                    ATE2.ci = c(0.2645,0.3757),
                                                                    ACME1.ci = c(0.0751, 0.1509),
                                                                    ACME2.ci = c(0.139, 0.2147),
                                                                    estimand1.ci = c(1.1864, 2.2052),
                                                                    estimand2.ci = c(0.7456, 1.224)),
                                                               row.names = c("2.5%", "97.5%"), class = "data.frame"))
  set.seed(1)
  out = getCCM('nvotep','trt1','trt2','immorp',ICAapp,noInteraction = FALSE)
  expect_equal(round(out[["ATE1"]],4),0.182)
  expect_equal(round(out[["ATE2"]],4),0.2811)
  expect_equal(round(out[["ACME1"]],4),0.0964)
  expect_equal(round(out[["ACME2"]],4),0.1763)
  expect_equal(round(out[["estimand1"]],4),1.8289)
  expect_equal(round(out[["estimand2"]],4),1.1843)
  expect_equal(round(out[["confidenceIntervals"]],4),structure(list(ATE1.ci = c(0.1258, 0.2407),
                                                                    ATE2.ci = c(0.2274,0.3376),
                                                                    ACMET1.ci = c(0.0635, 0.1312),
                                                                    ACMET2.ci = c(0.138, 0.2156),
                                                                    estimand1.ci = c(1.3239, 2.6987),
                                                                    estimand2.ci = c(0.911,1.5571)),
                                                               row.names = c("2.5%", "97.5%"), class = "data.frame"))




})









