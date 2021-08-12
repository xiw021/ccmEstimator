context("test getCCM function")

test_that("whether the getCCM function gives correct output using ICAapp data",{
  set.seed(1)
  icadata = testdata = read.csv(file = system.file("extdata","ICAapp_data.csv",package = 'ccmEstimator'))
  out = getCCM('dapprp','trt1','trt2','immorp',icadata)
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
})
