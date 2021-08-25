context("test checkInteraction
      function")

test_that("whether the checkInteraction function gives correct output using ICAapp data",{
  data(ICAapp)
  final.dat <- checkData(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp)
  expect_equal(checkInteraction(final.dat,sigLevel = 0.05),0)
  final.dat <- checkData(Y = "nvotep", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp)
  expect_equal(checkInteraction(final.dat,sigLevel = 0.05),1)
})









