
<!-- README.md is generated from README.Rmd. Please edit that file -->
# ccmEstimator

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/xiw021/ccmEstimator.svg?branch=master)](https://travis-ci.com/xiw021/ccmEstimator) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/xw2510/ccmEstimator?branch=master&svg=true)](https://ci.appveyor.com/project/xw2510/ccmEstimator) [![Codecov test coverage](https://codecov.io/gh/xiw021/ccmEstimator/branch/master/graph/badge.svg)](https://codecov.io/gh/xiw021/ccmEstimator?branch=master) <!-- badges: end -->

The goal of ccmEstimator is to perform comparative causal mediation analysis to compare the mediation effects of different treatments via a common mediator.

## Installation

You can install the released version of ccmEstimator from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ccmEstimator")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ccmEstimator)
## basic example code
# load example data. ICAapp is a data frame with 14 variables. It is the application data used to illustrate comparative causal mediation analysis methods in Bansak (2020).
data(ICAapp)
set.seed(321, kind = "Mersenne-Twister", normal.kind = "Inversion")
ccm.results <-
   getCCM(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", data = ICAapp,
   noInteraction = TRUE, sigLevel = 0.05, boots = 1000)
#> Estimation assuming no interactions between treatments and mediator.
#> The estimated mediation effect for treatment 2 is 1.563 times larger than the mediation effect for treatment 1 
#>    (with 95% CI: [1.183,2.16])
#> The estimated proportion mediated for treatment 2 is 0.952 times the size of that for treatment 1 
#>    (with 95% CI: [0.755,1.236])
#> Please use summary() to view the results in more detail.

summary(ccm.results)
#> Call:
#> getCCM(Y = "dapprp", T1 = "trt1", T2 = "trt2", M = "immorp", 
#>     data = ICAapp, noInteraction = TRUE, sigLevel = 0.05, boots = 1000)
#> 
#>       point.estimate lower.ci   upper.ci 
#> ATE1  0.1949336      0.1444899  0.2509022
#> ATE2  0.3201672      0.2655617  0.3794662
#> ACME1 0.1132405      0.07559316 0.1502631
#> ACME2 0.1769843      0.1382323  0.2186979
#> 
#>                           point.estimate lower.ci  upper.ci
#> ACME2/ACME1               1.562907       1.182749  2.159821
#> (ACME2/ATE2)/(ACME1/ATE1) 0.9515747      0.7549031 1.235673
```
