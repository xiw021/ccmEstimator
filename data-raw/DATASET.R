## code to prepare `DATASET` dataset goes here

ICAapp = read.csv(file = system.file("extdata","ICAapp.csv",package = 'ccmEstimator'))
usethis::use_data(ICAapp, overwrite = TRUE)
