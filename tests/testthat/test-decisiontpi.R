test_that("Error Check", {
  expect_error(decisiontpi(pt = 0.3, e1 = 0.35, e2 = 0.05, x = 4 , n = 6, eta = 0.2, design = c("tpi", "mtpi", "mmtpi"), w = 0.2, a1 = 4, b1 = 5, a2 = 6, b2 = 7),
               "e1 and e2, two thresholds should be small compared to the target probability pt")
})

test_that("Error Check", {
  expect_error(decisiontpi(pt = 0.3, e1 = 0.02, e2 = 0.75, x = 4 , n = 6, eta = 0.2, design = c("tpi", "mtpi", "mmtpi"), w = 0.2, a1 = 4, b1 = 5, a2 = 6, b2 = 7),
               "e1 and e2, two thresholds should be small compared to the target probability pt")
})

test_that("Error Check", {
  expect_error(decisiontpi(pt = 1.001, e1 = 0.02, e2 = 0.05, x = 4 , n = 6, eta = 0.2, design = c("tpi", "mtpi", "mmtpi"), w = 0.2, a1 = 4, b1 = 5, a2 = 6, b2 = 7),
               "Target toxicity Probability should take values between 0 and 1")
})

test_that("Error Check", {
  expect_error(decisiontpi(pt = 0.3, e1 = 0.02, e2 = 0.05, x = 14 , n = 6, eta = 0.2, design = c("tpi", "mtpi", "mmtpi"), w = 0.2, a1 = 4, b1 = 5, a2 = 6, b2 = 7),
               "Number of patients experiencing DLT 's must be less than or equal to total number of patients treated in that cohort.", fixed = TRUE)
})

test_that("Warning Check", {
  expect_warning(decisiontpi(pt = 0.3, e1 = 0.02, e2 = 0.05, x = 4 , n = 6, eta = 0.2, design = c("tpi", "mtpi", "mmtpi"), w = 1, a2 = 6, b2 = 7),
               "You should put the parameter values for a1 and b1 instead of a2 and b2", fixed = TRUE)
})

