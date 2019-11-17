test_that("Error Check", {
  testthat::expect_error(weights_formulate(w = 0.5, x = 5, n = 10, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
                        "Please input model parameters  for both priors properly")
})
