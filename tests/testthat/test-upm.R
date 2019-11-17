
test_that("Error Check", {
  expect_error(UPM(w = 0.5, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
               "Please input model parameters for both priors properly")
})

test_that("Error Check", {
  expect_error(UPM(w = 1.2, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
               "w is weight taken on first prior (informative), which can lie between 0 and 1", fixed = TRUE)
})

test_that("Error Check", {
  expect_error(UPM(w = - 10^-5, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
               "w is weight taken on first prior (informative), which can lie between 0 and 1", fixed = TRUE)
})

test_that("Warning Message Check", {
  expect_warning(UPM(w = 1, a1 = 5, b1= 1, a2 = 2, b2 = 4),
                 "Check inputs for prior parameters, taking a1 and b1 as original parameters")
})



check_object = evaluate_promise(UPM(w = 0, a = 0.2, b = 0.8, a2 = 2, b2 = 4))
test_that("Warning Message Check", {
  expect_equal(check_object$warnings,
                 "You should put the parameter values for a1 and b1 instead of a2 and b2")
})

check_value = (pbeta(0.8, 2, 4) - pbeta(0.2, 2, 4)) / (0.8 - 0.2)
test_that("Output Check", {
  expect_equal(check_object$result,
                check_value , tol = 1^-14)
})
remove(check_object)




test_that("Error Check", {
  expect_error(UPM(w = 0.1, a1 = 2, b2 = 4),
               "Please input model parameters for both priors properly")
})

test_that("Error Check", {
  expect_error(UPM(w = 0.1, a1 = 2, b2 = 4),
               "Please input model parameters for both priors properly")
})


check_object = evaluate_promise(UPM(w = 0.1, a = -0.3, b = 0.7, a1 = 2, b1 = 6, a2 = 4, b2 = 4))

test_that("Warning message Check", {
  expect_equal(check_object$warnings,
               "Domain of Beta distribution is (0,1), changing a to 0", fixed = TRUE)
})

a_range = 0
b_range = 0.7
upm_value_check = (0.1 * (pbeta(b_range, shape1 = 2, shape2 = 6) - pbeta(a_range, shape1 = 2, shape2 = 6)) +
                  (1 - 0.1) * (pbeta(b_range, shape1 = 4, shape2 = 4) - pbeta(a_range, shape1 = 4, shape2 = 4)))/(b_range - a_range)

test_that("Output Check", {
  expect_equal(check_object$result,
              upm_value_check, tol = 10^-14)
})
