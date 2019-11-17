test_that("Error Check", {
  expect_error(weights_formulate(w = 0.5, x = 5, n = 10, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
                        "Please input model parameters  for both priors properly")
})

test_that("Error Check", {
  expect_error(weights_formulate(w = 1.2, x = 15, n = 10, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
               "Weight on informative prior can be at most 1")
})

test_that("Error Check", {
  expect_error(weights_formulate(w = - 10^-5, x = 5, n = 10, a1 = NULL, b1= 1, a2 = 2, b2 = 4),
               "Weight on a prior can not be negative")
})

test_that("Warning Message Check", {
  expect_warning(weights_formulate(w = 1, x = 5, n = 10, a1 = 5, b1= 1, a2 = 2, b2 = 4),
                 "Check inputs for prior parameters, taking a1 and b1 as original parameters")
})

test_that("Warning Message Check", {
  expect_warning(weights_formulate(w = 0, x = 5, n = 10, a2 = 2, b2 = 4),
                 "You should put the parameter values for a1 and b1 instead of a2 and b2")
})

test_that("Error Check", {
  expect_error(weights_formulate(w = 0.1, x = 5, n = 10, a1 = 2, b2 = 4),
               "Please input model parameters  for both priors properly")
})


test_that("Output Check", {
  expect_equal(weights_formulate(w = 1, x = 5, n = 10, a1 = 2, b1 = 4)$param_inform,
               c(7,9))
})


n = 15
x = sample.int(n, size = 1)
w = runif(1)
a1 = 5
a2 = 6
b1 = 0.9
b2 = 0.1

w_post = (w  * (beta(a1 + x, b1 + n - x) / beta(a1, b1))) /
  (w * (beta(a1 + x, b1 + n - x) / beta(a1, b1)) +
     (1 - w) * (beta(a2 + x, b2 + n - x) / beta(a2, b2)))
check_object = weights_formulate(w = w, x = x, n = n, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
test_that("Output Check", {
  expect_equal(check_object$weight,
               w_post, tol = 10^-16)
})
test_that("Output Check", {
  expect_equal(check_object$param_inform,
               c(a1 + x, b1 + n -x), tol = 10^-10)
})
test_that("Output Check", {
  expect_equal(check_object$param_noninform,
               c(a2 + x, b2 + n -x), tol = 10^-10)
})
