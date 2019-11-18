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

########################################
#                                      #
# Developing Example for Output Test   #
########################################
set.seed(2)
n = 16
pt = runif(1, min = 0.25, max = 0.35)
e1 = pt - runif(1, min = 0, max = pt)
e2 = runif(1, min = pt, max = 1) - pt
x = rbinom(1, n, prob = pt)
w = runif(1)
a1 = 2;b1 = 4;a2 = 5;b2 = 6
w_post = (w  * (beta(a1 + x, b1 + n - x) / beta(a1, b1))) /
  (w * (beta(a1 + x, b1 + n - x) / beta(a1, b1)) +
     (1 - w) * (beta(a2 + x, b2 + n - x) / beta(a2, b2)))
a1_post = a1 + x; a2_post = a2 + x
b1_post = b1 + n - x;b2_post = b2 + n - x

design = sample(c("tpi", "mtpi", "mmtpi"), 1, prob = c(0.2, 0.4, 0.4))
if(design %in% c("tpi", "mtpi"))
{
  cutpoints = c(0, pt - e1, pt + e2, 1)
  if(design == "mtpi")
  {
    breaks_lower = 0
    breaks_upper = 0
  }
}else
{
  gap = e1 + e2
  breaks_lower = floor((pt - e1) / gap)
  breaks_upper = floor((1 - pt - e2)  / gap)
  cutpoints = c(0, pt - e1 - gap * (breaks_lower : 0) , pt + e2 + gap * (0 : breaks_upper) , 1)
}

eta = 0.9
check_object = decisiontpi(pt = pt, e1 = e1, e2 = e2, x = x, n = n, eta = eta, design = design, w = w, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
if(check_object == "DU")
{
  test_that("Output Check",{
    expect_equal(as.numeric(w_post * pbeta(pt, a1_post, b1_post, lower.tail = FALSE) + (1 - w_post) * pbeta(pt, a2_post, b2_post, lower.tail = FALSE) >= eta),
                 1, tol = 10^-20)
  })
}else
{
  if(design == "tpi")
  {
    p_E = w_post * (pbeta(pt - e1, a1_post, b1_post)) + (1 - w_post) *(pbeta(pt - e1, a2_post, b2_post))
    p_S = w_post * (pbeta(pt + e2, a1_post, b1_post) - pbeta(pt - e1, a1_post, b1_post)) +
      (1 - w_post) *(pbeta(pt + e2, a2_post, b2_post) - pbeta(pt - e1, a2_post, b2_post))
    p_D = w_post * (1 - pbeta(pt + e2, a1_post, b1_post)) + (1 - w_post) *(1 - pbeta(pt + e2, a2_post, b2_post))
    if(check_object == "D"){
      test_that("Output Check",{expect_equal(max(c(p_E, p_S, p_D)), p_D, tol = 10^-20)})
    }else if(check_object == "E"){
      test_that("Output Check",{expect_equal(max(c(p_E, p_S, p_D)), p_E, tol = 10^-20)})
    }else{
      test_that("Output Check",{expect_equal(max(c(p_E, p_S, p_D)), p_S, tol = 10^-20)})
    }
  }else{
    cutlength = length(cutpoints)
    upm_E_check = c()
    upm_D_check = c()
    upm_S_check = c()
    for(i in 1:(breaks_lower + 1))
    {
      upm_E_check[i] = UPM(w = w_post, a = cutpoints[i], b = cutpoints[i + 1], a1 = a1_post, b1 = b1_post, a2 = a2_post, b2 = b2_post)
    }
    upm_S_check = UPM(w = w_post, a = pt - e1, b = pt + e2, a1 = a1_post, b1 = b1_post, a2 = a2_post, b2 = b2_post)
    for(i in 1:(breaks_upper + 1))
    {
      upm_D_check[i] = UPM(w = w_post, a = cutpoints[breaks_lower + 2 + i], b = cutpoints[breaks_lower + 3 + i], a1 = a1_post, b1 = b1_post, a2 = a2_post, b2 = b2_post)
    }
    max_upm = max(c(upm_S_check, upm_E_check, upm_D_check))
    if(check_object == "D"){
      test_that("Output Check",{expect_equal(max_upm, max(upm_D_check), tol = 10^-20)})
    }else if(check_object == "E"){
      test_that("Output Check",{expect_equal(max_upm, max(upm_E_check), tol = 10^-20)})
    }else{
      test_that("Output Check",{expect_equal(max_upm, upm_S_check, tol = 10^-20)})
    }
  }

}
