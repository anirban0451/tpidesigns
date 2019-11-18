##########################################
#                                        #
# Generation of Data for testing Purpose #
#                                        #
##########################################
set.seed(2)
nmax = sample(3:16, 1)
pt = runif(1, min = 0.25, max = 0.35)
e1 = pt - runif(1, min = 0, max = pt)
e2 = runif(1, min = pt, max = 1) - pt
w = runif(1)
a1 = 2;b1 = 4;a2 = 5;b2 = 6
design = sample(c("tpi", "mtpi", "mmtpi"), 1, prob = c(0.2, 0.4, 0.4))
eta = 0.9
test_that("Sample Size is conformable", {
  expect_error(tpitable(nmax = sample.int(2, 1), design = design, pt = pt, e1 = e1, e2 = e2, eta = eta, w, a1 = a1, b1 = b1, a2 = a2, b2 = b2),
               "Number of patients must be atleast 3")
})

check_objective =  as.matrix(tpitable(nmax = nmax, design = design, pt = pt, e1 = e1, e2 = e2, eta = eta, w, a1 = a1, b1 = b1, a2 = a2, b2 = b2), dimnames = NULL)
pos = sample(1:4, 1)
n = sample(3:nmax, 1)
n_pos = n - 2
out = check_objective[ ,n_pos]
decisions = c("E", "S", "D", "DU")
test_that("Value Check",{
  expect_equal(decisiontpi(pt = pt, e1 = e1, e2 = e2, x = out[4],n = n, w = w, eta = eta, a1 = a1, b1 = b1, a2 = a2, b2 = b2),
               "DU")
})
if(!is.na(out[pos]))
{
  test_that("Output Check", {
    expect_equal(decisiontpi(pt = pt, e1 = e1, e2 = e2, design = design, x = out[pos],n = n, w = w, eta = eta, a1 = a1, b1 = b1, a2 = a2, b2 = b2),
                 decisions[pos])
  })
}
