test_that("chain metrics",{


  capthist <- matrix(0, nrow = 22, ncol = 3)
  capthist[,1] <- 1:22
  capthist[,2] <- c(rep(1,9), rep(2,7), rep(3,6))
  capthist[,3] <- c(1,1,1,1,2,3,2,3,3,2,3,3,3,1,2,3,1,2,3,1,3,1)

  out <- chainMetrics(NULL, capthist, NULL, NULL, NULL)

  expect_equal(out, c(4, 8/3, 3))

})

test_that("Basic FT",{

  capthist <- matrix(0, nrow = 22, ncol = 3)
  capthist[,1] <- 1:22
  capthist[,2] <- c(rep(1,9), rep(2,7), rep(3,6))
  capthist[,3] <- c(1,1,1,1,2,3,2,3,3,2,3,3,3,1,2,3,1,2,3,1,3,1)

  expected <- matrix(0, nrow = 3, ncol = 3)
  expected[1, ] <- c(9,3,4)
  expected[2, ] <- c(4,2,1)
  expected[3, ] <- c(3,3,3)

  out <- FTMetrics(expected, capthist, 3, 3, NULL)

  expect_equal(out[1], 24 - 4*sqrt(2)*sqrt(3) - 6*sqrt(3))
  expect_equal(out[2], 16 - 6*sqrt(6))
  expect_equal(out[3], 54 - 14*sqrt(8) - 2*sqrt(5)*sqrt(8))



})

test_that("trans FT",{

  capthist <- matrix(0, nrow = 22, ncol = 3)
  capthist[,1] <- 1:22
  capthist[,2] <- c(rep(1,9), rep(2,7), rep(3,6))
  capthist[,3] <- c(1,1,1,1,2,3,2,3,3,2,3,3,3,1,2,3,1,2,3,1,3,1)

  expected <- matrix(0, nrow = 3, ncol = 3)
  expected[1, ] <- c(6,4,3)
  expected[2, ] <- c(9,1,4)
  expected[3, ] <- c(3,5,1)

  out <- transMetric(expected, capthist, 3, 3, NULL)

  expect_equal(out, 35 - 2*sqrt(6)*sqrt(3) - 10*sqrt(3))

})

test_that("trans FT prep_SCR",{

  s <- c(0,0)
  C <- matrix(0, nrow = 4, ncol = 2)
  C[,1] <- c(0,1,0,1)
  C[,2] <- c(0,0,1,1)
  duration <- 4
  pars <- c(1,1)
  times <- c(1,2,4)

  Hs <- hazard_SCR(s, C, duration, pars, times, NULL)$H[3, ]
  H <- sum(Hs)

  P_desired <- (Hs %*% t(Hs)) * (H - 1 + exp(-H)) / (H^2)

  S <- matrix(c(0,0), nrow = 1, ncol = 2)

  out <- transMetric_prep_SCR(S, C, duration, pars, hazard_SCR, NULL)

  expect_equal(out, P_desired)

})
