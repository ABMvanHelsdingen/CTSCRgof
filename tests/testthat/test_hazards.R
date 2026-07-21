
test_that("SCR hazard",{

  s <- c(0,0)
  C <- matrix(0, nrow = 4, ncol = 2)
  C[,1] <- c(0,1,0,1)
  C[,2] <- c(0,0,1,1)
  duration <- 4
  pars <- c(1,1)
  times <- c(1,2,3)
  out <- hazard_SCR(s,C,duration,pars,times,NULL)
  expect_equal(out$h[,1], rep(1,3))
  expect_equal(out$h[,2], rep(exp(-0.5),3))
  expect_equal(out$h[,3], rep(exp(-0.5),3))
  expect_equal(out$h[,4], rep(exp(-1),3))

  expect_equal(out$H[,1], 1 * 1:3)
  expect_equal(out$H[,2], exp(-0.5) * 1:3)
  expect_equal(out$H[,3], exp(-0.5) * 1:3)
  expect_equal(out$H[,4], exp(-1) * 1:3)
})

test_that("SESCR hazard",{

  s <- c(0,0)
  C <- matrix(0, nrow = 3, ncol = 2)
  C[,1] <- c(0,1,1)
  C[,2] <- c(0,0,1)
  duration <- 3
  pars <- c(2,1,1,0.5)
  times <- c(1,2,3)

  capthist <- matrix(0, nrow = 2, ncol = 3)
  capthist[,1] <- c(1,2)
  capthist[,2] <- c(1,1)
  capthist[,3] <- c(1,2)


  out <- hazard_SESCR(s,C,duration,pars,times,capthist)
  expect_equal(out$h[,1], c(2,2 +6*exp(-1), 2 - 2*exp(-1) + 8*exp(-3)))
  expect_equal(out$h[,2], c(2*exp(-0.5), 2*exp(-0.5) -2*exp(-1.5) + 8*exp(-3),
                            2*exp(-0.5) - 2*exp(-1.5) + 8*exp(-1)))
  expect_equal(out$h[,3], c(2*exp(-1), 2*exp(-1) - 2*exp(-2) + 8*exp(-5),
                            2*exp(-1) - 2*exp(-2) + 8*exp(-3)))

  H_desired <- matrix(0, nrow = 3, ncol = 3)
  H_desired[1, ] <- c(2, 2*exp(-0.5), 2*exp(-1))
  H_desired[2, ] <- 2* H_desired[1, ] + (1-exp(-1)) *
    c(6, 8*exp(-2) - 2*exp(-0.5), 8*exp(-4) - 2*exp(-1))
  H_desired[3, ] <- H_desired[1, ] + H_desired[2,] + (1-exp(-1)) *
    c(8*exp(-2) - 2, 8 - 2*exp(-0.5), 8*exp(-2) - 2*exp(-1))

  expect_equal(out$H, H_desired)
})

test_that("Cosine hazard",{

  s <- c(0,0)
  C <- matrix(0, nrow = 4, ncol = 2)
  C[,1] <- c(0,1,0,1)
  C[,2] <- c(0,0,1,1)
  duration <- 4*pi
  pars <- c(1,1,1)
  times <- c(1,2,3) * pi
  time_rates <- 1 + cos(2*pi*times)
  out <- hazard_cosine(s,C,duration,pars,times,NULL)
  expect_equal(out$h[,1], rep(1,3) * time_rates)
  expect_equal(out$h[,2], rep(exp(-0.5),3) * time_rates)
  expect_equal(out$h[,3], rep(exp(-0.5),3) * time_rates)
  expect_equal(out$h[,4], rep(exp(-1),3) * time_rates)
})
