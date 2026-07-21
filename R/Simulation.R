#' Simulates standard continuous-time SCR
#'
#' All individuals and terrain homogeneous
#' @inheritParams hazard_SCR
#' @param S matrix of activity centers
#' @returns individuals are numbered according to their row in `S`.
#' @export

simulate_SCR <- function(S, C, duration, pars){

  checkBasic(S, C, duration, NULL)
  if (!(all(pars > 0))){
    stop("All parameters must be positive")
  }

  K <- nrow(C)
  N <- nrow(S)
  n <- 0 # number of animals with observations

  lambda0 <- pars[1]
  sigma <- pars[2]


  # empty vectors for storing simulations
  times <- cameras <- animals <- numeric(0)

  # for each animal:
  for(i in 1:N){
    # Calculate mu for MHP
    mus <- numeric(K)
    for(k in 1:K){
      distance <- C[k, ] - S[i, ]
      d2 <- sum(distance^2)
      mus[k] <- exp(-d2 / (2 * sigma^2))
    }
    mus <- lambda0 * mus

    # Generate MHP
    result <- IHSEP::simPois(int = function(x){sum(mus)}, cens = duration, int.M = sum(mus))

    if (length(result) > 0){ # If animal is detected
      times <- append(times, result)
      animals <- append(animals, rep(i, length(result)))
      cameras <- append(cameras, sample(c(1:K), size = length(result), replace = TRUE,
                                        prob = mus/sum(mus)))
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  capthist <- cbind(times, animals, cameras)
  return(capthist = capthist)
}



#' Simulates from the SESCR model
#' @inheritParams hazard_SESCR
#' @inheritParams simulate_SCR
#' @returns individuals are numbered according to their row in `S`.
#' @export
simulate_SESCR <- function(S, C, duration, pars){

  checkBasic(S, C, duration, NULL)

  # Extract parameters
  lambda0 <- pars[1]
  sigma <- pars[2]
  beta <- pars[3]
  r <- pars[4]

  if (r > 1){
    warning("r should not usually exceed 1")
  }

  # Constants
  K <- nrow(C)
  N <- nrow(S)
  n <- 0 # number of animals with observations
  d <- r * sigma

  # Calculate spatial self-excitement
  sse <- matrix(0, nrow = K, ncol = K)
  for(x in 1:K){
    for(y in 1:x){
      distance <- C[x, ] - C[y, ]
      d2 <- sum(distance^2)
      sse[x,y] <- exp(-d2 / (2 * d^2)) * (sigma^2 / d^2)
      sse[y,x] <- sse[x,y] # sse is symmetric
    }
  }

  # empty vectors for storing simulations
  times <- cameras <- animals <- numeric(0)

  for(i in 1:N){ # for each individual
    # Calculate mu for MHP
    mus <- numeric(K)
    for(k in 1:K){
      distance <- C[k, ] - S[i, ]
      d2 <- sum(distance^2)
      mus[k] <- exp(-d2 / (2 * sigma^2))
    }

    mus <- lambda0 * mus

    # Generate MHP
    result <- SESCROgata(mu = mus, beta = beta,
                          d = sse, mu0 = lambda0, K = K, t = duration)
    if (length(result$times) > 0){ # If animal is detected
      times <- append(times, result$times)
      cameras <- append(cameras, result$streams)
      animals <- append(animals, rep(i, length(result$times)))
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  capthist <- cbind(times, animals, cameras)
  return(capthist = capthist)
}

#' Simulates continuous-time SCR with a cosine wave function
#'
#' All individuals and terrain homogeneous
#' @inheritParams simulate_SCR
#' @param pars vector of pars: \eqn{\lambda_0}, \eqn{\sigma} and \eqn{\phi}
#' @returns individuals are numbered according to their row in `S`.
#' @export

simulate_cosine <- function(S, C, duration, pars){

  checkBasic(S, C, duration, NULL)
  if (!(all(pars[1:2] > 0))){
    stop("Lambda0 and sigma must be positive")
  }
  if ((pars[3] > 1) || (pars[3] < 0)){
    stop("phi must be between 0 and 1 (inclusive)")
  }

  K <- nrow(C)
  N <- nrow(S)
  n <- 0 # number of animals with observations

  lambda0 <- pars[1]
  sigma <- pars[2]
  phi <- pars[3]


  # empty vectors for storing simulations
  times <- cameras <- animals <- numeric(0)

  # for each animal:
  for(i in 1:N){
    # Calculate mu for MHP
    mus <- numeric(K)
    for(k in 1:K){
      distance <- C[k, ] - S[i, ]
      d2 <- sum(distance^2)
      mus[k] <- exp(-d2 / (2 * sigma^2))
    }
    mus <- lambda0 * mus

    for(k in 1:K){

      # Generate cosine Poisson point process
      result <- IHSEP::simPois(int = function(x){mus[k] * (1 + (phi * cos(2*pi*x)))},
                               cens = duration, int.M = mus[k] * (1 + phi))

      if (length(result) > 0){ # If animal is detected at that camera
        times <- append(times, result)
        animals <- append(animals, rep(i, length(result)))
        cameras <- append(cameras, rep(k, length(result)))
      }
    }
  }
  times_order = order(times, decreasing = FALSE)
  times <- times[times_order]
  cameras <- cameras[times_order]
  animals <- animals[times_order]
  capthist <- cbind(times, animals, cameras)
  return(capthist = capthist)
}
