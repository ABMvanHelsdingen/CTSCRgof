#' Standard SCR Hazard function
#'
#' Calculates the hazard and integrals of the hazard for the simplest possible SCR model
#'
#' Returns a list with elements `h` and `H`, each being a matrix with a row for each value of times
#' and a column for each trap. The values of `h` are hazards and `H` are integrals of `h` since time zero.
#'
#' \eqn{\lambda_0} is the rate of detections at the activity center, and \eqn{\sigma} the spatial scale.
#'
#' @param s activity center of a single individual
#' @param C matrix of camera or trap locations
#' @param duration duration of survey
#' @param pars vector of pars: \eqn{\lambda_0} and \eqn{\sigma}, see details
#' @param times times at which to evaluate hazards and their integrals
#' @param capthist 3-column matrix giving times, ids and trap ids for each detection.
#'  `capthist` must be have continuous detection times in the first column,
#'  individual IDs numbered 1:m in the second, and trap IDs numbered 1:K in the third
#' @returns two matrices with hazards and integrals, see details
#' @export

hazard_SCR <- function(s, C, duration, pars, times, capthist = NULL){

  checkBasic(s, C, duration, times)
  if (!(all(pars > 0))){
    stop("All parameters must be positive")
  }

  K <- nrow(C)
  S <- t(replicate(K, s))

  lambda0 <- pars[1]
  sigma <- pars[2]


  loci = S - C
  dist2 = loci[,1]^2 + loci[,2]^2

  rates = lambda0 * exp(-dist2/ (2 * sigma^2))
  compensators = times %*% t(rates)
  h_out = matrix(rates, nrow = length(times), ncol = length(rates), byrow = TRUE)
  return(list(h = h_out, H = compensators))


}

#' SESCR Hazard function
#'
#' Hazard for the SESCR model of van Helsdingen & Jones-Todd, 2026
#'
#' Returns a list with components `h` and `H`, each being a matrix with a row for each value of times
#' and a column for each trap. The values of `h` are hazards and `H` are integrals of `h` since time zero.
#'
#' \eqn{\lambda_0} is the rate of detections at the activity center when there has been no prior detections. \eqn{\sigma} is the home range size.
#' \eqn{\beta} is the temporal decay rate for self-excitement and \eqn{r} the spatial decay.
#'
#' @inheritParams hazard_SCR
#' @param pars vector of pars: vector of pars: \eqn{\lambda_0}, \eqn{\sigma}, \eqn{\beta} and \eqn{r}, see details
#' @export

hazard_SESCR <- function(s, C, duration, pars, times, capthist){

  checkBasic(s, C, duration, times)
  if (!(all(pars > 0))){
    stop("All parameters must be positive")
  }

  # Extract parameters
  lambda0 <- pars[1]
  sigma <- pars[2]
  beta <- pars[3]
  r <- pars[4]

  if (r > 1){
    warning("r should not usually exceed 1")
  }

  # If individual is never detected
  if (nrow(capthist) == 0){
    return(hazard_SCR(s, C, duration, pars, times, NULL))
  }


  # Extract capthist
  tau <- capthist[,1]
  cameras <- capthist[,3]

  # If not all tau are within times, add them
  if (!all(tau %in% times)){
    times2 = unique(c(times, tau))
    times2 = sort(times2)
    indices = which(times2 %in% times)
    times = times2
  } else {
    indices = 1:length(times)
  }

  loci <- numeric(2)
  K <- nrow(C)
  S <- t(replicate(K, s))
  I <- length(times)

  loci = S - C
  dist2 = loci[,1]^2 + loci[,2]^2
  DF = exp(-dist2/ (2 * sigma^2)) # Detectability at each animal-camera combination

  #Matrix of spatial self-excitement between cameras
  sse = matrix(0, nrow = K, ncol = K)
  d = r * sigma
  for(i in 1:K){
    for(k in i:K){
      loci = C[i, ] - C[k, ]
      dist2 = loci[1]^2 + loci[2]^2
      sse[i,k] = exp(-dist2/(2*d^2)) / (r^2)
      sse[k,i] = sse[i,k] # Matrix is symmetric
    }
  }


  # A[j] determines the self-exciting effect at event j
  # The other variables are used to calculate the sums of the integral for each stream
  A = numeric(I)
  decays <- spikes <- matrix(0, nrow = I, ncol = K)
  Clast = 1 # Last camera that detected a given animal

  for(j in 2:I){
    if (times[(j-1)] %in% tau){
      # Update current camera
      Clast = cameras[which(tau == times[j-1])]
      A[j] = exp(-beta * (times[j] - times[(j - 1)])) # Markovian, no dependency on previous events

      for(k in 1:K){
        decays[j, k] = DF[k] * (1 - A[j])
        spikes[j, k] = sse[Clast, k] * (1 - A[j])
      }

    } else {
      # Last detection was not of the individual
      A[j] = exp(-beta * (times[j] - times[(j - 1)])) * A[(j-1)]

      for(k in 1:K){
        decays[j, k] = DF[k] * (A[(j-1)] - A[j])
        spikes[j, k] = sse[Clast, k] * (A[(j-1)] - A[j])
      }
    }
  }

  # Output hazards and integrals

  H <- h <- matrix(0, nrow = I, ncol = K)
  Clast <- 1
  for(j in 1:I){
    # Update last camera
    if (j >= 2 && times[(j-1)] %in% tau){
      Clast = cameras[which(tau == times[j-1])]
    }

    for(k in 1:K){
      # If A[j] = 0, then Clast=1 but is irrelevant
      h[j, k] <- ((1 - A[j])*DF[k]) + (A[j] * sse[Clast, k])
      h[j, k] <- lambda0 * h[j, k]
      H[j, k] <- sum(spikes[1:j,k] - decays[1:j,k]) * (lambda0/beta)
      H[j, k] <- H[j, k] + lambda0 * DF[k]*times[j]
    }
  }

  return(list(h = h[indices,], H = H[indices,]))

}

#' Cosine SCR Hazard function
#'
#' Calculates the hazard and integrals of the hazard for a cosine wave detection function.
#'
#' Returns a list with elements `h` and `H`, each being a matrix with a row for each value of times
#' and a column for each trap. The values of `h` are hazards and `H` are integrals of `h` since time zero.
#'
#'
#' \eqn{\lambda_0} is the rate of detections at the activity center, and \eqn{\sigma} the spatial scale. \cr
#' The hazard of standard SCR is multiplied by \eqn{1 + A\;\text{cos}(2\pi F(t - P))} \cr
#' \eqn{A} is the magnitude of the wave, 0 = standard SCR and 1=max possible. \cr
#' \eqn{F} is the frequency of the wave \cr
#' \eqn{P} is the phase shift. The maxima of the wave occur at P + \eqn{\frac{k}{F}} for k = 0,1,2...
#'
#' @inheritParams hazard_SCR
#' @param pars vector of pars: \eqn{\lambda_0}, \eqn{\sigma}, \eqn{A}, \eqn{F}, \eqn{P}
#' @returns two matrices with hazards and integrals, see details
#' @export

hazard_cosine <- function(s, C, duration, pars, times, capthist = NULL){

  checkBasic(s, C, duration, times)
  if (!(all(pars[1:2] > 0))){
    stop("Lambda0 and sigma must be positive")
  }
  if ((pars[3] > 1) || (pars[3] < 0)){
    stop("A must be between 0 and 1 (inclusive)")
  }
  if (pars[4] <= 0){
    stop("F must be positive")
  }
  if ((pars[5] <= 0) || (pars[5] > (1/pars[4]))){
    stop("P must be between 0 and 1/F")
  }

  K <- nrow(C)
  S <- t(replicate(K, s))

  lambda0 <- pars[1]
  sigma <- pars[2]
  A <- pars[3]
  Fr <- pars[4]
  P <- pars[5]



  loci = S - C
  dist2 = loci[,1]^2 + loci[,2]^2

  trap_rates = exp(-dist2/ (2 * sigma^2))
  time_rates = lambda0 * (1 + A * cos(2*pi*Fr*(times - P)))
  time_ints = lambda0 * (times + (A * sin(pi*Fr*times) * cos(pi*Fr*(2*P - times)) / (Fr * pi)))
  h = time_rates %*% t(trap_rates)
  H = time_ints %*% t(trap_rates)
  return(list(h = h, H = H))

}

