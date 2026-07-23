#' Get a starting estimate for activity centers
#'
#' A weighted average of the locations where an individual was observed
#' @inheritParams Hazards
#' @export

empiricalAC <- function(C, capthist){

  # cameras and ids must be 1-indexed


  ids <- capthist[,2]
  cameras <- capthist[,3]

  m = max(ids)
  K = nrow(C)

  # initial estimates for activity centers based on observed locations
  a = matrix(0, nrow = m, ncol= 3)
  for (i in 1:length(ids)){
    An = ids[i]
    k = cameras[i]
    a[An, 1:2] = a[An, 1:2] + C[k, ]
    a[An, 3] = a[An, 3] + 1
  }
  a[, 1:2] = a[, 1:2] / a[,3]
  a = a[, 1:2]

  return(a)
}

#' Check Capture History
#'
#' Checks if a supplied capthist matrix is in the expected format
#'
#' @inheritParams Hazards
#' @param K number of cameras or traps in a survey
#' @export

checkCaptHist <- function(capthist, duration, K){
  if (ncol(capthist) < 3){
    stop("capthist must have at least three columns")
  }

  times <- capthist[,1]
  ids <- capthist[,2]
  cameras <- capthist[,3]

  for (i in 2:length(times)) {
    if ((times[i] - times[i - 1]) < 1.e-14)
      stop("times must be in ascending order with no simultaneous events")
  }

  if (max(times) > duration){
    stop("No detections can occur after the end of the survey")
  }

  if (min(times) < 0){
    stop("No detections can occur before time zero")
  }

  if (max(ids) != length(unique(ids))){
    stop("individuals must be numbered 1 to m, with no unused integers")
  }

  if (max(cameras) > K){
    stop("Incorrect trap ids")
  }

  if (min(cameras) < 1){
    stop("traps must be numbered from 1")
  }
}

#' Checks Inputs
#'
#' Checks if inputs for a hazard or simulate function are correct
#'
#' @inheritParams hazard_SCR
#' @export

checkBasic <- function(s, C, duration, times){

  if (duration <= 0){
    stop("duration must be a postive scalar")
  }

  if (!is.null(times)){

    if (min(times) < 0){
      stop("times must be non-negative")
    }

    if (max(times) > duration){
      stop("Cannot evaluate hazard after survey has ended")
    }
  }

  # CHECK ARGUMENTS
  if (!((length(s) == 2) || (ncol(s)==2))){
    stop("Activity center(s) must have an x and y coordinate")
  }

  if(ncol(C) != 2){
    stop("C must be a Kx2 matrix")
  }

}

#' Freeman-Tukey statistic
#'
#' Returns the Freeman-Tukey statistic, without the coefficient of 4.
#' @param obs The observed number of events, either a vector of matrix
#' @param exp The expected number of events, with the same dimensions
#' @export

FT_Stat <- function(obs, exp){
  if (length(class(obs)) == 1){ # Vectors
    stat = 0
    L = length(obs)
    for(i in 1:L){
      stat = stat + (sqrt(obs[i]) - sqrt(exp[i]))^2
    }
  } else { # Matrix
    R = nrow(obs)
    C = ncol(obs)
    stat = 0
    for(i in 1:R){
      for(j in 1:C){
        stat = stat + (sqrt(obs[i,j]) - sqrt(exp[i,j]))^2
      }
    }
  }
  return(stat)
}

#' Calculate hazard integrals
#'
#' Output is the comps input of `timeRescaling()` and `residuals()` functions. For multiple individuals.
#' @inheritParams timeRescaling
#' @param times times at which to evaluate the hazard and its integral
#' @export

Hazards <- function(S, C, duration, pars, times, hazard, N, capthist){

  out <- list()
  tau <- c(times, duration)
  CH <- ncol(capthist)
  for(ani in 1:N){
    Ij = sum(capthist[,2] == ani)
    capthist_ani <- matrix(0, nrow = Ij, ncol = CH)
    capthist_ani[1:Ij,] <- capthist[capthist[,2] == ani,]
    # Get hazard integrals for the individual
    out[[ani]] = hazard(S[ani,], C, duration, pars, tau, capthist_ani)$H
  }
  return(out)
}

#' Calculate hazard integrals at end of survey,
#'
#' Returns hazard integrals for each pair of individual and camera
#' @inheritParams timeRescaling
#' @export

Integrals <- function(S, C, duration, pars, hazard, N, capthist){

  N <- nrow(S)
  K <- nrow(C)
  CH <- ncol(capthist)
  out <- matrix(0, nrow = N, ncol = K)
  if(!is.null(capthist)){
    for(ani in 1:N){
      Ij = sum(capthist[,2] == ani)
      capthist_ani <- matrix(0, nrow = Ij, ncol = CH)
      capthist_ani[1:Ij,] <- capthist[capthist[,2] == ani,]
      # Get hazard integrals for the individual
      out[ani, ] = hazard(S[ani,], C, duration, pars, duration, capthist_ani)$H
    }
  } else {
    for(ani in 1:N){
      # Get hazard integrals for the individual
      out[ani, ] = hazard(S[ani,], C, duration, pars, duration, NULL)$H
    }
  }
  return(out)
}

#' Transition Matrix for SESCR
#'
#' Returns a matrix showing the expected number of "transitions" between cameras, for a single individual
#' Only valid for the SESCR Model. Called by `transMetric_prep_SESCR`.
#' @inheritParams hazard_SESCR
#' @param lambda0 detection rate at activity center
#' @param beta temporal decay of self-excitement
#' @param Dratio spatial decay of self-excitement
#' @param sigma home range size
#' @export

trans_matrix <- function(s, C, lambda0, beta, Dratio, sigma){

  K <- nrow(C)
  d <- Dratio * sigma

  acs <- s[rep(1,K),]
  locis2 <- (acs[,1] - C[,1])^2 + (acs[,2] - C[,2])^2
  baseline <- exp(-locis2/(2*sigma^2))

  transitions <- matrix(0, nrow = K, ncol = K)

  # Set up points for numerical integration
  points <- c(seq(0.0005,0.9995,length.out=1000),seq(1.05,19.95,length.out=190),
              seq(20.5,99.5,length.out=80))
  weights <- c(rep(0.001,1000),rep(0.1,190),rep(1,80))
  decays <- exp(-beta*points)
  S <- length(weights)

  for(c in 1:K){

    sc <- matrix(c(C[c,]),nrow=1); scs <- sc[rep(1,K),]

    locis2S <- (scs[,1] - C[,1])^2 + (scs[,2] - C[,2])^2
    spike <- exp(-locis2S/(2*d^2)) / (Dratio^2)

    for(j in 1:K){

      if (locis2[j] > 5 * sigma^2){
        # If distance > 5*sigma, do not calculate
        transitions[c, j] <- 1e-10
        next
      }

      S1 <- spike[j]; H1 <- baseline[j]
      S2 <- sum(spike[-j]); H2 <- sum(baseline[-j])
      D1 <- H1 - S1; D2 = H2 - S2

      Lambda1_vec = lambda0*(H1*points + (decays - 1)*(D1)/beta)
      Lambda2_vec = lambda0*(H2*points + (decays - 1)*(D2)/beta)
      integrands = exp(-Lambda1_vec - Lambda2_vec) * (H1 - decays*(D1))

      transitions[c, j] <- sum(integrands * weights) * lambda0
    }
  }

  # normalize rows
  for(c in 1:K){
    transitions[c, ] <- transitions[c, ] / sum(transitions[c, ])
  }

  return(transitions)
}

# Chen (2016): https://web.archive.org/web/20240502132513/https://www.math.fsu.edu/~ychen/research/multiHawkes.pdf
# Laub (2014): https://laub.au/pdfs/honours_thesis.pdf
#' Helper function for simulate_SESCR()
#'
#' Helper function that simulates from SESCR for a single individual
#' Not designed to be called except within simSESCR()
#' @param mu baseline rate for each stream
#' @param beta temporal decay of self-excitement
#' @param d self-excitement over space (NxN)
#' @param mu0 volume of baseline rate across all streams
#' @param K number of cameras and streams
#' @param t duration of survey
#' @export

SESCROgata <- function(mu, beta, d, mu0, K, t){

  # Custom Ogata algorithm for one individual in an SESCR model
  # This is a thinning algorithm, we simulate event times assuming the maximum possible
  # intensities, then remove some points to account for the lower intensities that actually occur.
  # following Chen (2016) and the definition of A in Laub (2014)


  # calculate maximum intensity of each stream/camera
  max_lambda = numeric(K)
  for(i in 1:K){
    max_lambda[i] = max(mu[i], mu0 * d[i,i])
  }
  lambda_bar = sum(max_lambda)

  times = numeric(0); streams = numeric(0)
  s = 0; n = 0; A = 0
  lambda_s = mu
  last_event = 0 # stream of last event
  while (s < t){
    u <- stats::runif(1)
    w <- -1*log(u)/lambda_bar
    s <- s + w # candidate event
    if (n > 0){ # When n=0, no prior events
      # Markovian, no sum across all prior events
      A <- exp(-beta*(s-max(times)))
      #for(i in 1:K){
        # lambda_s is intensity at candidate point s
        #lambda_s[i] = ((1 - A) * mu[i]) + (A * mu0 * d[i, last_event])
      #}
      lambda_s = ((1 - A) * mu) + (A * mu0 * d[, last_event])
    }

    # We generate a random number D, if D multiplied by the sum of the maximum intensities
    # exceeds the sum of the current intensities, the candidate event is thinned.
    # Otherwise, the probability that an event is assigned to each stream is proportional to their intensities
    D <- stats::runif(1)
    if (D*lambda_bar <= sum(lambda_s)){
      k = 1
      while(D*lambda_bar > sum(lambda_s[1:k])){
        k = k + 1
      }
      n <- n + 1
      last_event <- k
      times <- append(times, s)
      streams <- append(streams, k)
    }
  }
  # Check if last event > t, and return output
  if (n == 0 || (n == 1 && max(times) > t)){
    return(list(times=numeric(0), streams=numeric(0)))
  } else if (max(times)<t){
    return(list(times=times, streams=streams))
  } else {
    return(list(times=times[1:(n-1)], streams=streams[1:(n-1)]))
  }
}

#' Converts individual ids
#'
#' Converts capthist from those produced by simulation functions
#' Converts capthist to the form required by CheckCaptHist()
#' In particular, individuals must be numbered 1:m without unused ids
#'
#' The simulate_* functions used the row numbers of the activity center matrix S to number ids.
#' This results in individuals with integers in the range of 1:N, but when there are undetected individuals, there are unused numbers.
#' The timeRescaling() and runPPC() functions require detected individuals to have the lowest numbers, starting at one.
#'
#' @param ids id number of the individuals captured at detections
#' @param N is total population size
#' @returns ids, giving the converted individual IDs, and new_order, how the individuals should be ordered for use with timeRescaling() and runPPC()
#' @export

convertCaptHist <- function(ids, N){

  unique <- sort(unique(ids))
  ids2 = ids

  for(i in 1:length(unique)){
    ids2[ids == unique[i]] <- i
  }

  new_order <- c(unique, setdiff(1:N, ids))

  return(list(ids = ids2, new_order = new_order))
}

#' Subsampling
#'
#' Subsamples a sequence of numbers that under a null hypothesis are exponentially distributed
#'
#' Size of subsample defaults to n^(2/3) as suggested by Reynaud-Bouret et al (2014)
#' @param values sequence of compensator differences
#' @param n size of subsample
#' @returns p-value from Kolmogorov-Smirnov test on subsample, null is Exp(1)
#' @export

subSampleTest <- function(values, n=NULL){

  indexes = sample(1:length(values), size = round(length(values)^(2/3)))
  ks <- stats::ks.test(values[indexes], "pexp")$p.value
  return(ks)
}

#' Make a visually meaningful residual plot
#' @import ggplot2
#' @param resids the residuals
#' @export

residualPlot <- function(resids){

  data <- data.frame(resid = resids)
  pl <- ggplot2::ggplot(data, aes(resid)) +
    ggplot2::stat_ecdf(geom = "point", pad = FALSE, size = 2) +
    ggplot2::geom_abline(slope=1,intercept=0) +
    ggplot2::coord_flip() +
    ylim(0,1) +
    xlim(0,1) +
    ggplot2::theme_classic() +
    xlab("Residuals") +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18)) +
    ylab("Quantiles")
  return(pl)
}
