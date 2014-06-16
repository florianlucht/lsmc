####################################################################################################
# heston
####################################################################################################
#' Heston model
#'
#' This method simulates discrete paths in the Heston model.
#' The price paths satisfy the following stochastic differential equations: 
#' \deqn{dS_t = S_t \left( r dt  + \sqrt{\nu_t} d W_t^{S}  \right) \ ,}
#'      {dS = S * ( r * dt  + sqrt{nu} * dW ) ,}
#' \deqn{d\nu_t = \kappa ( \theta - \nu_t ) dt + \xi \sqrt{\nu_t} dW_t^{\nu} \ ,}
#'      {d nu = kappa * ( theta - nu ) * dt + xi * sqrt{nu_t} * dB ,}
#' \deqn{\rho =  Kor\left(W_t^{S}, W_t^{\nu}\right) \ .}{rho =  Kor[ W, B ] .}
#'
#' @export
#'
#' @family models
#' 
#' @param initialValue initial value (scalar)
#' @param r            interest rate (scalar)
#' @param kappa        speed of adjustment (scalar)
#' @param theta        long run volatility (scalar)
#' @param rho          correlation (scalar)
#' @param xi           volatlility of volatility (scalar)
#' @param sigmaZero    initial volatility (scalar)
#' @param numDates     number of dates (scalar)
#' @param maturity     maturity (scalar)
#' @param numPaths     number of paths (scalar)
#' @param antithetic   logical. If \code{TRUE} antitethic paths are generated.
#' @param plot         logical. If \code{TRUE} paths are plotted.
#'
#' @return simulated paths as matrix with \code{numDates+1} rows and \code{numPaths} columns 
#'
#' @examples
#' # generate one simple path
#' S <- heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100)
#' 
#' # generate 100 paths with maturity T=10
#' S <- heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100, maturity=10, numPaths=100)
heston <- function(initialValue, r, kappa, theta, rho,
                   xi, sigmaZero, numDates, 
                   maturity=1, 
                   numPaths=1, 
                   stepSize=FALSE, 
                   antithetic=FALSE, 
                   plot=FALSE) {

  S <- matrix(0, numDates+1, numPaths)
  S[1, ] <- initialValue
  dt <- maturity / numDates
  X <- rep(initialValue, numPaths)
  
  if(antithetic) numPaths <- numPaths / 2
  
  if (identical(stepSize, FALSE)) {
    numSteps <- 1    
  } else {
    if((dt%%stepSize) != 0) stop('stepSize ist kein Teiler von dt')
    numSteps <- as.integer(dt / stepSize)
    dt <- stepSize
  }
  sqdt <- sqrt(dt)
  
  # current value of the volatility process
  V <- rep(sigmaZero, numPaths)
  
  for(n in 1:numDates)
  {
    for(i in 1:numSteps) {
      # generate Z,W with cov(Z,W)<-rho
      Z <- rnorm(numPaths)
      Y <- rnorm(numPaths)
      W <- rho * Z + sqrt(1 - rho^2) * Y
      
      truncatedV <- pmax(V,0)
      
      V <- V + kappa * (theta - truncatedV ) * dt + xi * sqrt(truncatedV) * sqdt * W
      
      if(antithetic) {
        Z <- c(Z, -Z)
        truncatedV <- c(truncatedV, truncatedV)
      }
      
      X <- X * exp((r - 0.5 * truncatedV) * dt + sqrt(truncatedV) * sqdt * Z)
    }    
    S[n+1, ] <- X
  }
  
  # verbose
  if (plot) {
    matplot(seq(0, maturity, maturity / numDates), S, type="l", xlab="Time") 
  }
  
  return(S)
}


####################################################################################################
# getHestonGenerator
####################################################################################################
#' Generator function of the Heston model 
#' 
#' This method returns a function which can be used to simulate one step in the Heston model for 
#' given parameters.  
#'
#' @export
#'
#' @param r            interest rate (scalar)
#' @param kappa        speed of adjustment (scalar)
#' @param theta        long run volatility (scalar)
#' @param rho          correlation (scalar)
#' @param xi           volatlility of volatility (scalar)
#' @param sigmaZero    initial volatility (scalar)
#' @param antithetic   logical. If \code{TRUE} antitethic paths are generated.
#' @param stepSize     step size in the Euler-Maruyama method for process simulation
#'
#' @return function with arguments:
#' \describe{
#'   \item{S}{current value(s) (scalar or vector)}
#'   \item{dt}{step size (scalar)}
#'   \item{notStopped}{optional index of paths which are not exercised (vector of logical). This 
#'                     argument is used internally to calculate the correct values if 
#'                     \code{antithetic=TRUE} is given.}
#'   \item{reset}{logical. If \code{TRUE} the volatility process will be reset to the initial 
#'                value on the next function call.}
#'
#' @examples
#' # get generator function
#' genHeston <- getHestonGenerator(0.05, 2, 0.4, -0.6, 0.4, 0.4)
#' # generate one step in the Heston model with step size 0.1 and current value 100
#' next <- genHeston(100, 0.1)
#' # reset volatility process
#' genHeston(NULL, NULL, reset=TRUE)
#' # genrate path with 100 steps 
#' S <- c(100, rep(NaN, 100))
#' for (i in 2:100) {
#'   S[i] <- genHeston(S[i-1], 0.1)
#' }
#' # heston(100, 0.05, 2, 0.4, -0.6, 0.4, 0.4, 100) returns the same output but faster!
getHestonGenerator  <- function(r, kappa, theta, rho, xi, sigmaZero, antithetic=FALSE, stepSize=FALSE) {

  # volatility process
  V <- FALSE
  
  return (function(S, dt, notStopped=NULL, reset=FALSE) {
    
    if(reset) {
      V <<- FALSE
      return()
    }
    
    if (is.list(S)) {
      numPaths <- length(S[[1]])
      numAssets <- length(S)
    } else {
      numPaths <- length(S)
      numAssets <- 1
    }
    
    if (identical(V, FALSE)) {
      if(antithetic) {
        initiaValue <- rep(sigmaZero, numPaths / 2)
      } else {
        initiaValue <- rep(sigmaZero, numPaths)
      }
      
      if (numAssets == 1) {
        V <<- initiaValue
      } else {
        V <<- list()
        for (d in 1:numAssets) {
          V[[d]] <<- initiaValue
        }
      }
    } 
    
    if (is.null(notStopped)) {
      notStopped <- rep(TRUE, numPaths)
    }
    
    if (identical(stepSize, FALSE)) {
      numSteps <- 1
    } else {
      if((dt%%stepSize) != 0) stop('stepSize ist kein Teiler von dt')
      numSteps <- as.integer(dt / stepSize)
      dt <- stepSize
    }
    
    numTotal <- length(notStopped)
    sqdt <- sqrt(dt)
    
    for (i in 1:numSteps) {
      for (d in 1:numAssets) {
        if (antithetic) {
          Z <- rnorm(numTotal/2)
          Y <- rnorm(numTotal/2)
          W <- rho * Z + sqrt(1 - rho^2) * Y
          Z <- c(Z, -Z)
          Z <- Z[notStopped]
        } else {
          Z <- rnorm(numTotal)
          Y <- rnorm(numTotal)
          W <- rho * Z + sqrt(1 - rho^2) * Y
          Z <- Z[notStopped]
        }
        
        if (numAssets == 1) {
          truncatedV <- pmax(V, 0)
          V <<- V + kappa * (theta- truncatedV ) * dt + xi * sqrt(truncatedV) * sqdt * W
          if(antithetic) {
            truncatedV <- c(truncatedV, truncatedV)
          }
          S <-  S * exp((r - 0.5 * truncatedV[notStopped]) * dt + sqrt(truncatedV[notStopped]) * sqdt * Z)
        } else {
          truncatedV <- pmax(V[[d]], 0)
          V[[d]]  <<- V[[d]]  + kappa * (theta- truncatedV ) * dt + xi * sqrt(truncatedV) * sqdt * W
          if(antithetic) {
            truncatedV <- c(truncatedV, truncatedV)
          }
          S[[d]] <-  S[[d]] * exp((r - 0.5 * truncatedV[notStopped]) * dt + sqrt(truncatedV[notStopped]) * sqdt * Z)
        }
      } 
    }    
    return(S)
  })
}