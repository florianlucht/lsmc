####################################################################################################
# merton
####################################################################################################
#' Merton model
#'
#' This method simulates discrete paths in the Merton model.
#' The price paths satisfy the following stochastic differential equations: 
#' \deqn{S_t = S_0 \exp \left[ \left( \mu-\frac{\sigma^2}{2} \right) t + \sigma W_t
#'        + \sum\limits_{i=1}^{N_t} Y_i \right] \ ,}
#'      {S_t = S_0 exp[ ( mu- 1/2 * sigma^2 ) * t + sigma * W_t + sum[ Y_i | i = 1,...,N ] ] }
#' \deqn{N_t \sim Pois(\lambda t) \ , \ Y_i \sim N(\mu, \delta) \ . }
#'      {N ~ Pois(lambda*t) , Y_i ~ N(mu, delta) . }
#' 
#' @export
#'
#' @family models
#' 
#' @param initialValue initial value (scalar)
#' @param r            interest rate (scalar)
#' @param sigma        volatility (scalar)
#' @param lambda       intensity/rate of the jumps (scalar)
#' @param mu           mean of the jumps (scalar)
#' @param delta        variance of the jumps (scalar)
#' @param numDates     number of dates (scalar)
#' @param maturity     maturity (scalar)
#' @param numPaths     number of paths (scalar)
#' @param antitheric   logical. If \code{TRUE} antitethic paths are generated.
#' @param plot         logical. If \code{TRUE} paths are plotted.
#'
#' @return simulated paths as matrix with \code{numDates+1} rows and \code{numPaths} columns 
#'
#' @examples
#' # generate one simple path
#' S <- merton(100, 0.05, 0.4, 1, 0, 0.4, 100)
#' 
#' # generate 100 paths with maturity T=10
#' S <- merton(100, 0.05, 0.4, 1, 0, 0.4, 100, maturity=10, numPaths=100)
merton <- function(initialValue, r, sigma, lambda, mu, delta, numDates,
                   maturity=1, 
                   numPaths=1, 
                   antithetic=FALSE, 
                   plot=FALSE) {
  # step size
  dt <- maturity/numDates
  
  kappa <- exp(mu + delta^2 / 2) - 1
 
  # paths
  if (antithetic) {
    Z <- rnorm(numDates * numPaths / 2)
    Z <- c(Z, -Z)
    J <- rpois(numDates * numPaths / 2,lambda * dt)
    J <- c(J, J)
    Y <- rnorm(numDates * numPaths / 2)
    J <- mu * J + delta * sqrt(J) * c(Y, -Y)
    rm(Y)
  } else  {
    Z <- rnorm(numDates * numPaths)
    J <- rpois(numDates * numPaths,lambda * dt)
    J <- mu * J + delta * sqrt(J) * rnorm(numDates * numPaths)
  }
  S <- matrix(exp((r - 0.5 * sigma^2 - lambda * kappa) * dt + sigma * sqrt(dt) * Z + J), numDates, numPaths)   
  S <- rbind(rep(initialValue, numPaths), S)

  S <- apply(S, 2, cumprod)  
  
  # verbose
  if (plot) {  
    matplot(seq(0, maturity, dt), S, type="l", xlab="Time")
  }
  
  return(S)	
}

####################################################################################################
# getMertonGenerator
####################################################################################################
#' Generator function of the Merton model 
#' 
#' This method returns a function which can be used to simulate one step in the Merton model
#' for given parameters.
#'
#' @export
#'
#' @param r            interest rate (scalar)
#' @param sigma        volatility (scalar)
#' @param lambda       intensity/rate of the jumps (scalar)
#' @param mu           mean of the jumps (scalar)
#' @param delta        variance of the jumps (scalar)
#' @param antithetic   logical. If \code{TRUE} antitethic paths are generated. 
#'
#' @return   function with arguments:
#' \describe{
#'   \item{S}{current value(s) (scalar or vector)}
#'   \item{dt}{step size (scalar)}
#'   \item{notStopped}{optional index of paths which are not exercised (vector of logical). This 
#'                     argument is used internally to calculate the correct values if 
#'                     \code{antithetic=TRUE} is given.}
#'  }    
#'
#' @examples
#' # get generator function
#' genMerton <- getMertonGenerator(0.05, 0.4, 1, 0, 0.4)
#' # generate one step in the Merton model with step size 0.1 and current value 100
#' next <- genMerton(100, 0.1)
#' # reset volatility process
#' genMerton(NULL, NULL, reset=TRUE)
#' # generate path with 100 steps 
#' S <- c(100, rep(NaN, 100))
#' for (i in 2:100) {
#'   S[i] <- genMerton(S[i-1], 0.1)
#' }
#' # merton(100,0.05, 0.4, 1, 0, 0.4, 100) returns the same output but faster!
getMertonGenerator  <- function(r, sigma, lambda, mu, delta, antithetic=FALSE) {
  return( function(S, dt, notStopped=NULL) { 
    
    if (is.list(S)) {
      numPaths <- length(S[[1]])
      numAssets <- length(S)
    } else {
      numPaths <- length(S)
      numAssets <- 1
    }
    
    if (is.null(notStopped)) {
      notStopped <- rep(TRUE, numPaths)
    }
    
    kappa <- exp(mu + delta^2 / 2) - 1
    
    for(d in 1:numAssets) {
      if (antithetic) {
        numTotal <- length(notStopped) / 2
        Z <- rnorm(numTotal)
        Z <- c(Z, -Z)
        Z <- Z[notStopped]
        N <- rpois(numTotal, lambda * dt)
        N <- c(N, N)
        Y <- rnorm(numTotal)
        J <- mu * N + delta * sqrt(N) * c(Y, -Y)
        J <- J[notStopped]
      } else {
        Z <- rnorm(numPaths)
        N <- rpois(numPaths,lambda * dt)
        J <- mu * N + delta * sqrt(N) * rnorm(numPaths)
      }
      
      if (numAssets==1) {
        S <- S * exp( (r - 0.5 * sigma^2 - lambda * kappa) * dt + sigma * sqrt(dt) * Z + J)
      } else {
        S[[d]] <- S[[d]] * exp( (r - 0.5 * sigma^2 - lambda * kappa) * dt + sigma * sqrt(dt) * Z + J)
      }
    } 
    return(S)
  })
}
