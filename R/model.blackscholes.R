####################################################################################################
# blackscholes
####################################################################################################
#' Black-Scholes model
#'
#' This method simulates discrete paths in the Black-Scholes model.
#' The price paths follow a geometric Brownian motion, hence they satisfy the following stochastic 
#' differential equation
#' \deqn{dS  = S ( r dt + \sigma dW )}{S = S * (r * dt + sigma * dW) } 
#' i.e. its explicit solution
#' \deqn{S_t  = S_0 \exp\left[ \left(r - \frac{1}{2}\sigma^2 \right)t + \sigma W_t \right] }
#' {S_t  = S_0 exp[ (r - 1/2 * sigma^2) * t + sigma * W_t ] .}
#'
#' @export
#'
#' @family models
#' 
#' @param initialValue initial value (scalar)
#' @param r            interest rate (scalar)
#' @param sigma        volatility (scalar)
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
#' S <- blackscholes(100, 0.05, 0.4, 100)
#' 
#' # generate 100 paths with maturity T=10
#' S <- blackscholes(100, 0.05, 0.4, 100, maturity=10, numPaths=100)
blackscholes <- function(initialValue, r, sigma, numDates, 
                         maturity   = 1, 
                         numPaths   = 1, 
                         antithetic = FALSE, 
                         plot       = FALSE) {
	# step size	
  dt <- maturity / numDates 

  # paths 
  if (antithetic) {
    Z <- rnorm(numDates * numPaths / 2)
    Z <- c(Z, -Z)
  } else  {
    Z <- rnorm(numDates * numPaths)
  }
  S <- matrix(exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z), numDates, 
              numPaths)   
  S <- rbind(rep(initialValue, numPaths), S)
  S <- apply(S, 2, cumprod)

  # verbose
  if (plot) {
    matplot(seq(0, maturity, dt), S, type="l", xlab="Time", ylab="S") 
  }
  return(S)	  
}


####################################################################################################
# getBlackscholesGenerator
####################################################################################################
#' Generator function of the Black-Scholes model 
#' 
#' This method returns a function which can be used to simulate one step in the Black-Scholes model
#' for given parameters.
#' 
#' @export
#'
#' @param r            interest rate (scalar)
#' @param sigma        volatility (scalar)
#' @param antithetic   logical. If \code{TRUE} antitethic paths are generated.
#'
#' @return function with arguments:
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
#' genBlackscholes <- getBlackScholesGenerator(0.05, 0.4)
#' # generate one step in the BS-model with step size 0.1 and current value 100
#' next <- genBlackscholes(100, 0.1)
#' # genrate path with 100 steps 
#' S <- c(100, rep(NaN, 100))
#' for (i in 2:100) {
#'   S[i] <- genBlackscholes(S[i-1], 0.1)
#' }
#' # blackscholes(100, 0.05, 0.4, 100) returns the same output but faster!
getBlackscholesGenerator <- function(r, sigma, antithetic=FALSE) {
  return( function(S, dt, notStopped=NULL) { 
    
    # single or multiple assets?
    if (is.list(S)) {
      numPaths <- length(S[[1]])
      numAssets <- length(S)
    } else {
      numPaths <- length(S)
      numAssets <- 1
    }
    
    # for the correct calculation in the case of antithetic paths the index of 
    # the unexercised paths is required
    if (is.null(notStopped)) {
      notStopped <- rep(TRUE, numPaths)
    }

    # iteration over the assets
    for(d in 1:numAssets) {
      if (antithetic) {
        Z <- rnorm(length(notStopped) / 2)
        Z <- c(Z, -Z)
        Z <- Z[notStopped]
      } else {
        Z <- rnorm(numPaths)
      }
      
      if(numAssets == 1) {
        S <- S * exp( (r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      } else {
        S[[d]] <- S[[d]] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      }
    }
    return(S)
  })
}