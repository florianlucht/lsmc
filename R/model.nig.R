####################################################################################################
# nig
####################################################################################################
#' Normal-inverse Gaussian model
#'
#' This method simulates discrete paths in the normal-inverse Gaussian model.
#' Dieses ist beschrieben durch
#' \deqn{ S_t = S_0 \exp\left[ (r-\omega) + \beta I_t + W_{I_t} \right] \ , }{S_t = S_0 exp[ (r - omega) + beta * I_t + W_{I_t} ] ,} 
#' \deqn{ I_t \sim IG(\delta t, \gamma) \ , \  \omega = \lambda \gamma - \lambda \sqrt{\alpha^2  (1+\beta)^2} \ .}{I_t ~ IG(delta*t , gamma) , omega = lambda * gamma - lambda * sqrt{ alpha^2 * (1 + beta)^2 } .} 
#'
#' @export
#'
#' @family models
#' 
#' @param initialValue initial value (scalar)
#' @param r            interest rate (scalar)
#' @param alpha        tail heavyness (scalar)
#' @param beta         asymmetry parameter (scalar)
#' @param delta        scale parameter (scalar)
#' @param numDates     number of dates (scalar)
#' @param maturity     maturity (scalar)
#' @param numPaths     number of paths (scalar)
#' @param antithetic   logical. If \code{TRUE} antitethic paths are generated.
#' @param plot         logical. If \code{TRUE} paths are plotted.
#'
#' @return Pfade als Matrix mit \eqn{numDates + 1}{numDates + 1} Zeilen und \eqn{numPaths} Spalten
#' 
#' @examples
#' # generate one simple path
#' S <- nig(100, 0.05, 20, -5, 1, 100)
#' 
#' # generate 100 paths with maturity T=10
#' S <- nig(100, 0.05, 20, -5, 1, 100, maturity=10, numPaths=100)
nig <- function(initialValue, r, alpha, beta, delta, numDates, 
                maturity=1, 
                numPaths=1, 
                plot=FALSE) {
  dt = maturity / numDates
  a = 1
  b = delta * sqrt(alpha^2 - beta^2)
  omega = delta * (sqrt(alpha^2 - (beta + 1)^2) - sqrt(alpha^2 - beta^2))
  
  S = matrix(.rinvg(numDates * numPaths, a * dt, b), numDates, numPaths)
  Z = matrix( rnorm(numDates * numPaths), numDates, numPaths)
  
  S = rbind(initialValue, exp((r + omega) * dt  + beta * delta^2 * S + delta * sqrt(S) * Z ))
  
  S = apply(S, 2, cumprod)
  
  # Ausgabe
  if(plot) {
    matplot(seq(0, maturity, maturity / numDates), S, type="l", xlab="Time")
  }
  
  return(S)
}

####################################################################################################
# getNigGenerator
####################################################################################################
#' Generator function of the normal-inverse Gaussian model
#' 
#' This method returns a function which can be used to simulate one step in the normal-inverse 
#' Gaussian model for given parameters.
#'
#' @export
#'
#' @param r             interest rate (scalar)
#' @param alpha         tail heavyness (scalar)
#' @param beta          asymmetry parameter (scalar) 
#' @param delta         scale parameter (scalar)
#'
#' @return function with arguments:
#' \describe{
#'   \item{S}{current value(s) (scalar or vector)}
#'   \item{dt}{step size (scalar)}
#'  }    
#'
#' @examples
#' # get generator function
#' genNig <- getNigGenerator(0.05, 20, -5, 1)
#' # generate one step in the NIG-model with step size 0.1 and current value 100
#' next <- genNig(100, 0.1)
#' # generate path with 100 steps 
#' S <- c(100, rep(NaN, 100))
#' for (i in 2:100) {
#'   S[i] <- genNig(S[i-1], 0.1)
#' }
#' # nig(100, 0.05, 20, -5, 1, 100) returns the same output but faster!
getNigGenerator  <- function(r, alpha, beta, delta) {
  
  a = 1
  b = delta * sqrt(alpha^2 - beta^2)
  omega = delta * (sqrt(alpha^2 - (beta + 1)^2) - sqrt(alpha^2 - beta^2))
  
  antithetic = FALSE
  
  return( function(S, dt) { 
    
    if (is.list(S)) {
      numPaths <- length(S[[1]])
      numAssets <- length(S)
    } else {
      numPaths <- length(S)
      numAssets <- 1
    }
      
    for(d in 1:numAssets) {
      if (antithetic) {
        Z <- rnorm(numPaths / 2)
        Z <- c(Z, -Z)
      } else {
        Z <- rnorm(numPaths)
        X = .rinvg(numPaths, a * dt, b)
      }
      
      if (numAssets == 1) {
        S <- S * exp((r + omega) * dt  + beta * delta^2 * X + delta * sqrt(X) * Z )
      } else {
        S[[d]] <- S[[d]] * exp((r + omega) * dt  + beta * delta^2 * X + delta * sqrt(X) * Z )
      }
    } 
    return(S)
  })
}

# randon number generator for the normal-inverse Gaussian distribution
.rinvg <- function(n,a,b)
{
  x = rnorm(n)
  x = x*x
  x = a/b + x/(2*b^2) - sqrt(4*a*b*x + x*x)/(2*b^2) 
  u = runif(n) 
  t = (u<=a/(a+x*b))  
  x = x*t + a^2/(b^2*x)*(!t) 
  return(x)  
}