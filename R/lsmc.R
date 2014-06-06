####################################################################################################
# lsmc
####################################################################################################
#' Longstaff-Schwartz Algorithm
#'
#' \coder(lsmc) implements the least-squares monte-carlo alogrithm of Longstaff and Schwartz to
#' valuing american options.
#'
#'
#' @param X                  matrix where each column represents a price path of the underlying.
#'                           In the case of multiple assets, X is a list with length of number of 
#'                           assets. Each element corresponds to the paths of each asset as matrix.
#'                           If \code{pathGeneratorFct} is specified, \code{X} (scalar) will be used
#'                           as initial value. 
#' @param payoffFct          payoff function with arguments \code{(S, ...)} where \code{S} is a 
#'                           vector of prices. In the case of multiple assets, \code{S} is a list 
#'                           of vectors.
#' @param interestRate       interest rate (scalar)
#' @param maturity           maturity (scalar)
#' @param payoffParameter    optional list of arguments for the payoff function (list) 
#' @param numExerciseDates   optional number of exercise dates (scalar). If the value is not 
#'                           specified, the number of rows of X minus one is used. 
#' @param basisFct           optional basis function used for regression (function)
#' @param numBasisFct        number of basis functions (scalar). Required if \code{basisFct} is used.
#' @param regressionCoef     optional precalculated regression coefficients
#'                           (\code{numExerciseDates-1} x \code{numBasisFct} matrix)
#'                           (leads to low biased valuation)
#' @param regressionMethod   optional string to specify the method for solving the least squares
#'                           problem. Currently \code{method = "qr"} for QR decomposition (default)
#'                           and \code{method = "svd"} for singular value decomposition are supported.
#' @param onlyInTheMoney     logical. If \code{TRUE} (default) only paths which are in-the-money are
##'                          used.
#' @param pathGeneratorFct   optional function to generate price paths. If this is specified, 
#'                           \code{X} is used as initial value.
#' @param pathInfo           an optional additional path information needed for valuation, e.g.
#'                           average price for asian option (scalar or vector)
#' @param verbose            logical. If \code{TRUE} extra information will be printed.
#'
#' @return list:
#' \describe{
#'   \item{V}{vector of option values}
#'   \item{regressionCoef}{regression coefficients as \code{numExerciseDates}x\code{numBasisFct}
#'         matrix}
#'   \item{exerciseDates}{vector of exercise dates for each path}
#'   \item{elapsedTime}{elapsed time for calculation}
#' }    
#' 
#' @examples  
#' T <- 1
#' r <- 0.05
#' sigma <- 0.4
#' initalValue <- 100            
#'
#' N <- 100
#' M <- 10000 
#'
#' S <- blackscholes(initialValue, r, sigma, numDates=N, maturity=T, numPaths=M)
#'
#' strike <- 100
#' payoffParameter <- list(strike=strike)
#' payoff <- put 
#'  
#' result <- lsmc(S, payoff, r, T, payoffParameter)
#' cat("option value: ", result$V)
lsmc <- function(X, payoffFct, interestRate, maturity,
                 payoffParameter  = NULL,
                 numExerciseDates = NULL, 
                 basisFct         = NULL,  
                 numBasisFct      = NULL,
                 regressionCoef   = NULL, 
                 regressionMethod = 'qr', 
                 onlyInTheMoney   = TRUE,
                 pathGeneratorFct = NULL,
                 pathInfo         = NULL, 
                 verbose          = FALSE) {
  # environment for private methods
  environment(.lsmc.assets) <- environment()
  environment(.lsmc.payoff) <- environment()
  environment(.lsmc.backward) <- environment()
  environment(.lsmc.forward) <- environment()
  environment(.lsmc.inputValidation) <- environment()

  # validate arguments
  .lsmc.inputValidation()

  # check for multiple asset and on-the-fly path generation
  multiAsset <- is.list(X)
  onTheFly <- !is.null(pathGeneratorFct)

  # number of paths M
  if (multiAsset) {
    numOfAssets <- length(X)
    M <- ifelse(onTheFly, length(X[[1]]), dim(X[[1]])[2])
  } else {
    M <- ifelse(onTheFly, length(X), dim(X)[2])
  } 

  # number of exercise dates N and
  # number of steps per execise date dN
  if (is.null(numExerciseDates)) {
    N <- ifelse(multiAsset, dim(X[[1]])[1] , dim(X)[1]) - 1
    dN <- 1
  } else {
    N <- numExerciseDates
    if (onTheFly) {
      dN <- 1
    } else {
      dN <- (ifelse(multiAsset, dim(X[[1]])[1] , dim(X)[1]) - 1) %/% N 
    }
  }

  # standard basis functions
  if(is.null(basisFct)) {
    if (multiAsset) {
      basisFct <- .lsmc.basisFctMultiAsset
      numBasisFct <- 2 * numOfAssets + 2
    } else {
      basisFct <- .lsmc.basisFct
      numBasisFct <- 4
    }
  }
  
  # ermittelte Ausuebungszeitpunkte
  exerciseDates <- rep(NA,M)
  
  # verbose
  if (verbose) {
    cat('\nM', M )
    cat('\nN', N )
    cat('\nT', maturity )
    cat('\nr', interestRate )  
    cat('\n\ncurrent step: ')
  }  
  
  #-----------------------------------------------------------------------------
  # calculation 
  #-----------------------------------------------------------------------------
  # start time measuring 
  ptm <- proc.time()

  # european option
  if (N == 1) {
    currentDate <- N * dN + 1
    S <- .lsmc.assets(currentDate)
    V <- exp(-interestRate * maturity) * .lsmc.payoff(S, currentDate)
    exerciseDates[V>0] <- currentDate
    return(list(V = mean(V), regressionCoef=NULL, exerciseDates=exerciseDates, 
                elapsedTime = (proc.time() - ptm)[1] ))
  } 

  # if no regression coefficients are avaiable ... 
  if (is.null(regressionCoef)) {  
    result <- .lsmc.backward()
    regressionCoef <- result$regressionCoef
    exerciseDates <- result$exerciseDates
  # else (low biased valuation)
  } else {
    result <- .lsmc.forward()
    exerciseDates <- result$exerciseDates
  }
  
  # option value
  V <- mean(result$V)

  # stop time measuring 
  elapsedTime <- (proc.time() - ptm)[1]
  
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  if (verbose) {
    cat('\n\nV        ',V)
    cat(  '\nLaufzeit ', elapsedTime, "s")
  }
  
  return(list(V = V, 
              regressionCoef = regressionCoef, 
              exerciseDates = exerciseDates, 
              elapsedTime = elapsedTime ))
}


####################################################################################################
# .lsmc.backward
####################################################################################################
.lsmc.backward <- function() {

  regressionCoef <- matrix(NA, N-1, numBasisFct)

  # index of used paths
  itm <- rep(TRUE, M)
  
  currentDate <- N * dN + 1
  S <- .lsmc.assets(currentDate)
  V <- exp(-interestRate * maturity) * .lsmc.payoff(S, currentDate)
  exerciseDates[V>0] <- currentDate

  for(n in seq(N-1, 1))
  {
    if (verbose) cat(n, '')
    
    currentDate <- n*dN+1
      
    S <- .lsmc.assets(currentDate)
    P <- .lsmc.payoff(S, currentDate)
    
    if (onlyInTheMoney == TRUE) itm <- (P > 0)
    
    if (sum(itm) > numBasisFct)
    {
      H <- basisFct(S, P)
      
      if (regressionMethod == "qr") {
        regressionCoef[n, ] <- tryCatch({
          solve( t(H[itm, ]) %*% H[itm, ] ) %*% t(H[itm, ]) %*% V[itm]
        }, error=function(e) { 
          warning('Step ',n,': ',e) 
          rep(NA,numBasisFct)
        })        
      } else if (regressionMethod == "svd") {
        SVD <- svd(H[itm, ])          
        regressionCoef[n, ]   = SVD$v %*% diag(1/SVD$d) %*% t(SVD$u) %*% V[itm]
      }
      else 
        stop("method not known") 
      
      if (!any(is.na(regressionCoef[n, ]))) {

        # continuation value
        C <- H %*% regressionCoef[n, ]  
        
        P <- exp(-(n * maturity / N) * interestRate) * P
        
        stopped <- (P >= C & P > 0)

        V[stopped] <- P[stopped]
        exerciseDates[stopped] <- currentDate
      }
    }      
  }  
  return(list(V = V, regressionCoef = regressionCoef, exerciseDates = exerciseDates))
}


####################################################################################################
# .lsmc.forward
####################################################################################################
.lsmc.forward <- function() {
  # index of paths wich are not exercised 
  notStopped <- rep(TRUE, M)
  
  V <- rep(0, M)

  if(onTheFly) {
    dt <- maturity/N
    S <- X
    #stopped <- rep(FALSE, length(S))
  }

  if (any(names(formals(payoffFct)) == 'reset'))
    payoffFct(reset=TRUE)

  for (n in seq(1, N-1)) { 
    if (verbose) cat(n, '')

    currentDate <- n*dN+1
    
    # current price
    if(onTheFly) {
      if (any(names(formals(pathGeneratorFct)) == 'notStopped')) {
        S <- pathGeneratorFct(S, dt, notStopped=notStopped)
      } else {
        S <- pathGeneratorFct(S, dt)
      }
    } else {
      S <- .lsmc.assets(currentDate, notStopped)
    }

    # intrinsic value 
    P <- .lsmc.payoff(S, currentDate, notStopped)

    if (!any(is.na(regressionCoef[n, ]))) {
      
      # contination value
      H <- basisFct(S,P)
      C <- H %*% regressionCoef[n, ]   
      
      P <- exp(-(n*maturity/N) * interestRate) * P
      
      # index of paths wich were currently exercised
      stopped <- ( P >= C & P > 0 )



      # current value
      V[notStopped][stopped] <- P[stopped]
      exerciseDates[notStopped][stopped] <- currentDate 
      notStopped[notStopped][stopped] <- FALSE
      if (onTheFly) {
        S <- .lsmc.util.getSub(S, !stopped)
      }
    }
    
    # stop if all paths have been exercised
    if (sum(notStopped)==0) break
  }

  # exercise in the last step 
  currentDate <- N * dN + 1
  if (onTheFly) {
    if (any(names(formals(pathGeneratorFct)) == 'notStopped')) {
      S <- pathGeneratorFct(S, dt, notStopped=notStopped)
    } else {
      S <- pathGeneratorFct(S, dt)
    }
  } else {
    S <- .lsmc.assets(currentDate, notStopped)
  }
  P <- .lsmc.payoff(S, currentDate, notStopped)
  stopped <- (P > 0)
  V[notStopped][stopped] <- exp(-interestRate * maturity) * P[stopped]
  exerciseDates[notStopped][stopped] <- currentDate

  return(list(V = V, exerciseDates = exerciseDates))
}


####################################################################################################
# .lsmc.assets
####################################################################################################
.lsmc.assets <- function(currentDate, notStopped=NULL) {
  .lsmc.util.getSub(X, currentDate, notStopped)
}


####################################################################################################
# .lsmc.payoff
####################################################################################################
.lsmc.payoff <- function(S, currentDate, notStopped=NULL) {
  args <- list()
  args$x <- S

  if (is.list(payoffParameter)) {
    name <- names(payoffParameter)
    for(i in 1:length(payoffParameter)) {
      args[[name[i]]] <- payoffParameter[[i]]
    }
  }

  if (any(names(formals(payoffFct)) == 'path')) {
    if (is.null(pathInfo)) {
      args$path <- .lsmc.util.getSub(X, 1:currentDate, notStopped)
    } else {
      args$path <- .lsmc.util.getSub(pathInfo, currentDate, notStopped)
    }
  }

  if (any(names(formals(payoffFct)) == 'notStopped')) {
    args$notStopped <- notStopped
  }

  return(do.call(payoffFct, args))
}


####################################################################################################
# .lsmc.basisFct
####################################################################################################
.lsmc.basisFct <- function(S, P) {
  return(matrix(c( rep(1,length(P)), P, P^2, P^3), length(P), 4))
}


####################################################################################################
# .lsmc.basisFctMultiAsset
####################################################################################################
.lsmc.basisFctMultiAsset <- function(S, P) {
  numOfAssets <- length(S)
  numOfPath <- length(S[[1]])

  H <- matrix(1,numOfPath,2*numOfAssets+2)
  
  tempData <- matrix(unlist(S), numOfPath, numOfAssets)
 
  H[, seq(2, numOfAssets+1)] <- tempData
  H[, seq(numOfAssets+2, 2*numOfAssets+1)] <- tempData*tempData      
  H[, 2*numOfAssets+2] <- P

  return(H)
}


####################################################################################################
# .lsmc.inputValidation
####################################################################################################
.lsmc.inputValidation <-function() {
  # Argument X
  if( ! (is.list(X) | is.vector(X) | is.matrix(X)))
    stop("Invalid Argument: X wrong type") 

  # TODO not yet implemented
} 


####################################################################################################
# .lsmc.util.getSub
####################################################################################################
.lsmc.util.getSub <- function(x,row=TRUE,col=TRUE) {
  if (is.null(row)) row <- TRUE
  if (is.null(col)) col <- TRUE

  if (is.list(x)) {
    val <- list()
    for (d in 1:length(x)) {
      if (is.vector(x[[d]])) {
        val[[d]] <- x[[d]][row]
      } else {
        val[[d]] <- x[[d]][row, col]
      }
    }
    return(val)
  } else if (is.vector(x)) {
    return(x[row])
  } else {
    return(x[row,col])
  }
}