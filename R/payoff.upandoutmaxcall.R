####################################################################################################
# upAndOutMaxCall
####################################################################################################
#' Up-And-Out Max-Call option
#'
#' Option with payoff
#' \deqn{(\max\{S_t^{(1)},\dots,S_t^{(d)}\} - K )^+  \cdot 1_{\max\{S_u^{(d)} \mid 0\leq u \leq t \ , \ 1 \leq l \leq d \} \leq B }}{ max{ max{ S^1 , ... , S^d } - K , 0 }  * ( S_u^n < B) fuer alle 0 <= u <=t und alle S.}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x        market value(s) of the underlying(s)  (scalar or vector)
#' @param strike   strike (scalar)
#' @param barrier  barrier (scalar)
#' @param path     previous path of the underlying asset's price (vector)
#'
#' @return Optionswert(e)
upAndOutMaxCall <- function(x, strike, barrier, path) {

  numPath <- length(x[[1]])
  numAssets <- length(x)
  
  value <- rep(0, numPath)
  notKnockedOut <- rep(TRUE, numPath)

  if (is.list(path)) {
    for (d in 1:numAssets) {
      notKnockedOut <- notKnockedOut &  (apply(path[[d]], 2, max) < barrier)
    }
  } else {
    notKnockedOut <- (path < barrier)
  }

  if(any(notKnockedOut)) {
    x <- matrix(unlist(x), numPath, numAssets)
    x <- x[notKnockedOut, ]
    if(is.matrix(x)) {
      maxi <- apply(x, 1, max)
    } else {
      maxi <- max(x)
    }      
    value[notKnockedOut] <-  pmax(0, maxi - strike)
  }
      
  return(value)
}

####################################################################################################
# getUpAndOutMaxCall
####################################################################################################
#' Generate Up-And-Out Max-Call option
#'
#' This method returns a function which represents a Up-And-Out Max-Call option with payoff
#' \deqn{(K - \min_{i= 1,\dots,d} x_i)^+}{max[(K - min[x_i]), 0]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param strike   sktrike (scalar)
#' @param barrier  barrier (scalar)
#' @param forward  logical. If \code{TRUE} the function keeps the knock-out information.
#'
#' @return   Methode mit folgenden Arguementen:
#' \describe{
#'   \item{x}{Basiswert(e) (Vektor oder einzelner Wert)}
#'   \item{para}{TODO}
#'   \item{notStopped}{Index der nichtausgeuebten Pfade (siehe lsmc Funktion)}
#'   \item{reset}{falls TRUE wird beim naechsten Aufruf die Knocked-out Information verworfen}
#'  } 
getUpAndOutMaxCall <- function(strike, barrier, forward=FALSE) {
  
  strike <- strike
  notKnockedOut <- NULL 
  
  if (forward) {
    return( function(x, para, notStopped, reset=FALSE) {

      if (reset) {
        notKnockedOut <<- NULL 
        return()
      }

      numPath <- length(x[[1]])
      numAssets <- length(x)
      
      if(is.null(notKnockedOut)) {
        notKnockedOut <<- rep(TRUE, numPath)
      }
      
      value <- rep(0, numPath)
      
      x <- matrix(unlist(x), numPath, numAssets)
        if(is.matrix(x)) {
          maxi <- apply(x, 1, max)
        } else {
          maxi <- max(x)
        }     

      notKnockedOut[notStopped] <<-  notKnockedOut[notStopped] & (maxi < barrier)
      value[notKnockedOut[notStopped]] <-  pmax(0, maxi[notKnockedOut[notStopped]] - strike)
      
      return(value)
    })
    
  } else {
    f <- function(x, path) {
      return(upAndOutMaxCall(x, strike, barrier, path))
    }
    return(f)
  }
}