####################################################################################################
# minPut
####################################################################################################
#' Min-Put option
#'
#' Option with payoff
#' \deqn{(K - \min_{i= 1,\dots,d} x_i)^+}{max[(K - min[x_i]), 0]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x       market value(s) of the underlying(s)  (scalar or vector)
#' @param strike  strike (scalar)
#'
#' @return Optionswert(e)
minPut <- function(x, strike) {
  numPath <- length(x[[1]])
  numAssets <- length(x)
  x <- matrix(unlist(x), numPath, numAssets)
  mini <- apply(x, 1, min)
  
  pmax(0, strike-mini)
}

####################################################################################################
# getMinPut
####################################################################################################
#' Generate Min-Put option
#'
#' This method returns a function which represents a Min-Put option with payoff
#' \deqn{(K - \min_{i= 1,\dots,d} x_i)^+}{max[(K - min[x_i]), 0]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param strike  strike (scalar)
#'
#' @return function with arguments:
#' \describe{
#'   \item{x}{market value(s) of the underlying(s)  (scalar or vector)}
#' }  
getMinPut <- function(strike) {
  return (minPut <- function(x) {
    numPath <- length(x[[1]])
    numAssets <- length(x)
    x <- matrix(unlist(x), numPath, numAssets)
    mini <- apply(x, 1, min)
  
    pmax(0, strike - mini) })
}
