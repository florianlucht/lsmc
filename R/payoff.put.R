####################################################################################################
# put
####################################################################################################
#' Put option
#'
#' Option with payoff
#' \deqn{(K-x)^+}{max[(K - x), 0]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x       market value(s) of the underlying(s)  (scalar or vector)
#' @param strike  strike (scalar)
#'
#' @return option value(s)
put <- function(x, strike) {
  pmax(0, strike-x)
}

####################################################################################################
# getPut
####################################################################################################
#' Generate Put option
#'
#' This method returns a function which represents a Put option with payoff
#' \deqn{(K-x)^+}{max[(K-x), 0]}
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
getPut <- function(strike) {
  return (function(x) {
    pmax(0, strike - x)
  })
}