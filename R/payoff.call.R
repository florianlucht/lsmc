####################################################################################################
# call
####################################################################################################
#' Call option
#'
#' Option with payoff
#' \deqn{(x-K)^+}{max[(x - K), 0]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x       market value(s) of the underlying(s)  (scalar or vector)
#' @param strike  strike (scalar)
#'
#' @return option value(s)
call <- function(x, strike) {
  pmax(0, x - strike)
}


####################################################################################################
# getCall
####################################################################################################
#' Generate Call option
#'
#' This method returns a function which represents a Call option with payoff
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
getCall <- function(strike) {
  return (function(x){
    pmax(0, x - strike)
  })
}