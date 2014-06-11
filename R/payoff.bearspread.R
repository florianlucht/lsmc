####################################################################################################
# bearspread
####################################################################################################
#' Bear-Spread option
#'
#' Option with payoff
#' \deqn{\min\left\{K_2, (K_1-x)^+\right\}}{min[ K2, max[(K1 - x), 0] ]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x       market value(s) of the underlying(s)  (scalar or vector)
#' @param strike  parameter values \code{K1} and \code{K2} (vector)
#'
#' @return option value(s)
bearspread <- function(x, strike) {
  pmin(strike[2],pmax(0,strike[1] - x))
}


####################################################################################################
# getBearspread
####################################################################################################
#' Generate Bear-Spread option
#'
#' This method returns a function which represents a Bear-Spread option with payoff
#' \deqn{\min\left\{K_2, (K_1-x)^+\right\}}{min[ K2, max[(K1 - x), 0] ]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param strike parameter values \code{K1} and \code{K2} (vector)
#'
#' @return function with arguments:
#' \describe{
#'   \item{x}{market value(s) of the underlying(s)  (scalar or vector)}
#' }    
getBearspread <- function(strike) {
  K1 <- strike[1]
  K2 <- strike[2]
  return (function (x) {
  	pmin(K2, max(0, K1 - x))
  })
}