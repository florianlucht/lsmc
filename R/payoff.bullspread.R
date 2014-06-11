####################################################################################################
# bullspread
####################################################################################################
#' Bull-Spread option
#'
#' Option with payoff
#' \deqn{\min\left\{K_2, (x-K_1)^+\right\}}{min[ K2, max[(x - K1), 0] ]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param x       market value(s) of the underlying(s)  (scalar or vector)
#' @param strike  parameter values \code{K1} and \code{K2} (vector)
#'
#' @return option value(s)
bullspread <- function(x, strike)
{
  min(strike[2], max(0, x - strike[1]))
}


####################################################################################################
# getBullspread
####################################################################################################
#' Generate Bull-Spread option
#'
#' This method returns a function which represents a Bull-Spread option with payoff
#' \deqn{\min\left\{K_2, (x-K_1)^+\right\}}{min[ K2, max[(x - K1), 0] ]}
#'
#' @export
#'
#' @family payoff
#' 
#' @param strike  parameter values \code{K1} and \code{K2} (vector)
#'
#' @return function with arguments:
#' \describe{
#'   \item{x}{market value(s) of the underlying(s)  (scalar or vector)}
#'  } 
getBullspread <- function(strike) {
  K1 <- strike[1]
  K2 <- strike[2]
  return (function (x) {
  	 pmin(K2, pmax(0, x - K1))
  })
}