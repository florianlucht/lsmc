\name{getBullspread}
\alias{getBullspread}
\title{Generate Bull-Spread option}
\usage{
getBullspread(strike)
}
\arguments{
  \item{strike}{parameter values \code{K1} and \code{K2}
  (vector)}
}
\value{
function with arguments: \describe{ \item{x}{market
value(s) of the underlying(s) (scalar or vector)} }
}
\description{
This method returns a function which represents a
Bull-Spread option with payoff \deqn{\min\left\{K_2,
(x-K_1)^+\right\}}{min[ K2, max[(x - K1), 0] ]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getCall}};
\code{\link{getMinPut}}; \code{\link{getPut}};
\code{\link{minPut}}; \code{\link{put}}
}

