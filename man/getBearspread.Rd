\name{getBearspread}
\alias{getBearspread}
\title{Generate Bear-Spread option}
\usage{
getBearspread(strike)
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
Bear-Spread option with payoff \deqn{\min\left\{K_2,
(K_1-x)^+\right\}}{min[ K2, max[(K1 - x), 0] ]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBullspread}}; \code{\link{getCall}};
\code{\link{getMinPut}}; \code{\link{getPut}};
\code{\link{getUpAndOutMaxCall}}; \code{\link{minPut}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

