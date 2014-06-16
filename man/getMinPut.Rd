\name{getMinPut}
\alias{getMinPut}
\title{Generate Min-Put option}
\usage{
getMinPut(strike)
}
\arguments{
  \item{strike}{strike (scalar)}
}
\value{
function with arguments: \describe{ \item{x}{market
value(s) of the underlying(s) (scalar or vector)} }
}
\description{
This method returns a function which represents a Min-Put
option with payoff \deqn{(K - \min_{i= 1,\dots,d}
x_i)^+}{max[(K - min[x_i]), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getPut}};
\code{\link{getUpAndOutMaxCall}}; \code{\link{minPut}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

