\name{getPut}
\alias{getPut}
\title{Generate Put option}
\usage{
getPut(strike)
}
\arguments{
  \item{strike}{strike (scalar)}
}
\value{
function with arguments: \describe{ \item{x}{market
value(s) of the underlying(s) (scalar or vector)} }
}
\description{
This method returns a function which represents a Put
option with payoff \deqn{(K-x)^+}{max[(K-x), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getUpAndOutMaxCall}}; \code{\link{minPut}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

