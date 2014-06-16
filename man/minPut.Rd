\name{minPut}
\alias{minPut}
\title{Min-Put option}
\usage{
minPut(x, strike)
}
\arguments{
  \item{x}{market value(s) of the underlying(s) (scalar or
  vector)}

  \item{strike}{strike (scalar)}
}
\value{
Optionswert(e)
}
\description{
Option with payoff \deqn{(K - \min_{i= 1,\dots,d}
x_i)^+}{max[(K - min[x_i]), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{getUpAndOutMaxCall}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

