\name{call}
\alias{call}
\title{Call option}
\usage{
call(x, strike)
}
\arguments{
  \item{x}{market value(s) of the underlying(s) (scalar or
  vector)}

  \item{strike}{strike (scalar)}
}
\value{
option value(s)
}
\description{
Option with payoff \deqn{(x-K)^+}{max[(x - K), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{getBearspread}};
\code{\link{getBullspread}}; \code{\link{getCall}};
\code{\link{getMinPut}}; \code{\link{getPut}};
\code{\link{getUpAndOutMaxCall}}; \code{\link{minPut}};
\code{\link{put}}; \code{\link{upAndOutMaxCall}}
}

