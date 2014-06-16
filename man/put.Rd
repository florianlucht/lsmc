\name{put}
\alias{put}
\title{Put option}
\usage{
put(x, strike)
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
Option with payoff \deqn{(K-x)^+}{max[(K - x), 0]}
}
\seealso{
Other payoff: \code{\link{bearspread}};
\code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{minPut}}
}

