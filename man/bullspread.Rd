\name{bullspread}
\alias{bullspread}
\title{Bull-Spread option}
\usage{
bullspread(x, strike)
}
\arguments{
  \item{x}{market value(s) of the underlying(s) (scalar or
  vector)}

  \item{strike}{parameter values \code{K1} and \code{K2}
  (vector)}
}
\value{
option value(s)
}
\description{
Option with payoff \deqn{\min\left\{K_2,
(x-K_1)^+\right\}}{min[ K2, max[(x - K1), 0] ]}
}
\seealso{
Other payoff: \code{\link{bearspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{minPut}};
\code{\link{put}}
}

