\name{bearspread}
\alias{bearspread}
\title{Bear-Spread option}
\usage{
bearspread(x, strike)
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
(K_1-x)^+\right\}}{min[ K2, max[(K1 - x), 0] ]}
}
\seealso{
Other payoff: \code{\link{bullspread}}; \code{\link{call}};
\code{\link{getBearspread}}; \code{\link{getBullspread}};
\code{\link{getCall}}; \code{\link{getMinPut}};
\code{\link{getPut}}; \code{\link{minPut}};
\code{\link{put}}
}

