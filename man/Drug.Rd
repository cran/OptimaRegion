\name{Drug}
\alias{Drug}
\docType{data}
\title{
Mixture-amount experiment
}
\description{
A pharmaceutical mixture-amount experiment in two components
}
\usage{data("Drug")}
\format{
  A data frame with 360 observations on the following 3 variables.
  \describe{
    \item{\samp{Component_1}}{Component 1 amount (mg)}
    \item{\samp{Component_2}}{Component 2 amount (mg)}
    \item{\code{Percent}}{Percent of cells killed (response)}
  }
}


\examples{

plot(Drug[,1:2])
}
\keyword{datasets}
