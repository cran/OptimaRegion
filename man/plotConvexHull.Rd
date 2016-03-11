\name{plotConvexHull}
\alias{plotConvexHull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes and displays the convex hull of a set of 2-dimensional points
}
\description{
Given a vector of 2-dimensional coordinates, computes and displays the convex hull formed by these points.
}
\usage{
plotConvexHull(xin, LB, UB, xlab, ylab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xin}{
n x 2 vector of coordinate points
}
  \item{LB}{
2x1 vector of lower bounds for the x,y region where the convex hull is to be plot. Required.
}
  \item{UB}{
2x1 vector of lower bounds for the x,y region where the convex hull is to be plot. Required.
}
  \item{xlab}{
Label for x axis.
}
  \item{ylab}{
Label for y axis.
}
}

\value{
An integer vector giving the indices of the unique points lying on the convex hull, in clockwise order. (The first will be returned for duplicate points.). Same as returned by chull(xin).
}
\author{
E. del Castillo
}
\note{
Function also plots the convex hull upon return.
}


\examples{
\dontrun{
## Generate some random 2-dimensional point set

X<-cbind(runif(100,-2,2),runif(100,-2,2))

## Compute and plot convex hull

plotConvexHull(xin=X, LB=c(-4,-4), UB=c(4,4), xlab='X', ylab='Y')
}
}
