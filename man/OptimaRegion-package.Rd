\name{OptimaRegion-package}
\alias{OptimaRegion-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{OptimaRegion
%%  ~~function to do ... ~~
}
\description{
OptimaRegion is a package for the computation of confidence regions on the location of the optima of response surface models. Both parametric (quadratic polynomial) and nonparametric (thin plate spline) models are supported. The confidence regions obtained do not rely on any distributional assumption, such as Normality of the response, and are obtained by bootstrapping. The resulting regions are both valid and unbiased, and have a size that rapidly decreases as
the sample size increases. Regions are obtained both numerically (as a set of points) and graphically, as the convex hull of the points. Functionality for the
computation of a bootstrap confidence interval on the distance between the
optima of two different response surfaces is included.
}


\references{
Del Castillo, E., Hunt, J., Rapkin, J., and Zarmehri, S. , "Confidence
regions for the location of response surface optima: the R package OptimaRegion".
}
\author{Enrique del Castillo \email{exd13@psu.edu}, John Hunt
\email{J.Hunt@exeter.ac.uk} and James Rapkin \email{jr297@exeter.ac.uk}
%%  ~~who you are~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Response Surface Methodology }
\keyword{ Design of Experiments }% __ONLY ONE__ keyword per line
\keyword{ Thin Plate Splines}
