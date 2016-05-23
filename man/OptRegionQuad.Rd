\name{OptRegionQuad}
\alias{OptRegionQuad}

\title{Computes Confidence Regions of Optima of Quadratic Polynomial Models}

\description{
Computes and displays the confidence region on the location of the optima of a
quadratic response surface in 2 factors using bootstrapping.
}
\usage{OptRegionQuad(X, y, nosim=200, alpha=0.05, LB, UB,
triangularRegion=FALSE, vertex1=NULL, vertex2=NULL,
maximization=TRUE,
xlab="Protein eaten, mg",
ylab="Carbohydrates eaten, mg",
outputPDFFile="CRplot.pdf")}

\arguments{
  \item{X}{nx2 matrix with the values of the 2 regressors (experimental factors)
  in the n observations}
  \item{y}{nx1 vector of response value observations}
  \item{nosim}{number of simulations (default=200)}
  \item{alpha}{confidence level (0<alpha<1; default=0.05)}
  \item{LB,UB}{2x1 vectors of lower and upper bounds for search region where
  optima may lie}
  \item{triangularRegion}{logical: if TRUE it will constrain the optimum points
  to lie inside a triangle defined by the coordinates (0,0), and those in
  "vertex1", and "vertex2", see below (in addition to being constrained to lie
  inside the region defined by LB and UB). NOTE: use TRUE when the treatments
  form a triangular experimental region in shape. If FALSE, maxima will only be
  constrained to lie inside the rectangular region defined by LB and UB. Default
  is FALSE.}
  \item{vertex1,vertex2}{2x1 vectors with coordinates defining two of the 3
  vertices of a triangular region. Must be provided if triangularRegion is TRUE
  (NOTE: vertices numbered clockwise, vertex0=c(0,0) always)}
  \item{maximization}{logical: if TRUE (default) it maximizes it FALSE it
  minimizes}
  \item{xlab}{text label for x axis in confidence region plot (default:
  "Protein eaten, mg")}
  \item{ylab}{text label for y axis in confidence region plot (default:
  "Carbohydrates eaten, mg")}
  \item{outputPDFFile}{name of the PDF file where the CR plot is saved
  (default: "CRplot.pdf")}

}
\details{Computes and displays an approximated 100*(1-alpha) percent confidence
region (CR) for the linear-constrained maximum of a quadratic polynomial
regression model in 2 controllable factors. Grey region on output plot is the
approximate CR. The CR is computed as the convex hull of the coordinates of the
optima found from bootstrapping nosim quadratic polynomial regressions to the
data (therefore, it is an approximate CR). The mean value of the optimum is
shown as a red point, and a smoothed contour plot of the X,y data obtained via
thin plate splines is shown as well.

  Usage assuming all default options:

  out<-OptRegionQuad(X=X,y=y,LB=LB,UB=UB)

}
\value{

\item{meanPoint }{a 2x1 vector with the coordinates of the mean optimum point}
\item{xin }{an mx2 matrix with the x,y coordinates of all simulated
  #points that belong to the confidence region (dim(m) is (1-alpha)*nosim)}
%% ...
}
\references{
Del Castillo, E., Hunt, J., Rapkin, J., and Zarmehri, S. , "Confidence
regions for the location of response surface optima: the R package OptimaRegion".
}
\author{Enrique del Castillo \email{exd13@psu.edu}, John Hunt
\email{J.Hunt@exeter.ac.uk} and James Rapkin \email{jr297@exeter.ac.uk}
%%  ~~who you are~~
}

\note{Upon completion, a PDF file containing the CR plot with name as set in
ouputPDFFile is created. To help visualization of the CR and the experimental response data, the plot also shows contours of a smoothed
approximation to the response values, note this is not the quadratic
polynomial that was fit (it is actually a TPS model).
%%  ~~further notes~~
}

\examples{
\dontrun{
## Example 1: randomly generated 2-variable response surface data

X<-cbind(runif(100,-2,2),runif(100,-2,2))

y<-as.matrix(72-11.78*X[,1]+0.74*X[,2]-7.25*X[,1]^2-7.55*X[,2]^2-4.85*X[,1]*X[,2]+
rnorm(100,0,8))

## Find a 95 percent confidence region for the maximum of a quadratic polynomial
## fitted to these data

out<-OptRegionQuad(X=X,y=y,nosim=200,LB=c(-2,-2),UB=c(2,2), xlab="X1",ylab="X2")

## Example 2: a mixture-amount experiment in two components (Drug dataset) with
## non-normal data. Note triangular experimental region. Resulting 95\%
## confidence region is pushed against the constraint and results in a
## "thin line"

out<-OptRegionQuad(X=Drug[,1:2],y=Drug[,3],nosim=500,LB=c(0,0),UB=c(0.08,11),
 xlab="Component 1 (mg.)",ylab="Component 2 (mg.)",triangularRegion = TRUE,
 vertex1 = c(0.02,11),vertex2 = c(0.08,1.8),outputPDFFile="Mixture_plot.pdf")
}}

