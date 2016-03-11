\name{OptRegionTps}
\alias{OptRegionTps}

\title{Computes Confidence Regions of Optima of Thin Plate Spline Models}

\description{
Computes and displays the confidence region on the location of the optima of a thin plate spline surface in 2 factors using bootstrapping.
}
\usage{OptRegionTps(X, y, lambda=0.04, nosim=1000, alpha=0.05, LB, UB,
triangularRegion=FALSE, vertex1=NULL, vertex2=NULL, maximization=TRUE,
xlab="Protein eaten, mg",
ylab="Carbohydrate eaten, mg",
outputPDFFile="CRplot.pdf",outputOptimaFile="Optima.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{nx2 matrix with the values of the 2 regressors (experimental factors) in the n observations}
  \item{y}{nx1 vector of response value observations}
  \item{lambda}{penalization parameter (larger values implies more smoothing). Default is 0.04}
  \item{nosim}{number of simulations (default=1000)}
  \item{alpha}{confidence level (0<alpha<1; default=0.05)}
  \item{LB,UB}{2x1 vectors of lower and upper bounds for search region where optima may lie}
  \item{triangularRegion}{logical: if TRUE it will constrain the optimum points to lie inside a triangle defined by the coordinates (0,0), and those in "vertex1", and "vertex2", see below (in addition to being constrained to lie inside the region defined by LB and UB). NOTE: use TRUE when the treatments form a triangular experimental region in shape. If FALSE, maxima will only be constrained to lie inside the rectangular region defined by LB and UB. Default is FALSE.}
  \item{vertex1,vertex2}{2x1 vectors with coordinates defining two of the 3 vertices of a triangular region. Must be provided if triangularRegion is TRUE (NOTE: vertices numbered clockwise, vertex0=c(0,0) always)}
  \item{maximization}{logical: if TRUE (default) it maximizes it FALSE it minimizes}
  \item{xlab}{text label for x axis in confidence region plot (default: "Protein eaten, mg")}
  \item{ylab}{text label for y axis in confidence region plot (default: "Carbohydrates eaten, mg")}
  \item{outputPDFFile}{name of the PDF file where the CR plot is saved (default: "CRplot.pdf")}
  \item{outputOptimaFile}{name of the text file containing the coordinates of all the optima found (same information as in output vector xin, see below)}

}
\details{This program approximates the confidence region (CR) of the location of the optimum of a Thin Plate Spline (TPS) in 2 regressors x constrained inside  a rectangular region defined by LB and UB. If triangularRegion=TRUE it will also contrain the optimum to lie inside the experimental region assumed to be well approximated by a triangle. The CR is generated pointwise by bootstrapping the residuals of a TPS fit to the given (X,y) data, refitting Tps models, and solving the corresponding constrained maximization (or minimization) problems. The confidence region is approximated by the convex hull of all the solutions (x1*,x2*) found.

  Usage assuming all default options:

  out<-OptRegionTps(X=X,y=y,LB=LB,UB=UB)

}
\value{

\item{meanPoint }{a 2x1 vector with the coordinates of the mean optimum point}
\item{xin }{an mx2 matrix with the x,y coordinates of all simulated
  #points that belong to the confidence region (dim(m) is (1-alpha)*nosim)}
%% ...
}
\references{
Del Castillo, E., Hunt, J., and Rapkin, J. , "Data-depth bootstrapped confidence
regions for the optima of parametric and nonparametric response surfaces".
}
\author{Enrique del Castillo, \email{exd13@psu.edu}
%%  ~~who you are~~
}
\note{Upon completion, a PDF file containing the CR plot with name as set in ouputPDFFile  is created. A text file named as set in outputOptimaFile with all xin values is created too.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
## Example 1: randomly generated 2-variable response surface data

X<-cbind(runif(100,-2,2),runif(100,-2,2))

y<-as.matrix(72-11.78*X[,1]+0.74*X[,2]-7.25*X[,1]^2-7.55*X[,2]^2-4.85*X[,1]*X[,2]+
rnorm(100,0,8))

## Find a 95 percent confidence region for the maximum of a Thin Plate Spline
model fitted to these data

out<-OptRegionTps(X=X,y=y,nosim=200,LB=c(-2,-2),UB=c(2,2), xlab="X1",ylab="X2",
outputPDFFile="CR_plot.pdf")

## Example 2: a mixture-amount experiment in two components (Drug dataset) with
non-normal data. Note triangular experimental region. Resulting 95p confidence
region of the maxima of a TPS model has area > 0. Contrast with region for
quadratic polynomial model. Note: 500 bootstrap iterations may take a few minutes.

out<-OptRegionTps(X=Drug[,1:2],y=Drug[,3],nosim=500,lambda=0.05,LB=c(0,0),UB=
 c(0.08,11), xlab="Component 1 (mg.)",ylab="Component 2 (mg.)",triangularRegion
 = TRUE,vertex1 = c(0.02,11),vertex2 = c(0.08,1.8),outputPDFFile=
 "Mixture_plot.pdf")
}
}
