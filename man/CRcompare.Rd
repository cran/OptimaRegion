\name{CRcompare}
\alias{CRcompare}

\title{Computes bootstrapped confidence intervals for the distance between two
response surface optima}

\description{
Computes bootstrapped confidence intervals for the mean and median
distance between the optima of two response surface models. Models can be thin
plate splines or quadratic polynomials.
}
\usage{CRcompare(X1,y1,X2,y2,lambda=0.04,responseType='TPS',nosim1and2=200,
alpha=0.05,LB1,UB1,triangularRegion1=FALSE, vertex11=NULL, vertex21=NULL,
maximization1=TRUE,outputPDFFile1="CR_plot1.pdf",outputOptimafile1="Optima1.txt"
,LB2,UB2,triangularRegion2=FALSE, vertex12=NULL, vertex22=NULL, maximization2=
TRUE,outputPDFFile2="CR_plot2.pdf",outputOptimafile2="Optima2.txt",xlab1and2=
"Protein eaten (mg)",ylab1and2="Carbohydrates eaten (mg)")
}

\arguments{
\item{X1}{nx2 matrix with the values of the 2 regressors (experimental factors)
  in the n observations associated with response 1}
  \item{y1}{nx1 vector of response no. 1 value observations}
  \item{X2}{nx2 matrix with the values of the 2 regressors (experimental factors)
  in the n observations associated with response 2}
  \item{y2}{nx1 vector of response no. 2 value observations}
  \item{lambda}{smoothing penalty if a TPS model is selected (default=0.04)}
  \item{responseType}{use 'TPS' if fitting thin plate spline responses, 'Quad'
  if fitting quadratic polynomials}
  \item{nosim1and2}{number of simulations (default=200)}
  \item{alpha}{confidence level (0<alpha<1; default=0.05)}
  \item{LB1,UB1}{2x1 vectors of lower and upper bounds for search region where
  optima may lie for response 1}
  \item{triangularRegion1}{logical: if TRUE it will constrain the optimum points
  of response 1 to lie inside a triangle defined by the coordinates (0,0), and
  those in "vertex11", and "vertex21", see below). NOTE: use TRUE when the
  treatments form a triangular experimental region in shape. If FALSE, maxima
  will only be constrained to lie inside the rectangular region defined by LB1
  and UB1. Default is FALSE.}
  \item{vertex11,vertex21}{2x1 vectors with coordinates defining two of the 3
  vertices of a triangular region for searching the optima of response 1.
  Must be provided if triangularRegion is TRUE
  (NOTE: vertices numbered clockwise, vertex0=c(0,0) always)}
  \item{maximization1}{logical: if TRUE (default) it maximizes response 1, if
  FALSE it minimizes it}
  \item{outputPDFFile1}{name of the PDF file where the CR plot of response 1 is
  saved (default: "CRplot1.pdf")}
  \item{outputOptimafile1}{name of text file for saving the coordinates of
  the optima of response 1 (default: "Optima1.txt")}
  \item{LB2,UB2}{2x1 vectors of lower and upper bounds for search region where
  optima may lie for response 2}
  \item{triangularRegion2}{logical: if TRUE it will constrain the optimum points
  of response 2 to lie inside a triangle defined by the coordinates (0,0), and
  those in "vertex12", and "vertex22", see below (in addition to being constrained
  to lie inside the region defined by LB and UB). NOTE: use TRUE when the
  treatments form a triangular experimental region in shape. If FALSE, maxima
  will only be constrained to lie inside the rectangular region defined by LB2
  and UB2. Default is FALSE.}
  \item{vertex12}{see vertex22}
  \item{vertex22}{2x1 vectors with coordinates defining two of the 3
  vertices of a triangular region for searching the optima of response 2.
  Must be provided if triangularRegion is TRUE
  (NOTE: vertices numbered clockwise, vertex0=c(0,0) always)}
  \item{maximization2}{logical: if TRUE (default) it maximizes response 2 if
  FALSE it minimizes it}
  \item{outputPDFFile2}{name of the PDF file where the CR plot of response 2 is
  saved (default: "CRplot2.pdf")}
   \item{outputOptimafile2}{name of text file for saving the coordinates of
  the optima of response 2 (default: "Optima2.txt")}
  \item{xlab1and2}{text label for x axis in confidence region plot for the
  optima of each response (default: "Protein eaten, mg")}
  \item{ylab1and2}{text label for y axis in confidence region plot for the
  optima of each response (default: "Carbohydrates eaten, mg")}
}
\details{
Computes distribution-free bootstrapped confidence intervals on the mean and
median distance between the optima of two different responses. The responses can
be Thin Plate Spline models or Quadratic polynomial models. Program calls
OptRegionTps.R or OptRegionQuad.R to compute confidence regions on the optima of
each response, next computes all pairwise distances between points in each CR,
and finally bootstraps the distances to compute bca bootstrapped confidence
intervals for the mean and median distance.

 Usage assuming all default options:

 out<-CRcompare(X1=X1,y1=y1,X2=X2,y2=y2,LB1=LB1,UB1=UB1,LB2=LB2,UB2=UB2)
}
\value{
\item{dist}{vector of distances between pairs of points taken from each set of
optima}
\item{mean}{mean of all pairwise distances}
\item{median}{median of all pairwise distances}
\item{ciMean}{(1-alpha)*100\% confidence interval for the mean of the distances
using bca bootstrapping}
\item{ciMedian}{(1-alpha)*100\% confidence interval for the mean of the distances
using bca bootstrapping.

Note: ciMean and ciMedian are vectors with 5 columns, containing the
significance level, the next two containing the indices of the order statistics
used in the calculations and the final two the calculated endpoints of the CI's.
}
}
\references{
Del Castillo, E., Hunt, J., Rapkin, J., and Zarmehri, S. , "Confidence
regions for the location of response surface optima: the R package OptimaRegion".
}
\author{Enrique del Castillo \email{exd13@psu.edu}, John Hunt
\email{J.Hunt@exeter.ac.uk} and James Rapkin \email{jr297@exeter.ac.uk}
%%  ~~who you are~~
}
\note{Upon completion, PDF files containing the CR plots of the optima of each
response  with name as set in ouputPDFFile1 and outPDFFile2 are created. Two
text files named as set in outputOptimaFile1 and outputOptimaFile2 are also
created. These contain all the individual optima of each response in each CR.}


\examples{
\dontrun{
## Example: two randomly generated data sets, quadratic polynomial responses.

X1<-cbind(runif(100,-2,2),runif(100,-2,2))

y1<-as.matrix(72-11.78*X1[,1]+0.74*X1[,2]-7.25*X1[,1]^2-7.55*X1[,2]^2-4.85*X1[,1]*X1[,2]+
                  rnorm(100,0,8))

X2<-cbind(runif(100,-2,2),runif(100,-2,2))

y2<-as.matrix(72-11.78*X2[,1]+0.74*X2[,2]-7.25*X2[,1]^2-7.55*X2[,2]^2-4.85*X2[,1]*X2[,2]+
                  rnorm(100,0,8))

out<-CRcompare(X1=X1,y1=y1,X2=X2,y2=y2,responseType='Quad',nosim1and2=200,
alpha=0.05,LB1=c(-2,-2),UB1=c(2,2),LB2=c(-2,-2),UB2=c(2,2))
}
}


