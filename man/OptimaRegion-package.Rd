% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptimaRegion.R
\name{OptimaRegion-package}
\alias{OptimaRegion-package}
\title{OptimaRegion package description}
\description{
OptimaRegion is a package for the computation of confidence regions on the location
of the optima of response surface models
\insertCite{DelCastilloCR}{OptimaRegion}.
Both parametric (quadratic polynomial) and
nonparametric (thin plate spline) models are supported. The confidence regions
obtained do not rely on any distributional assumption, such as Normality of the
response, and are obtained by bootstrapping. The resulting regions are both valid and
unbiased, and have a size that rapidly decreases as the sample size increases.
Regions are obtained both numerically (as a set of points) and graphically, as the
convex hull of the points. Functionality for the computation of a bootstrap
confidence interval on the distance between the optima of two different response
surfaces is included.
}
\details{
The OptimaRegion package provides five main functions:
OptRegionQuad, OptRegionTps, GloptiPolyRegion, CRcompare and GloptipolyR.

OptRegionQuad computes distribution-free bootstrapped confidence regions for the
location of the optima of a quadratic polynomial model in 2 regressors.
OptRegionTps computes distribution-free bootstrapped confidence regions for the
location of the optima of a Thin Plate Spline model in 2 regressors.
GloptiPolyRegion computes distribution-free bootstrapped CRs for the location of
global optima for polynomial models up to cubic order in up to 5 regressors.
CRcompare computes bootstrapped confidence intervals for the distance between the
optima of two different response surface models, either quadratic polynomials or
thin plate spline models.
GloptipolyR is an R implementation of the “Gloptipoly” algorithm
\insertCite{lasserre2001global}{OptimaRegion}
for global optimization of polynomial equations subject to bounds.
}
\references{
{
  \insertAllCited{}
}
}
