#' Confidence region for optima of Thin Plate Spline Models (2 regressors)
#'
#' Computes and displays an approximated (1 - alpha)*100% confidence region (CR) for
#' the linear-constrained maximum of a penalized Thin Plate Spline (TPS) model
#' in 2 controllable factors
#' \insertCite{DelCastilloCR}{OptimaRegion}.
#' Generates a PDF file with a graph displaying the CR.
#' Grey region on output plot is the approximate CR. The mean coordinates (centroid)
#' of the optima is shown as a red point.
#'
#' This program approximates the confidence region (CR) of the location of the optimum
#' of a Thin Plate Spline (TPS) in 2 regressors x constrained inside  a rectangular
#' region defined by LB and UB. If triangularRegion=TRUE it will also contrain the
#' optimum to lie inside the experimental region assumed to be well approximated by
#' a triangle. The CR is generated pointwise by bootstrapping the residuals of a TPS fit
#' to the given (X,y) data, refitting Tps models, and solving the corresponding
#' constrained maximization (or minimization) problems. The confidence region is
#' approximated by the convex hull of all the optimal solutions found. The CR
#' computation is based on the "CS" bootstrapping approach for building a confidence set
#' of a parametric function described in
#' \insertCite{WoutersenHam2013;textual}{OptimaRegion}.
#' This version of the program
#' uses nonparametric bootstrapping confidence regions to get the Confidence region of
#' the Tps parameters,using the notion of data depth according to
#' \insertCite{yeh1997balanced;textual}{OptimaRegion}.
#' Hence, this version does not rely on the normality assumption of the data.
#' The TPS models are fit using the "fields" R package
#' \insertCite{fieldsPackage}{OptimaRegion}
#' and its "Tps" function.
#'
#' @inheritParams OptRegionQuad
#' @param lambda penalization parameter (larger values implies more smoothing).
#'               Default is 0.04
#' @param outputOptimaFile name of the text file containing the coordinates of all
#'        the optima found (same information as in output vector xin, see below)
#' @return Upon completion, a PDF file containing the CR plot with name as set in
#'         ouputPDFFile is created and a text file with all optima in the CR is created too.
#'         Also, the function returns a list containing
#'         the following 2 components:
#'         \describe{
#'           \item{meanPoint}{a 2x1 vector with the coordinates of the mean optimum
#'                              point (displayed as a red dot in the CR plot in output
#'                              PDF file)}
#'           \item{xin}{an mx2 matrix with the x,y coordinates of all simulated
#'                        points that belong to the confidence region (dim(m) is
#'                        (1-alpha)*nosim)}
#'         }
#' @inheritSection OptRegionQuad Author(s)
#' @references{
#'  \insertAllCited{}
#' }
#' @examples
#' \dontrun{
#' # Example 1: randomly generated 2-variable response surface data
#' X <- cbind(runif(100, -2, 2), runif(100, -2, 2))
#' y <- as.matrix(72 - 11.78 * X[, 1] + 0.74 * X[, 2] - 7.25 * X[, 1]^2 -
#'   7.55 * X[, 2]^2 - 4.85 * X[, 1] * X[, 2] + rnorm(100, 0, 8))
#' # Find a 95 percent confidence region for the maximum of a Thin Plate Spline
#' # model fitted to these data
#' out <- OptRegionTps(
#'   X = X, y = y, nosim = 200, LB = c(-2, -2), UB = c(2, 2),
#'   xlab = "X1", ylab = "X2"
#' )
#'
#' # Example 2: a mixture-amount experiment in two components (Drug dataset) with
#' # non-normal data. Note triangular experimental region. Resulting 95p confidence
#' # region of the maxima of a TPS model has area > 0. Contrast with region for
#' # quadratic polynomial model. Note: 500 bootstrap iterations may take a few minutes.
#' out <- OptRegionTps(
#'   X = Drug[, 1:2], y = Drug[, 3], nosim = 500, lambda = 0.05, LB = c(0, 0),
#'   UB = c(0.08, 11), xlab = "Component 1 (mg.)", ylab = "Component 2 (mg.)",
#'   triangularRegion = TRUE, vertex1 = c(0.02, 11), vertex2 = c(0.08, 1.8),
#'   outputPDFFile = "Mixture_plot.pdf"
#' )
#' }
#' @importFrom Rdpack reprompt
#' @importFrom grDevices chull dev.off heat.colors pdf
#' @importFrom graphics contour image lines par plot points polygon
#' @importFrom stats fitted lm resid vcov
#' @export
OptRegionTps <- function(X, y, lambda = 0.04, nosim = 1000, alpha = 0.05, LB, UB,
                         triangularRegion = FALSE, vertex1 = NULL, vertex2 = NULL,
                         maximization = TRUE,
                         xlab = "Protein eaten, mg", ylab = "Carbohydrate eaten, mg",
                         outputPDFFile = "CRplot.pdf", outputOptimaFile = "Optima.txt") {
  # Check that X matrix has k=2 factors
  k <- dim(X)[2]
  if ((k > 2) | (k < 2)) stop("Error. Number of factors must equal to 2")

  # If experimental region was specified as triangular, compute the parameteres defining the 3 lines that approximate its shape.
  epsi <-0.001
  if (triangularRegion) {
    x11p <- vertex1[1] # user defined vertices; vertex 1 and 2 are clockwise on the plane; third vertex is (0,0)
    x21 <- vertex1[2]
    x12 <- vertex2[1]
    x22 <- vertex2[2]
    m1 <- x21 / (x11p+epsi) # to avoid error if vertex1[1] =0
    m2 <- x22 / (x12+epsi)
    m3 <- (x21 - x22) / (x11p - x12)
    bintercept <- x22 - m3 * x12
  } else {
    m1 <- 0
    m2 <- 0
    m3 <- 0
    bintercept <- 0
  }

  # Load some libraries
  # t<-.libPaths()
  # library("nloptr", lib.loc=t)
  # library("lhs", lib.loc=t)
  # library("fields",lib.loc=t)
  # library("DepthProc",lib.loc=t) #better package than "depth"
  # library("scatterplot3d",lib.loc=t)
  # library("MASS",lib.loc=t)

  # First eliminate replicates and compute y-averages
  model <- fields::Tps(X, y, lambda = lambda)
  out <- fields::Krig.replicates(model)
  X <- out$xM
  y <- out$yM

  # Find parameters and covariance matrix of the fitted Thin Plate Spline (TPS) model
  model <- fields::Tps(X, y, lambda = lambda)

  fit <- model$fitted.values
  nn <- length(y)
  p <- 3 # true for 2-dimensional TPS's
  A <- fields::Krig.Amatrix(model)
  df <- 2 * sum(diag(A)) - sum(diag(A %*% t(A))) # degrees of freedom according to Ruppert, Wand and Carroll 2003
  # Compute studentized residuals and M matrix, where w2=M*Y is an (n-p) vector; these are the only really free to vary parameters in the penalized spline model
  T <- fields::Krig.null.function(model$knots, m = 2) # TPS basis functions vector[1,x1,x2] for a TPS in 2 dimensions
  Sigma <- fields::Rad.cov(model$knots, model$knots, p = 2, m = 2) # radial basis covariance; optimal smoothed vector is Yhat= T*d + Sigma*c
  Q <- qr.Q(model$matrices$qr.T, complete = TRUE) # extract QR decomposition of T. This is done by the Tps(Krig) function
  F2 <- Q[, (p + 1):nn] # needed since c=F2*w2
  M <- solve(t(F2) %*% (Sigma + lambda * diag(nn)) %*% F2) %*% t(F2)
  w2Hat <- M %*% y
  MMt <- M %*% t(M)
  den <- sqrt(diag((diag(nn) - A) %*% (diag(nn) - A))) # effective degrees of freedom for the residuals; close but not equal to I -tr(A)
  res <- model$residuals / den # residuals corrected for small bias, following Kauermann et al. 2009
  # Nonparametric bootstrapping of the residuals used to simulate w2 vectors (and hence, c and d)
  w2vectors <- matrix(nrow = nosim, ncol = nn - p)
  dstar <- matrix(nrow = nosim, ncol = p)
  w2vectorsStd <- matrix(nrow = nosim, ncol = nn - p)
  print("Bootstrapping...")
  for (i in 1:nosim) {
    indices <- sample(1:nn, replace = TRUE)
    ystar <- fit + res[indices]
    # Now fit Tps model to simulated data
    modelstar <- fields::Tps(X, ystar, lambda = lambda)
    w2vectors[i, ] <- M %*% ystar
    dstar[i, ] <- t(modelstar$d) # save also the d vectors corresponding to each simulated w2; it is easier this way
    # The smoothing A matrix in Yhat=A*Y is the same as using the original data, as it only depends on X and lambda and we are only changing (simulating) y
    resstar <- modelstar$residuals / den # following Kauermann et al. 2009
    sigma2star <- norm(resstar) / df
    Vi <- sigma2star * MMt # =Cov(w2 hat)
    e <- eigen(Vi, symmetric = TRUE)
    # find square root matrix and standardize w2 vectors following Yeh and Singh, J.Royal Stat. Soc. 59(3), pp. 639-652, 1997
    Ssqrtinv <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors))
    w2vectorsStd[i, ] <- Ssqrtinv %*% (w2vectors[i, ] - w2Hat) * sqrt(nn - p)
    print(i)
  }
  deep <- vector(length = nosim)
  w2vectorsStdMat <- as.matrix(w2vectorsStd)
  # Calculate data depth of each parameter vector
  d <- DepthProc::depthTukey(w2vectorsStdMat, w2vectorsStdMat, ndir = 3000)
  order.d <- order(d) # order only based on Tukey's depth
  # find number of points should be in the alpha percent confidence region
  ind.alpha <- alpha * nosim + 1
  # provide the indices of points in the confidence region
  indices <- order.d[(ind.alpha:nosim)]
  w2vectors_In <- w2vectors[indices, ]
  dstar_In <- dstar[indices, ]
  # Plot simulated betas inside and outside the CR--uncomment to check some simulated model parameter in a matrix scatterplot
  # library("car",lib.loc=t)
  # library("lattice",lib.loc=t)
  # x11()
  # scatterplotMatrix(w2vectors[,(1:10)],groups=deep>TenPerc,col=c("red","blue"))
  # x11()
  # optimize model subject to bounds
  l <- dim(w2vectors_In)
  ################## OPTIMIZATION
  # First create initial points for optimizer
  delta1 <- 0.1 * (UB[1] - LB[1])
  delta2 <- 0.1 * (UB[2] - LB[2])
  X0 <- matrix(nrow = 5, ncol = k)
  if (triangularRegion) {
    X0[1, ] <- c(LB[1] + delta1, LB[2] + delta2)
    X0[2, ] <- vertex1 - c(0, delta2)
    X0[3, ] <- vertex2 - c(delta1, 0)
    X0[4, ] <- c(mean(X0[1:3, 1]), mean(X0[1:3, 2]))
    # make sure points are inside LB and UB
    X0[, 2] <- apply(cbind(X0[, 2], UB[2]), 1, min)
    X0[, 2] <- apply(cbind(X0[, 2], LB[2]), 1, max)
    X0[, 1] <- apply(cbind(X0[, 1], UB[1]), 1, min)
    X0[, 1] <- apply(cbind(X0[, 1], LB[1]), 1, max)
    noinitial <- 4
    numberFeasible <- 4
  } else {
    X0[1, ] <- c(LB[1] + delta1, LB[2] + delta2)
    X0[2, ] <- c(LB[1] + delta1, UB[2] - delta2)
    X0[3, ] <- c(UB[1] - delta1, LB[2] + delta2)
    X0[4, ] <- c(UB[1] - delta1, UB[2] - delta2)
    X0[5, ] <- c(mean(X0[1:4, 1]), mean(X0[1:4, 2]))
    noinitial <- 5
    numberFeasible <- 5
  }
  xin <- matrix(nrow = l[1], ncol = k)
  # Main loop
  print("Computing optima...")
  for (m in 1:l[1]) {
    # Pick a pair (c,d) that correspond to a w2 inside the alpha% data depth region
    cCoef <- F2 %*% w2vectors_In[m, ]
    dCoef <- dstar_In[m, ]
    best <- 1e50
    for (j in 1:numberFeasible) {
      # Setup optimization problem depending on whether the experimental region has been defined as triangular or not
      # Note that we are using an optimizer that does not require the gradient of the objective function --this could be modified if expression of grad(fhat) are available for a TPS predictor
      if (triangularRegion) {
        local_opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-3)
        opts <- list("algorithm" = "NLOPT_LN_AUGLAG", print_level = 0, "local_opts" = local_opts)
        out <- nloptr::nloptr(X0[j, ], eval_f = computefTps, eval_grad_f = NULL, eval_jac_g_ineq = NULL, lb = LB, ub = UB, eval_g_eq = NULL, eval_g_ineq = constraintsTps, opts = opts, cCoef = cCoef, dCoef = dCoef, maximization = maximization, m1 = m1, m2 = m2, m3 = m3, bintercept = bintercept, knots = model$knots, center = model$transform$x.center, scale = model$transform$x.scale)
      } else {
        out <- nloptr::nloptr(X0[j, ], eval_f = computefTps, eval_grad_f = NULL, lb = LB, ub = UB, eval_g_eq = NULL, eval_g_ineq = NULL, opts = list("algorithm" = "NLOPT_LN_COBYLA", print_level = 0, xtol_rel = 1e-03), cCoef = cCoef, dCoef = dCoef, maximization = maximization, m1 = m1, m2 = m2, m3 = m3, bintercept = bintercept, knots = model$knots, center = model$transform$x.center, scale = model$transform$x.scale)
      }
      if ((out$objective < best) & (out$status > 0)) {
        best <- out$objective
        bestSol <- out$sol
        bestStatus <- out$status
      }
      # print(c('in=',X0[j,],'out=',bestSol))
    } # endfor j
    # save best solution found among all tries for simulated parameter set m
    xin[m, ] <- bestSol
    print(c(m, best, xin[m, ], bestStatus))
  } # endfor m

  # Plot CR and thin plate spline fit to the experimental data on output file
 pdf(file = outputPDFFile, 5.5, 5.5) # comment out to have output to screen
  #x11() # uncomment to have output to screen
  # library("splancs",lib.loc=t)
  # library("maptools",lib.loc=t)
  # library("Hmisc",lib.loc=t)
  # Draw Convex Hull of optima (approximates the CR)
  coords <- plotConvexHull(xin, LB, UB, xlab, ylab)
  par(new = TRUE)
  par(cex.axis = 1.35, cex.lab = 1.5)
  par(xaxt = "n", yaxt = "n")
  # draw centroid of CR
  centroid <- apply(xin, 2, mean)
  points(centroid[1], centroid[2], col = "red", pch = 19)
  par(new = TRUE)
  par(cex.axis = 1.35, cex.lab = 1.5)
  par(xaxt = "n", yaxt = "n")
  # Draw contour plot of Tps fitted to available data
  tpsfit <- fields::Tps(X, y, lambda = lambda)
  surface <- fields::predictSurface(tpsfit)
  image(surface, lwd = 2, col = heat.colors(0,alpha=1), cex.axis = 1.35, cex.lab = 1.5, xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]))
  contour(surface,  drawlabels = TRUE, lwd = 2, cex.axis = 1.35, cex.lab = 1.5, xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]),add = TRUE)
  # par(new=TRUE)
  # par(cex.axis=1.35, cex.lab=1.5)
  # par(xaxt='n', yaxt='n')
  # Draw arrows
  # arrows(29.51,59.12+1.96,29.51,59.12-1.96,code=3,angle=90,length=0.015,col="black",lwd=2)
  # par(new=TRUE)
  # arrows(29.51+1.30,59.12,29.51-1.30,59.12,code=3,angle=90,length=0.015,col="black",lwd=2)
  dev.off()
  # return optimal points in text file
  write(t(xin), file = outputOptimaFile, ncolumns = 2)
  # detach()
  return(list(meanPoint = centroid, xin = xin))
} # end main program


constraintsTps <- function(x, cCoef, dCoef, maximization, m1, m2, m3, bintercept, knots, center, scale) {
  # Computes the constraints limiting the triangular-like experimental region (approximated with a triangle, so 3 constraints)
  z <- vector(length = 3)
  z[1] <- x[2] - m1 * x[1]
  z[2] <- m2 * x[1] - x[2]
  z[3] <- x[2] - bintercept - m3 * x[1]
  return(z)
}

plotConvexHull <- function(xin, LB, UB, xlab = "X", ylab = "Y") {
  # Plots the convex hull of the points in vector xin
  plot(0, 0, col = "white", xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]), xlab = xlab, ylab = ylab)
  hpts_original <- chull(xin)
  hpts_closed <- c(hpts_original, hpts_original[1])
  lines(xin[hpts_closed, ], col = "blue")
  polygon(xin[hpts_closed, 1], xin[hpts_closed, 2], col = "grey")
  return(hpts_original)
}


computefTps <- function(x, cCoef, dCoef, maximization, m1, m2, m3, bintercept, knots, center, scale) {
  # Returns -f(x)--a Tps (predicted value at x), a scalar
  k <- length(x)
  xscaled <- scale(t(x), center, scale) # need to do this
  to <- fields::fields.mkpoly(xscaled, m = 2) # a TPS in 2 dimensions
  Psi <- fields::Rad.cov(xscaled, knots, p = 2, m = 2)
  f <- to %*% dCoef + Psi %*% cCoef # prediction Yhat(x0)
  if (maximization) {
    return(-f)
  } else {
    return(f)
  }
}
