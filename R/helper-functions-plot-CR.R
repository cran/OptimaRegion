# Draw a 2D plot for the optima confidence region
#
# @param boot_optima numeric matrix of shape ((1 - alpha)B, 2); it contains the
#                    (1 - alpha)B bootstrap optima projected onto a 2D plane
# @param boost_optimum numeric vector of shape (1, 2); the bootsted optimum
#                      computed by taking the column average of boot_optima
# @param xlab string; it specifies the lable of the horizontal axis
# @param ylab string; it specifies the lable of the vertical axis
# @param xlim numeric vector of shape (1, 2); it specifies the lower and upper
#             limits of the horizontal axis
# @param ylim numeric vector of shape (1, 2); it specifies the lower and upper
#             limits of the vertical axis
# @param main string; it specifies the title of the output figure
# @return a figure displaying the confidence region of the true optimum,
#         along with the boostrap optima and the boosted optimum,
#         projected onto a 2D plane
draw_2D_CR <- function(boot_optima,
                       boost_optimum, # typo, should be bagged optimum
                       xlab, ylab, xlim, ylim) {
  # get the indices of the points that are on the boundary of the convex hull
  id_cvx_hull <- chull(boot_optima)
  id_cvx_hull <- c(id_cvx_hull, id_cvx_hull[1]) # add 1st point to get a loop
  # plot the boundary of the convex hull of the bootstrap optima
  plot(
    boot_optima[id_cvx_hull, ],
    col = "black", cex = 0.01, pch = 16,
    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim
  )
  lines(boot_optima[id_cvx_hull, ], col = "black")
  # plot the boosted optimum
  points(boost_optimum[1], boost_optimum[2], col = "red", cex = 1, pch = 16)
}
# Draw pair-wise projected 2D plots for the optima confidence region
#
# @inheritParams draw_2D_CR
# @inheritParams GloptiPolyRegion
# @return a figure displaying the confidence region of the true optimum,
#         projected onto each pairwise-variable planes
#' @importFrom graphics plot.new
#' @importFrom grDevices dev.new
draw_2D_CRs <- function(boot_optima, boost_optimum, lb, ub,
                        for_dev = TRUE, axes_labels) {
  if (for_dev) dev.new()
  k <- ncol(boot_optima)
  par(mfrow = c(k - 1, k - 1))
  for (i in 1:(k - 1)) { # each row of the sub-figures
    for (j in 2:k) { # sub-figure in each row
      if (i < j) {
        if (is.null(axes_labels)) {
          xlab <- paste("x", j)
          ylab <- paste("x", i)
        } else {
          if(!is.null(axes_labels) && length(axes_labels) != k) {
            stop("Incorrect number of labels!")
          }
          xlab <- axes_labels[j]
          ylab <- axes_labels[i]
        }
        draw_2D_CR(
          boot_optima = boot_optima[, c(j, i)],
          boost_optimum = boost_optimum[c(j, i)],
          xlab = xlab, ylab = ylab,
          xlim = c(lb[j], ub[j]), ylim = c(lb[i], ub[i])
        )
      } else {
        plot.new()
      }
    }
  }
}
