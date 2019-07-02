# Fit polynomial models
#
# Formulas are entered via the polym function.
# Polym automatically formulates the full polynomial basis given its degree
# and names coefficients with 0/1 indicators, e.g. the name for term
# x1*x2^2 in a 3 variables model is "1.2.0".
# This naming pattern can be used to position the coefficients in the
# multi-dimentional arrays for the GloptiPolyR solver, with the help of
# string-manipulation using regular expression
#
# @inheritParams GloptiPolyRegion
# @return an object of class "lm"
fit_polym <- function(X, y, degree) {
  data <- data.frame(X, y)
  if (ncol(X) == 2) {
    colnames(data) <- c("x1", "x2", "y")
    model <- lm(y ~ polym(x1, x2, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 3) {
    colnames(data) <- c("x1", "x2", "x3", "y")
    model <- lm(y ~ polym(x1, x2, x3, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 4) {
    colnames(data) <- c("x1", "x2", "x3", "x4", "y")
    model <- lm(y ~ polym(x1, x2, x3, x4, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 5) {
    colnames(data) <- c("x1", "x2", "x3", "x4", "x5", "y")
    model <- lm(y ~ polym(x1, x2, x3, x4, x5, degree = degree, raw = TRUE), data = data)
  } else {
    stop("The function only takes 2 - 5 factors.")
  }
  # return
  model
}
# Get positions for monomial coefficients
#
# @param coefficients_name string vector of shape (1, p); it specifies the
#                     coefficient names following the polym pattern, e.g.,
#                     the name for x1*x2^2 in a 3-variable model is "1.2.0"
# @param k integer scalor; it specifies the number of variables
# @return integer matrix of shape (p, k); its (i, j) element speficies the
#         (power + 1) value of the jth variable in the ith monomial term,
#         (power + 1) accommodating the zero power; its ith row specifies
#         the position of the coefficient of the ith nomomial term in the
#         multi-dimensional array of the GloptiPolyR solver
#' @importFrom magrittr "%>%"
coef_name_to_array_index <- function(coefficients_name, k) {
  array_index_string <- stringr::str_extract(coefficients_name, "(\\d\\.)+[\\d]")
  array_index_number <- matrix(NA, length(array_index_string), k)
  array_index_number[1, ] <- 1
  for (i in 2:length(array_index_string)) {
    array_index_number[i, ] <- array_index_string[i] %>%
      stringr::str_split("\\.") %>%
      unlist() %>%
      as.numeric() + 1
  }
  # return
  array_index_number
}
# Optimize fitted polynomial functions via GloptiPolyR
#
# @param coefficients numeric vector of shape (1, p); it specifies the the
#                     coefficients of an "lm" objected formulated with the
#                     polym function
# @param k integer scalor; it specifies the number of variables
# @inheritParams GloptiPolyRegion
# @return the optimal solution and its corresponding objective value
Gloptipolym <- function(coefficients, k, degree, lb, ub, maximization) {
  Ps <- list() # argument for GloptiPolyR, a list of lists
  # Objective function ------------------------------------------------------
  P <- list()
  c <- array(0, dim = rep(degree + 1, k))
  # get position indices for the coefficients of the objective function
  id <- coef_name_to_array_index(names(coefficients), k = k)
  # put coefficient values into the multi-dimensional array
  for (i in 1:nrow(id)) {
    eval(parse(text = paste(
      "c[", toString(id[i, ]),
      "] <- coefficients[", i, "]"
    )))
    # example 1: eval(parse(text = "1+1")) -> 2
    # example 2: toString(id[1,]) -> "1, 1, 1"
  }
  if (maximization) { # assume GloptiPolyR only supports "min"
    P$c <- -c
  } else {
    P$c <- c
  }
  P$t <- "min" # specify attribute
  Ps[[1]] <- P # add to list
  # Constraint functions ----------------------------------------------------
  for (i in 1:k) { # loop through variables
    # Lower bound constraint
    P <- list()
    c <- array(0, dim = rep(degree + 1, k))
    # specify coefficient 1 of the variable
    index_for_c <- rep(1, k)
    index_for_c[i] <- 2
    eval(parse(text = paste("c[", toString(index_for_c), "] <- 1")))
    # specify the constraint constant
    eval(parse(text = paste("c[", toString(rep(1, k)), "] <- -lb[", i, "]")))
    P$c <- c
    P$t <- ">=" # specify attribute
    Ps[[2 * i]] <- P # add to list
    # Upper bound constraint
    P <- list()
    c <- array(0, dim = rep(degree + 1, k))
    # specify coefficient 1 of the variable
    index_for_c <- rep(1, k)
    index_for_c[i] <- 2
    eval(parse(text = paste("c[", toString(index_for_c), "] <- 1")))
    # specify the constraint constant
    eval(parse(text = paste("c[", toString(rep(1, k)), "] <- -ub[", i, "]")))
    P$c <- c
    P$t <- "<="
    Ps[[2 * i + 1]] <- P
  }
  # Call GloptiPolyR --------------------------------------------------------
  res <- GloptiPolyR(Ps)
  solution <- res$solution
  if (maximization) { # assume GloptiPolyR only supports "min"
    objective <- -res$objective
  } else {
    objective <- res$objective
  }
  # return
  list(solution = solution, objective = objective)
}
