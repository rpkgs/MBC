#' Random orthogonal rotation
#'
#' Generate a `k`-dimensional random orthogonal rotation matrix.
#'
#' @param k the number of dimensions.
#' 
#' @references 
#' 1. Mezzadri, F. 2007. How to generate random matrices from the
#' classical compact groups, Notices of the American Mathematical Society,
#' 54:592–604.
#' 
#' @export
rot.random <- function(k) {
  rand <- matrix(rnorm(k * k), ncol = k)
  QRd <- qr(rand)
  Q <- qr.Q(QRd)
  R <- qr.R(QRd)
  diagR <- diag(R)
  rot <- Q %*% diag(diagR / abs(diagR))
  return(rot)
}

# Multivariate goodness-of-fit scoring function

#' Energy distance score
#'
#' Calculate the energy distance score measuring the statistical discrepancy
#' between samples `x` and `y` from two multivariate distributions.
#'
#' @param x numeric matrix.
#' @param y numeric matrix.
#' @param scale.x logical indicating whether data should be standardized based
#' on `x`.
#' @param n.cases the number of sub-sampled cases; `NULL` uses all data.
#' @param alpha distance exponent in (0,2]
#' @param method method used to weight the statistics
#' @references 
#' 1. Székely, G.J. and M.L. Rizzo, 2004. Testing for equal
#' distributions in high dimension, InterStat, November (5).
#'
#' 2. Székely, G.J. and M.L. Rizzo, 2013. Energy statistics: statistics based on
#' distances. Journal of Statistical Planning and Inference, 143(8):1249-1272.
#' doi:10.1016/j.jspi.2013.03.018
#' 
#' 3. Rizzo, M.L. and G.L. Székely, 2016. Energy distance. Wiley Interdisciplinary
#' Reviews: Computational Statistics, 8(1):27-38.
#' @export
escore <- function(x, y, scale.x = FALSE, n.cases = NULL, alpha = 1, method = "cluster") {
  n.x <- nrow(x)
  n.y <- nrow(y)
  if (scale.x) {
    x <- scale(x)
    y <- scale(y, center = attr(x, "scaled:center"), scale = attr(
      x,
      "scaled:scale"
    ))
  }
  if (!is.null(n.cases)) {
    n.cases <- min(n.x, n.y, n.cases)
    x <- x[sample(n.x, size = n.cases), , drop = FALSE]
    y <- y[sample(n.y, size = n.cases), , drop = FALSE]
    n.x <- n.cases
    n.y <- n.cases
  }
  edist(rbind(x, y),
    sizes = c(n.x, n.y), distance = FALSE,
    alpha = alpha, method = method
  )[1] / 2
}
