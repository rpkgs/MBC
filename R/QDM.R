
#' Univariate bias correction via quantile delta mapping
#'
#' Univariate bias correction based on the quantile delta mapping `QDM`
#' version of the quantile mapping algorithm from Cannon et al. (2015).
#' `QDM` constrains model-projected changes in quantiles to be preserved
#' following bias correction by quantile mapping.
#'
#' @param o.c vector of observed samples during the calibration period.
#' @param m.c vector of model outputs during the calibration period.
#' @param m.p vector of model outputs during the projected period.
#' @param ratio logical value indicating if samples are of a ratio quantity
#' (e.g., precipitation).
#' @param trace numeric value indicating the threshold below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric value of a threshold used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric value indicating the maximum proportional change
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric value of a trace threshold used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ECBC logical value indicating whether `mhat.p` outputs should be
#' ordered according to `o.c` ranks, i.e., as in the empirical copula-bias
#' correction (ECBC) algorithm.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' @return a list of with elements consisting of: \item{mhat.c}{vector of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{vector of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[MBCp], [MBCr], [MRS], [escore]`
#' 
#' @references
#' 1. Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias
#' correction of simulated precipitation by quantile mapping: How well do
#' methods preserve relative changes in quantiles and extremes? Journal of
#' Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
#' 
#' @export QDM
QDM <- function(o.c, m.c, m.p, ratio = FALSE, trace = 0.05, trace.calc = 0.5 * trace,
                jitter.factor = 0, n.tau = NULL, ratio.max = 2, ratio.max.trace = 10 * trace,
                ECBC = FALSE, ties = "first", subsample = NULL, pp.type = 7) {
  # If jitter.factor > 0, apply a small amount of jitter to accommodate ties
  # due to limited measurement precision
  if (jitter.factor == 0 &&
    (length(unique(o.c)) == 1 ||
      length(unique(m.c)) == 1 ||
      length(unique(m.p)) == 1)) {
    jitter.factor <- sqrt(.Machine$double.eps)
  }
  if (jitter.factor > 0) {
    o.c <- jitter(o.c, jitter.factor)
    m.c <- jitter(m.c, jitter.factor)
    m.p <- jitter(m.p, jitter.factor)
  }
  # For ratio data, treat exact zeros as left censored values less than
  # trace.calc
  if (ratio) {
    epsilon <- .Machine$double.eps
    o.c[o.c < trace.calc] <- runif(sum(o.c < trace.calc),
      min = epsilon,
      max = trace.calc
    )
    m.c[m.c < trace.calc] <- runif(sum(m.c < trace.calc),
      min = epsilon,
      max = trace.calc
    )
    m.p[m.p < trace.calc] <- runif(sum(m.p < trace.calc),
      min = epsilon,
      max = trace.calc
    )
  }
  # Calculate empirical quantiles
  n <- length(m.p)
  if (is.null(n.tau)) n.tau <- n
  tau <- seq(0, 1, length = n.tau)
  if (!is.null(subsample)) {
    quant.o.c <- rowMeans(apply(
      replicate(
        subsample,
        sample(o.c, size = length(tau))
      ),
      2, quantile,
      probs = tau, type = pp.type
    ))
    quant.m.c <- rowMeans(apply(
      replicate(
        subsample,
        sample(m.c, size = length(tau))
      ),
      2, quantile,
      probs = tau, type = pp.type
    ))
    quant.m.p <- rowMeans(apply(
      replicate(
        subsample,
        sample(m.p, size = length(tau))
      ),
      2, quantile,
      probs = tau, type = pp.type
    ))
  } else {
    quant.o.c <- quantile(o.c, tau, type = pp.type)
    quant.m.c <- quantile(m.c, tau, type = pp.type)
    quant.m.p <- quantile(m.p, tau, type = pp.type)
  }
  # Apply quantile delta mapping bias correction
  tau.m.p <- approx(quant.m.p, tau, m.p, rule = 2)$y
  if (ratio) {
    approx.t.qmc.tmp <- approx(tau, quant.m.c, tau.m.p, rule = 2)$y
    delta.m <- m.p / approx.t.qmc.tmp
    delta.m[(delta.m > ratio.max) &
      (approx.t.qmc.tmp < ratio.max.trace)] <- ratio.max
    mhat.p <- approx(tau, quant.o.c, tau.m.p, rule = 2)$y * delta.m
  } else {
    delta.m <- m.p - approx(tau, quant.m.c, tau.m.p, rule = 2)$y
    mhat.p <- approx(tau, quant.o.c, tau.m.p, rule = 2)$y + delta.m
  }
  mhat.c <- approx(quant.m.c, quant.o.c, m.c, rule = 2)$y
  # For ratio data, set values less than trace to zero
  if (ratio) {
    mhat.c[mhat.c < trace] <- 0
    mhat.p[mhat.p < trace] <- 0
  }
  if (ECBC) {
    # empirical copula coupling/Schaake shuffle
    if (length(mhat.p) == length(o.c)) {
      mhat.p <- sort(mhat.p)[rank(o.c, ties.method = ties)]
    } else {
      stop("Schaake shuffle failed due to incompatible lengths")
    }
  }
  list(mhat.c = mhat.c, mhat.p = mhat.p)
}
