################################################################################
# MBC-QDM.R - Multivariate bias correction based on quantile delta mapping
# and iterative application of Cholesky decomposition rescaling (MBCp and MBCr)
# Multivariate bias correction based on quantile delta mapping and the 
# N-dimensional pdf transform (MBCn)
# Alex J. Cannon (alex.cannon@canada.ca)
################################################################################

# library(Matrix)
# library(energy)
# library(FNN)

# Quantile delta mapping bias correction for preserving changes in quantiles
# Note: QDM is equivalent to the equidistant and equiratio forms of quantile
# mapping (Cannon et al., 2015).

# Cannon, A.J., Sobie, S.R., and Murdock, T.Q. 2015. Bias correction of
#  simulated precipitation by quantile mapping: How well do methods preserve
#  relative changes in quantiles and extremes? Journal of Climate,
#  28: 6938-6959. doi:10.1175/JCLI-D-14-00754.1

# o = vector of observed values; m = vector of modelled values
# c = current period;  p = projected period
# ratio = TRUE --> preserve relative trends in a ratio variable
# trace = 0.05 --> replace values less than trace with exact zeros
# trace.calc = 0.5*trace --> treat values below trace.calc as censored
# jitter.factor = 0.01 --> jitter to accommodate ties
# n.tau = NULL --> number of empirical quantiles (NULL=sample length)
# ratio.max = 2 --> maximum delta when values are less than ratio.max.trace
# ratio.max.trace = 10*trace --> values below which ratio.max is applied
# ECBC = TRUE --> apply Schaake shuffle to enforce o.c temporal sequencing
# subsample = NULL --> use this number of repeated subsamples of size n.tau
#  to calculate empirical quantiles (e.g., when o.c, m.c, and m.p are of
#  very different size)
# pp.type = 7 --> plotting position type used in quantile
# tau.m-p = F.m-p(x.m-p)
# delta.m = x.m-p {/,-} F.m-c^-1(tau.m-p)
# xhat.m-p = F.o-c^-1(tau.m-p) {*,+} delta.m



################################################################################
# Multivariate goodness-of-fit scoring function



#' Energy distance score
#' 
#' Calculate the energy distance score measuring the statistical discrepancy
#' between samples `x` and `y` from two multivariate distributions.
#' 
#' 
#' @param x numeric matrix.
#' @param y numeric matrix.
#' @param scale.x logical indicating whether data should be standardized based
#' on `x`.
#' @param n.cases the number of sub-sampled cases; `NULL` uses all data.
#' @param alpha distance exponent in (0,2]
#' @param method method used to weight the statistics
#' @references Székely, G.J. and M.L. Rizzo, 2004. Testing for equal
#' distributions in high dimension, InterStat, November (5).
#' 
#' Székely, G.J. and M.L. Rizzo, 2013. Energy statistics: statistics based on
#' distances. Journal of Statistical Planning and Inference, 143(8):1249-1272.
#' doi:10.1016/j.jspi.2013.03.018
#' 
#' Rizzo, M.L. and G.L. Székely, 2016. Energy distance. Wiley Interdisciplinary
#' Reviews: Computational Statistics, 8(1):27-38.
#' @export escore
escore <- function (x, y, scale.x = FALSE, n.cases = NULL, alpha = 1, method = "cluster") {
    n.x <- nrow(x)
    n.y <- nrow(y)
    if (scale.x) {
        x <- scale(x)
        y <- scale(y, center = attr(x, "scaled:center"), scale = attr(x,
            "scaled:scale"))
    }
    if (!is.null(n.cases)) {
        n.cases <- min(n.x, n.y, n.cases)
        x <- x[sample(n.x, size = n.cases), , drop = FALSE]
        y <- y[sample(n.y, size = n.cases), , drop = FALSE]
        n.x <- n.cases
        n.y <- n.cases
    }
    edist(rbind(x, y), sizes = c(n.x, n.y), distance = FALSE,
        alpha = alpha, method = method)[1]/2
}

################################################################################
# Multivariate bias correction based on iterative application of quantile
# mapping (or ranking) and multivariate rescaling via Cholesky
# decomposition of the covariance matrix. Results in simulated marginal 
# distributions and Pearson or Spearman rank correlations that match
# observations
# Cannon, A.J. 2016. Multivariate Bias Correction of Climate Model Outputs:
#  Matching Marginal Distributions and Inter-variable Dependence Structure.
#  Journal of Climate, doi:10.1175/JCLI-D-15-0679.1



#' Multivariate linear rescaling using Cholesky decomposition
#' 
#' Multivariate linear bias correction based on Cholesky decomposition of the
#' covariance matrix following Scheuer and Stoller (1962) and Bürger et al.
#' (2011). Bias correction matches the multivariate mean and covariance
#' structure.
#' 
#' 
#' @param o.c matrix of observed samples during the calibration period.
#' @param m.c matrix of model outputs during the calibration period.
#' @param m.p matrix of model outputs during the projected period.
#' @param o.c.chol precalculated Cholesky decomposition of the `o.c`
#' covariance matrix; `NULL` calculates internally.
#' @param o.p.chol precalculated Cholesky decomposition of the target
#' `o.p` covariance matrix; `NULL` defaults to `o.c.chol`.
#' @param m.c.chol precalculated Cholesky decomposition of the `m.c`
#' covariance matrix; `NULL` calculates internally.
#' @param m.p.chol precalculated Cholesky decomposition of the `m.p`
#' covariance matrix; `NULL` calculates internally.
#' @return a list of with elements consisting of: \item{mhat.c}{matrix of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{matrix of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[MBCp], [MBCr]`
#' @references Scheuer, E.M. and D.S. Stoller, 1962. On the generation of
#' normal random vectors. Technometrics, 4(2):278-281.
#' 
#' Bürger, G., J. Schulla, and A.T. Werner, 2011. Estimates of future flow,
#' including extremes, of the Columbia River headwaters. Water Resources
#' Research, 47(10):W10520. doi:10.1029/2010WR009716
#' @export MRS
MRS <- function(o.c, m.c, m.p, o.c.chol=NULL, o.p.chol=NULL, m.c.chol=NULL,
         m.p.chol=NULL){
    # Center based on multivariate means
    o.c.mean <- colMeans(o.c)
    m.c.mean <- colMeans(m.c)
    m.p.mean <- colMeans(m.p)
    o.c <- sweep(o.c, 2, o.c.mean, '-')
    m.c <- sweep(m.c, 2, m.c.mean, '-')
    m.p <- sweep(m.p, 2, m.p.mean, '-')
    # Cholesky decomposition of covariance matrix
    # If !is.null(o.p.chol) --> projected target
    if(is.null(o.c.chol)) o.c.chol <- chol(cov(o.c))
    if(is.null(o.p.chol)) o.p.chol <- chol(cov(o.c))
    if(is.null(m.c.chol)) m.c.chol <- chol(cov(m.c))
    if(is.null(m.p.chol)) m.p.chol <- chol(cov(m.c))
    # Bias correction factors
    mbcfactor <- solve(m.c.chol) %*% o.c.chol
    mbpfactor <- solve(m.p.chol) %*% o.p.chol
    # Multivariate bias correction
    mbc.c <- m.c %*% mbcfactor
    mbc.p <- m.p %*% mbpfactor
    # Recenter and account for change in means
    mbc.c <- sweep(mbc.c, 2, o.c.mean, '+')
    mbc.p <- sweep(mbc.p, 2, o.c.mean, '+')
    mbc.p <- sweep(mbc.p, 2, m.p.mean-m.c.mean, '+')
    list(mhat.c=mbc.c, mhat.p=mbc.p)
}



#' Multivariate bias correction (Spearman rank correlation)
#' 
#' Multivariate bias correction that matches marginal distributions using
#' [QDM()] and the Spearman rank correlation dependence structure
#' following Cannon (2016).
#' 
#' 
#' @param o.c matrix of observed samples during the calibration period.
#' @param m.c matrix of model outputs during the calibration period.
#' @param m.p matrix of model outputs during the projected period.
#' @param iter maximum number of algorithm iterations.
#' @param cor.thresh if greater than zero, a threshold indicating the change in
#' magnitude of Spearman rank correlations required for convergence.
#' @param ratio.seq vector of logical values indicating if samples are of a
#' ratio quantity (e.g., precipitation).
#' @param trace numeric values indicating thresholds below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric values of thresholds used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric values indicating the maximum proportional changes
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric values of trace thresholds used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param qmap.precalc logical value indicating if `m.c` and `m.p`
#' are outputs from `QDM`.
#' @param silent logical value indicating if algorithm progress should be
#' reported.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' @return a list of with elements consisting of: \item{mhat.c}{matrix of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{matrix of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[QDM], [MBCp], [MRS], [MBCn]
#' [escore]`
#' @references Cannon, A.J., 2016. Multivariate bias correction of climate
#' model output: Matching marginal distributions and inter-variable dependence
#' structure. Journal of Climate, 29:7045-7064. doi:10.1175/JCLI-D-15-0679.1
#' 
#' Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate,
#' 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
#' @export MBCr
MBCr <-
# Multivariate quantile mapping bias correction (Spearman correlation)
# Cannon, A.J., 2016. Multivariate bias correction of climate model outputs: 
#  matching marginal distributions and inter-variable dependence structure.
#  Journal of Climate, 29(19):7045–7064. doi:10.1175/JCLI-D-15-0679.1
function(o.c, m.c, m.p, iter=20, cor.thresh=1e-4,
         ratio.seq=rep(FALSE, ncol(o.c)), trace=0.05,
         trace.calc=0.5*trace, jitter.factor=0, n.tau=NULL, ratio.max=2,
         ratio.max.trace=10*trace, ties='first', qmap.precalc=FALSE,
         silent=FALSE, subsample=NULL, pp.type=7){
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max) == 1)
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    # Ordinal ranks of observed and modelled data
    o.c.r <- apply(o.c, 2, rank, ties.method=ties)
    m.c.r <- apply(m.c, 2, rank, ties.method=ties)
    m.p.r <- apply(m.p, 2, rank, ties.method=ties)
    m.c.i <- m.c.r
    if(cor.thresh > 0){
        # Spearman correlation to assess convergence
        cor.i <- cor(m.c.r)
        cor.i[is.na(cor.i)] <- 0
    } else{
        cor.diff <- 0
    }
    # Iterative MBC/reranking
    o.c.chol <- o.p.chol <- as.matrix(chol(nearPD(cov(o.c.r))$mat))
    for(i in seq(iter)){
        m.c.chol <- m.p.chol <- as.matrix(chol(nearPD(cov(m.c.r))$mat))
        fit.mbc <- MRS(o.c=o.c.r, m.c=m.c.r, m.p=m.p.r, o.c.chol=o.c.chol,
                       o.p.chol=o.p.chol, m.c.chol=m.c.chol, m.p.chol=m.p.chol)
        m.c.r <- apply(fit.mbc$mhat.c, 2, rank, ties.method=ties)
        m.p.r <- apply(fit.mbc$mhat.p, 2, rank, ties.method=ties)
        if(cor.thresh > 0){
            # Check on Spearman correlation convergence
            cor.j <- cor(m.c.r)
            cor.j[is.na(cor.j)] <- 0
            cor.diff <- mean(abs(cor.j-cor.i))
            cor.i <- cor.j
        }
        if(!silent){
            cat(i, mean(m.c.r==m.c.i), cor.diff, '')
        }
        if(cor.diff < cor.thresh) break
        if(identical(m.c.r, m.c.i)) break
        m.c.i <- m.c.r
    }
    if(!silent) cat('\n')
    for(i in seq(ncol(o.c))){
        # Replace ordinal ranks with QDM outputs
        m.c.r[,i] <- sort(m.c.qmap[,i])[m.c.r[,i]]
        m.p.r[,i] <- sort(m.p.qmap[,i])[m.p.r[,i]]
    }
    list(mhat.c=m.c.r, mhat.p=m.p.r)
}



#' Multivariate bias correction (Pearson correlation)
#' 
#' Multivariate bias correction that matches marginal distributions using
#' [QDM()] and the Pearson correlation dependence structure following
#' Cannon (2016).
#' 
#' 
#' @param o.c matrix of observed samples during the calibration period.
#' @param m.c matrix of model outputs during the calibration period.
#' @param m.p matrix of model outputs during the projected period.
#' @param iter maximum number of algorithm iterations.
#' @param cor.thresh if greater than zero, a threshold indicating the change in
#' magnitude of Pearson correlations required for convergence.
#' @param ratio.seq vector of logical values indicating if samples are of a
#' ratio quantity (e.g., precipitation).
#' @param trace numeric values indicating thresholds below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric values of thresholds used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric values indicating the maximum proportional changes
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric values of trace thresholds used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param qmap.precalc logical value indicating if `m.c` and `m.p`
#' are outputs from `QDM`.
#' @param silent logical value indicating if algorithm progress should be
#' reported.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate initial empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' @return a list of with elements consisting of: \item{mhat.c}{matrix of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{matrix of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[QDM], [MBCr], [MRS], [MBCn]
#' [escore]`
#' @references Cannon, A.J., 2016. Multivariate bias correction of climate
#' model output: Matching marginal distributions and inter-variable dependence
#' structure. Journal of Climate, 29:7045-7064. doi:10.1175/JCLI-D-15-0679.1
#' 
#' Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate,
#' 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
#' @export MBCp
MBCp <-
# Multivariate quantile mapping bias correction (Pearson correlation)
# Cannon, A.J., 2016. Multivariate bias correction of climate model outputs: 
#  matching marginal distributions and inter-variable dependence structure.
#  Journal of Climate, 29(19):7045–7064. doi:10.1175/JCLI-D-15-0679.1
function(o.c, m.c, m.p, iter=20, cor.thresh=1e-4,
         ratio.seq=rep(FALSE, ncol(o.c)), trace=0.05, trace.calc=0.5*trace,
         jitter.factor=0, n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
         ties='first', qmap.precalc=FALSE, silent=FALSE, subsample=NULL,
         pp.type=7){
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max) == 1)
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    m.c <- m.c.qmap
    m.p <- m.p.qmap
    # Pearson correlation to assess convergence
    if(cor.thresh > 0){
        cor.i <- cor(m.c)
        cor.i[is.na(cor.i)] <- 0
    }
    o.c.chol <- o.p.chol <- as.matrix(chol(nearPD(cov(o.c))$mat))
    # Iterative MBC/QDM
    for(i in seq(iter)){
        m.c.chol <- m.p.chol <- as.matrix(chol(nearPD(cov(m.c))$mat))
        fit.mbc <- MRS(o.c=o.c, m.c=m.c, m.p=m.p, o.c.chol=o.c.chol,
                       o.p.chol=o.p.chol, m.c.chol=m.c.chol, m.p.chol=m.p.chol)
        m.c <- fit.mbc$mhat.c
        m.p <- fit.mbc$mhat.p
        for(j in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,j], m.c=m.c[,j], m.p=m.p[,j], ratio=FALSE,
                            n.tau=n.tau, pp.type=pp.type)
            m.c[,j] <- fit.qmap$mhat.c
            m.p[,j] <- fit.qmap$mhat.p
        }
        # Check on Pearson correlation convergence
        if(cor.thresh > 0){
            cor.j <- cor(m.c)
            cor.j[is.na(cor.j)] <- 0
            cor.diff <- mean(abs(cor.j-cor.i))
            cor.i <- cor.j
        } else{
            cor.diff <- 0
        }
        if(!silent) cat(i, cor.diff, '')
        if(cor.diff < cor.thresh) break
    }
    if(!silent) cat('\n')
    # Replace with shuffled QDM elements
    for(i in seq(ncol(o.c))){
        m.c[,i] <- sort(m.c.qmap[,i])[rank(m.c[,i], ties.method=ties)]
        m.p[,i] <- sort(m.p.qmap[,i])[rank(m.p[,i], ties.method=ties)]
    }
    list(mhat.c=m.c, mhat.p=m.p)
}

################################################################################
# Multivariate bias correction based on iterative application of random
# orthogonal rotation and quantile mapping (N-dimensional pdf transfer)
# Pitié, F., Kokaram, A.C., and Dahyot, R. 2005. N-dimensional probability
#  density function transfer and its application to color transfer.
#  In Tenth IEEE International Conference on Computer Vision, 2005. ICCV 2005.
#  (Vol. 2, pp. 1434-1439). IEEE.
# Pitié, F., Kokaram, A.C., and Dahyot, R. 2007. Automated colour grading
#  using colour distribution transfer. Computer Vision and Image Understanding,
#  107(1), 123-137.

# Random orthogonal rotation


#' Random orthogonal rotation
#' 
#' Generate a `k`-dimensional random orthogonal rotation matrix.
#' 
#' 
#' @param k the number of dimensions.
#' @references Mezzadri, F. 2007. How to generate random matrices from the
#' classical compact groups, Notices of the American Mathematical Society,
#' 54:592–604.
#' @export rot.random
rot.random <- function(k) {
  rand <- matrix(rnorm(k * k), ncol=k)
  QRd <- qr(rand)
  Q <- qr.Q(QRd)
  R <- qr.R(QRd)
  diagR <- diag(R)
  rot <- Q %*% diag(diagR/abs(diagR))
  return(rot)
}

# Multivariate quantile mapping bias correction (N-dimensional pdf transfer)
# Cannon, A.J., 2018. Multivariate quantile mapping bias correction: An 
#  N-dimensional probability density function transform for climate model
#  simulations of multiple variables. Climate Dynamics, 50(1-2):31-49.
#  doi:10.1007/s00382-017-3580-6


#' Multivariate bias correction (N-pdft)
#' 
#' Multivariate bias correction that matches the multivariate distribution
#' using [QDM()] and the N-dimensional probability density function
#' transform (N-pdft) following Cannon (2018).
#' 
#' 
#' @param o.c matrix of observed samples during the calibration period.
#' @param m.c matrix of model outputs during the calibration period.
#' @param m.p matrix of model outputs during the projected period.
#' @param iter maximum number of algorithm iterations.
#' @param ratio.seq vector of logical values indicating if samples are of a
#' ratio quantity (e.g., precipitation).
#' @param trace numeric values indicating thresholds below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric values of thresholds used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric values indicating the maximum proportional changes
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric values of trace thresholds used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param qmap.precalc logical value indicating if `m.c` and `m.p`
#' are outputs from `QDM`.
#' @param rot.seq use a supplied list of random rotation matrices. `NULL`
#' generates on the fly.
#' @param silent logical value indicating if algorithm progress should be
#' reported.
#' @param n.escore number of cases used to calculate the energy distance when
#' monitoring convergence.
#' @param return.all logical value indicating whether results from all
#' iterations are returned.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate initial empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' @return a list of with elements consisting of: \item{mhat.c}{matrix of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{matrix of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[QDM], [MBCp], [MBCr], [MRS],
#' [escore], [rot.random]`
#' @references Cannon, A.J., 2018. Multivariate quantile mapping bias
#' correction: An N-dimensional probability density function transform for
#' climate model simulations of multiple variables. Climate Dynamics,
#' 50(1-2):31-49. doi:10.1007/s00382-017-3580-6
#' 
#' Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate,
#' 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
#' 
#' Pitié, F., A.C. Kokaram, and R. Dahyot, 2005. N-dimensional probability
#' density function transfer and its application to color transfer. In Tenth
#' IEEE International Conference on Computer Vision, 2005. ICCV 2005. (Vol. 2,
#' pp. 1434-1439). IEEE.
#' 
#' Pitié, F., A.C. Kokaram, and R. Dahyot, 2007. Automated colour grading using
#' colour distribution transfer. Computer Vision and Image Understanding,
#' 107(1):123-137.
#' @export MBCn
MBCn <- function(o.c, m.c, m.p, iter=30, ratio.seq=rep(FALSE, ncol(o.c)),
         trace=0.05, trace.calc=0.5*trace, jitter.factor=0, n.tau=NULL,
         ratio.max=2, ratio.max.trace=10*trace, ties='first',
         qmap.precalc=FALSE, rot.seq=NULL, silent=FALSE, n.escore=0,
         return.all=FALSE, subsample=NULL, pp.type=7){
    if(!is.null(rot.seq)){
        if(length(rot.seq)!=iter){
            stop('length(rot.seq) != iter')
        }
    }
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max) == 1)
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    # Energy score (rescaled)
    escore.iter <- rep(NA, iter+2)
    if(n.escore > 0){
        n.escore <- min(nrow(o.c), nrow(m.c), n.escore)
        escore.cases.o.c <- unique(suppressWarnings(matrix(seq(nrow(o.c)),
                                   ncol=n.escore)[1,]))
        escore.cases.m.c <- unique(suppressWarnings(matrix(seq(nrow(m.c)),
                                   ncol=n.escore)[1,]))
        escore.iter[1] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                 y=m.c[escore.cases.m.c,,drop=FALSE],
                                 scale.x=TRUE)
        if(!silent) cat('RAW', escore.iter[1], ': ')        
    }
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    m.c <- m.c.qmap
    m.p <- m.p.qmap
    # Energy score (QDM)
    if(n.escore > 0){
        escore.iter[2] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                 y=m.c[escore.cases.m.c,,drop=FALSE],
                                 scale.x=TRUE)
        if(!silent) cat('QDM', escore.iter[2], ': ')        
    }
    # Standardize observations
    m.iter <- vector('list', iter)
    o.c.mean <- colMeans(o.c)
    o.c.sdev <- apply(o.c, 2, sd)
    o.c.sdev[o.c.sdev < .Machine$double.eps] <- 1    
    o.c <- scale(o.c, center=o.c.mean, scale=o.c.sdev)
    # Standardize model
    m.c.p <- rbind(m.c, m.p)
    m.c.p.mean <- colMeans(m.c.p)
    m.c.p.sdev <- apply(m.c.p, 2, sd)
    m.c.p.sdev[m.c.p.sdev < .Machine$double.eps] <- 1
    m.c.p <- scale(m.c.p, center=m.c.p.mean, scale=m.c.p.sdev)    
    Xt <- rbind(o.c, m.c.p)
    for(i in seq(iter)){
        if(!silent) cat(i, '')
        # Random orthogonal rotation
        if(is.null(rot.seq)){
            rot <- rot.random(ncol(o.c))
        } else{
            rot <- rot.seq[[i]]
        }
        Z <- Xt %*% rot
        Z.o.c <- Z[1:nrow(o.c),,drop=FALSE]
        Z.m.c <- Z[(nrow(o.c)+1):(nrow(o.c)+nrow(m.c)),,drop=FALSE]
        Z.m.p <- Z[(nrow(o.c)+nrow(m.c)+1):nrow(Z),,drop=FALSE]
        # Bias correct rotated variables using QDM
        for(j in seq(ncol(Z))){
            Z.qdm <- QDM(o.c=Z.o.c[,j], m.c=Z.m.c[,j], m.p=Z.m.p[,j],
                         ratio=FALSE, jitter.factor=jitter.factor[j],
                         n.tau=n.tau, pp.type=pp.type)
            Z.m.c[,j] <- Z.qdm$mhat.c
            Z.m.p[,j] <- Z.qdm$mhat.p
        }
        # Rotate back
        m.c <- Z.m.c %*% t(rot)
        m.p <- Z.m.p %*% t(rot)
        Xt <- rbind(o.c, m.c, m.p)
        # Energy score (MBCn)
        if(n.escore > 0){
            escore.iter[i+2] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                       y=m.c[escore.cases.m.c,,drop=FALSE],
                                       scale.x=TRUE)
            if(!silent) cat(escore.iter[i+2], ': ')
        }
        if(return.all){
            m.c.i <- sweep(sweep(m.c, 2, attr(m.c.p, 'scaled:scale'), '*'), 2,
                           attr(m.c.p, 'scaled:center'), '+')
            m.p.i <- sweep(sweep(m.p, 2, attr(m.c.p, 'scaled:scale'), '*'), 2,
                           attr(m.c.p, 'scaled:center'), '+')
            m.iter[[i]] <- list(m.c=m.c.i, m.p=m.p.i)
        }
    }
    if(!silent) cat('\n')
    # Rescale back to original units
    m.c <- sweep(sweep(m.c, 2, m.c.p.sdev, '*'), 2, m.c.p.mean, '+')
    m.p <- sweep(sweep(m.p, 2, m.c.p.sdev, '*'), 2, m.c.p.mean, '+')
    # Replace npdft ordinal ranks with QDM outputs
    for(i in seq(ncol(o.c))){
        m.c[,i] <- sort(m.c.qmap[,i])[rank(m.c[,i], ties.method=ties)]
        m.p[,i] <- sort(m.p.qmap[,i])[rank(m.p[,i], ties.method=ties)]
    }
    names(escore.iter)[1:2] <- c('RAW', 'QM')
    names(escore.iter)[-c(1:2)] <- seq(iter)
    list(mhat.c=m.c, mhat.p=m.p, escore.iter=escore.iter, m.iter=m.iter)
}

################################################################################
# Multivariate bias correction based on application of the nearest
# neighbour algorithm to ordinal ranks.
# Vrac, M., 2018. Multivariate bias adjustment of high-dimensional climate
#   simulations: the Rank Resampling for Distributions and Dependences (R2D2)
#   bias correction. Hydrology and Earth System Sciences, 22:3175-3196.
#   doi:10.5194/hess-22-3175-2018


# Vrac et al., (2018)


#' Multivariate bias correction (R2D2)
#' 
#' Multivariate bias correction that matches the multivariate distribution
#' using [QDM()] and the R2D2 algorithm following Vrac (2018).
#' 
#' @param o.c matrix of observed samples during the calibration period.
#' @param m.c matrix of model outputs during the calibration period.
#' @param m.p matrix of model outputs during the projected period.
#' @param ref.column index of the reference column used for the 1D nearest
#' neighbour matching
#' @param ratio.seq vector of logical values indicating if samples are of a
#' ratio quantity (e.g., precipitation).
#' @param trace numeric values indicating thresholds below which values of a
#' ratio quantity (e.g., `ratio=TRUE`) should be considered exact zeros.
#' @param trace.calc numeric values of thresholds used internally when handling
#' of exact zeros; defaults to one half of `trace`.
#' @param jitter.factor optional strength of jittering to be applied when
#' quantities are quantized.
#' @param n.tau number of quantiles used in the quantile mapping; `NULL`
#' equals the length of the `m.p` series.
#' @param ratio.max numeric values indicating the maximum proportional changes
#' allowed for ratio quantities below the `ratio.max.trace` threshold.
#' @param ratio.max.trace numeric values of trace thresholds used to constrain
#' the proportional change in ratio quantities to `ratio.max`; defaults to
#' ten times `trace`.
#' @param ties method used to handle ties when calculating ordinal ranks.
#' @param qmap.precalc logical value indicating if `m.c` and `m.p`
#' are outputs from `QDM`.
#' @param subsample use `subsample` draws of size `n.tau` to
#' calculate initial empirical quantiles; if `NULL`, calculate normally.
#' @param pp.type type of plotting position used in `quantile`.
#' @return a list of with elements consisting of: \item{mhat.c}{matrix of bias
#' corrected `m.c` values for the calibration period.}
#' \item{mhat.p}{matrix of bias corrected `m.p` values for the projection
#' period.}
#' @seealso `[QDM], [MBCp], [MBCr], [MRS],
#' [MBCn]`
#' @references Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias
#' correction of simulated precipitation by quantile mapping: How well do
#' methods preserve relative changes in quantiles and extremes? Journal of
#' Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
#' 
#' Vrac, M., 2018. Multivariate bias adjustment of high-dimensional climate
#' simulations: the Rank Resampling for Distributions and Dependences (R2D2)
#' bias correction. Hydrology and Earth System Sciences, 22:3175-3196.
#' doi:10.5194/hess-22-3175-2018
#' @export R2D2
R2D2 <- function(o.c, m.c, m.p, ref.column = 1, ratio.seq = rep(FALSE,
    ncol(o.c)), trace = 0.05, trace.calc = 0.5 * trace, jitter.factor = 0,
    n.tau = NULL, ratio.max = 2, ratio.max.trace = 10 * trace,
    ties = "first", qmap.precalc = FALSE, subsample = NULL,
    pp.type = 7)
{
    if ((length(o.c) != length(m.c)) || (length(o.c) != length(m.p))){
        stop("R2D2 requires data samples of equal length")
    }
    if (length(trace.calc) == 1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if (length(trace) == 1)
        trace <- rep(trace, ncol(o.c))
    if (length(jitter.factor) == 1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if (length(ratio.max) == 1)
        ratio.max <- rep(ratio.max, ncol(o.c))
    if (length(ratio.max.trace) == 1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if (!qmap.precalc) {
        for (i in seq(ncol(o.c))) {
            fit.qmap <- QDM(o.c = o.c[, i], m.c = m.c[, i], m.p = m.p[,
                i], ratio = ratio.seq[i], trace.calc = trace.calc[i],
                trace = trace[i], jitter.factor = jitter.factor[i],
                n.tau = n.tau, ratio.max = ratio.max[i],
                ratio.max.trace = ratio.max.trace[i],
                subsample = subsample, pp.type = pp.type)
            m.c.qmap[, i] <- fit.qmap$mhat.c
            m.p.qmap[, i] <- fit.qmap$mhat.p
        }
    }
    # Calculate ordinal ranks of observations and m.c.qmap and m.p.qmap
    o.c.r <- apply(o.c, 2, rank, ties.method = ties)
    m.c.r <- apply(m.c.qmap, 2, rank, ties.method = ties)
    m.p.r <- apply(m.p.qmap, 2, rank, ties.method = ties)
    # 1D rank analog selection based on ref.column
    nn.c.r <- rank(knnx.index(o.c.r[,ref.column],
                   query=m.c.r[,ref.column], k=1), ties.method='random')
    nn.p.r <- rank(knnx.index(o.c.r[,ref.column],
                   query=m.p.r[,ref.column], k=1), ties.method='random')
    # Shuffle o.c.r ranks based on 1D rank analogs
    new.c.r <- o.c.r[nn.c.r,,drop=FALSE]
    new.p.r <- o.c.r[nn.p.r,,drop=FALSE]
    # Reorder m.c.qmap and m.p.qmap
    r2d2.c <- m.c.qmap
    r2d2.p <- m.p.qmap
    for (i in seq(ncol(o.c))) {
        r2d2.c[,i] <- sort(r2d2.c[,i])[new.c.r[,i]]
        r2d2.p[,i] <- sort(r2d2.p[,i])[new.p.r[,i]]
    }    
    list(mhat.c = r2d2.c, mhat.p = r2d2.p)
}
