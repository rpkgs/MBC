################################################################################
# MBC-QDM.R - Multivariate bias correction based on quantile delta mapping
# and iterative application of Cholesky decomposition rescaling (MBCp and MBCr)
# Multivariate bias correction based on quantile delta mapping and the
# N-dimensional pdf transform (MBCn)
# Alex J. Cannon (alex.cannon@canada.ca)
################################################################################

#' @keywords internal
#' @importFrom FNN knnx.index
#' @importFrom energy edist
#' @importFrom Matrix nearPD
#' @importFrom stats approx cor cov quantile rnorm runif sd
"_PACKAGE"


#' Multivariate Bias Correction of Climate Model Outputs
#' 
#' Calibrate and apply multivariate bias correction algorithms for climate model
#' simulations of multiple climate variables. Three iterative methods are
#' supported: (i) MBC Pearson correlation (MBCp), (ii) MBC rank correlation
#' (MBCr), and (iii) MBC N-dimensional probability density function transform
#' (MBCn). The first two, MBCp and MBCr (Cannon, 2016), match marginal
#' distributions and inter-variable dependence structure. Dependence structure
#' can be measured either by the Pearson correlation ([MBCp()]) or by the
#' Spearman rank correlation ([MBCr()]). The energy distance score ([escore()])
#' is recommended for model selection. The third, [MBCn()] (Cannon, 2018), which
#' operates on the full multivariate distribution, is more flexible and can be
#' considered to be a multivariate analogue of univariate quantile mapping. All
#' aspects of the observed distribution are transferred to the climate model
#' simulations. In each of the three methods, marginal distributions are
#' corrected by the change-preserving quantile delta mapping ([QDM()]) algorithm
#' (Cannon et al., 2015). Finally, an implementation of the Rank Resampling for
#' Distributions and Dependences (R2D2) method introduced by Vrac (2018) is also
#' included.
#' 
#' An example application of the three MBC methods using the `cccma` dataset can
#' be run via:
#' 
#' `example(MBC, run.dontrun=TRUE)`
#' 
#' Note: because empirical quantiles and their changes are used by
#' [QDM()], sample sizes of the observed, model calibration, and
#' model projection datasets should be approximately equal.
#' 
#' \tabular{ll}{ Package: \tab MBC\cr Type: \tab Package\cr License: \tab
#' GPL-2\cr LazyLoad: \tab yes\cr }
#' 
#' @name MBC-package
#' @aliases MBC-package MBC
#' @docType package
#' @seealso `[QDM], [MBCp], [MBCr], [MBCn], [R2D2], [escore], [rot.random], [cccma]`
#' 
#' @references 
#' 1. Cannon, A.J., 2018. Multivariate quantile mapping bias correction: An
#' N-dimensional probability density function transform for climate model
#' simulations of multiple variables. Climate Dynamics, 50(1-2):31-49.
#' doi:10.1007/s00382-017-3580-6
#' 
#' 2. Cannon, A.J., 2016. Multivariate bias correction of climate model output:
#' Matching marginal distributions and inter-variable dependence structure.
#' Journal of Climate, 29:7045-7064. doi:10.1175/JCLI-D-15-0679.1
#' 
#' 3. Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate, 28:6938-6959.
#' doi:10.1175/JCLI-D-14-00754.1
#' 
#' 4. Francois, B., M. Vrac, A.J. Cannon, Y. Robin, and D. Allard, 2020.
#' Multivariate bias corrections of climate simulations: Which benefits for
#' which losses? Earth System Dynamics, 11:537-562. doi:10.5194/esd-11-537-2020
#' 
#' 5. Vrac, M., 2018. Multivariate bias adjustment of high-dimensional climate
#' simulations: the Rank Resampling for Distributions and Dependences (R2D2)
#' bias correction. Hydrology and Earth System Sciences, 22:3175-3196.
#' doi:10.5194/hess-22-3175-2018 
#' 
#' @keywords package 
#' @example R/examples/ex-all.R
NULL
