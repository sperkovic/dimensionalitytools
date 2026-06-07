############
# Package documentation
############

#' bifactorTools: Efficient Bifactor Model Indices
#'
#' Computes bifactor model indices from fitted lavaan or mirt objects, or from
#' a raw loading matrix. Designed for efficiency in simulation contexts (no model
#' refitting) while remaining accessible for single-dataset applied use.
#'
#' @section Main function:
#' \code{\link{bifactor_indices}} is the primary entry point. It accepts a fitted
#' model or loading matrix, extracts the loading structure once, and computes all
#' indices from that extracted matrix.
#'
#' @section Index set:
#' \itemize{
#'   \item ECV (Explained Common Variance) - general and per subscale
#'   \item PUC (Percent Uncontaminated Correlations)
#'   \item Omega hierarchical, subscale, and total
#'   \item H (Construct Replicability) - general and per subscale
#'   \item FD (Factor Determinacy) - per factor
#'   \item Cronbach's alpha - total and per subscale
#'   \item ARPB (Average Relative Parameter Bias) - optional, requires unidimensional model
#' }
#'
#' @section Note on H computation:
#' H is computed using total communality residuals (psi_i = 1 - sum_f lambda_fi^2),
#' consistent with Reise (2012). This differs from Dueber (2017), who uses per-factor
#' residuals. The total communality approach is more conservative and appropriate when
#' working from the full standardized loading matrix.
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based
#' tool to calculate various indices relevant to bifactor CFA models.
#' doi:10.13023/edp.tool.01
#'
#' Grice, J. W. (2001). Computing and evaluating factor scores. Psychological
#' Methods, 6(4), 430-450.
#'
#' Reise, S. P. (2012). The rediscovery of bifactor measurement models.
#' Multivariate Behavioral Research, 47(5), 667-696.
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @docType package
#' @name bifactorTools
"_PACKAGE"
