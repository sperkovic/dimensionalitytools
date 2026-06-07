############
# Main dispatcher
# Single entry point for all bifactor indices
############


#' Compute bifactor model indices
#'
#' @param object A fitted lavaan object, mirt SingleGroupClass object, or a numeric
#'   matrix of standardized factor loadings (items x factors, general factor first column).
#' @param uni_object Optional. A second fitted model or loading vector for ARPB computation.
#'   Can be a lavaan object, mirt object, or numeric vector of unidimensional loadings.
#'   Default NULL skips ARPB.
#' @param output Character. Either "list" (default) or "dataframe". The dataframe format
#'   returns a single-row data frame with named columns, suitable for direct use in
#'   SimDesign analyse functions.
#' @param standardized Logical. If TRUE (default), uses standardized loadings.
#'   If FALSE, uses unstandardized loadings for omega, H, FD, and alpha. ECV is
#'   always computed from the standardized solution where possible; a warning is
#'   issued if standardized loadings are unavailable.
#' @param gen_col Integer or character. Column index or name of the general factor
#'   in the loading matrix. Default 1.
#' @param general_only_items Integer vector of item row indices that load only on the
#'   general factor. Used for PUC computation. Default NULL assumes standard bifactor
#'   structure where all items belong to one specific factor (Rodriguez et al., 2016).
#' @param tol Numeric. Threshold below which loadings are treated as zero for
#'   membership assignment and PUC counting. Default 1e-6.
#'
#' @return If output = "list": a named list containing all computed indices.
#'   If output = "dataframe": a single-row data frame with named columns.
#'   Subscale-level indices are stored as named vectors in list output, and
#'   expanded to separate columns (named e.g. ECV_subscale_SF1) in dataframe output.
#'
#' @details
#' All indices are computed from a loading matrix extracted once from the input object.
#' No model refitting occurs. This makes the function efficient in simulation loops
#' where models are already fitted.
#'
#' Index set computed:
#' \itemize{
#'   \item ECV_general, ECV_subscale: Explained common variance (Rodriguez et al., 2016)
#'   \item PUC: Percent uncontaminated correlations (Rodriguez et al., 2016)
#'   \item omega_h: Omega hierarchical
#'   \item omega_s: Omega subscale (per specific factor)
#'   \item omega_t: Omega total
#'   \item H_general, H_subscale: Construct replicability (Reise, 2012)
#'   \item FD: Factor determinacy (Grice, 2001)
#'   \item alpha_total, alpha_subscale: Cronbach's alpha
#'   \item ARPB, SRPB, ARPB_items: Average relative parameter bias (if uni_object supplied)
#' }
#'
#' H is computed using total communality residuals (psi_i = 1 - sum_f lambda_fi^2),
#' consistent with Reise (2012). This differs from Dueber (2017), who uses per-factor
#' residuals.
#'
#' @examples
#' # Raw matrix input (fastest for simulation)
#' Lambda <- matrix(c(.82, .10,  0,   0,
#'                    .77, .35,  0,   0,
#'                    .79, .32,  0,   0,
#'                    .66, .39,  0,   0,
#'                    .51,  0,  .71,  0,
#'                    .56,  0,  .43,  0,
#'                    .68,  0,  .13,  0,
#'                    .60,  0,  .50,  0,
#'                    .83,  0,   0,  .47,
#'                    .60,  0,   0,  .27,
#'                    .78,  0,   0,  .28,
#'                    .55,  0,   0,  .75),
#'                  ncol = 4, byrow = TRUE)
#' colnames(Lambda) <- c("G", "SF1", "SF2", "SF3")
#'
#' # List output
#' bifactor_indices(Lambda)
#'
#' # Dataframe output (SimDesign-friendly)
#' bifactor_indices(Lambda, output = "dataframe")
#'
#' @references
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
#' @export
bifactor_indices <- function(object,
                             uni_object         = NULL,
                             output             = c("list", "dataframe"),
                             standardized       = TRUE,
                             gen_col            = 1,
                             general_only_items = NULL,
                             tol                = 1e-6) {

  output <- match.arg(output)

  ############
  # 1. Extract loading matrix
  ############
  Lambda <- extract_loadings(object, standardized = standardized)
  gen_col <- identify_general(Lambda, gen_col)

  ############
  # 2. Extract unique variances
  ############
  if (is.matrix(object) || is.data.frame(object)) {
    Psi <- extract_psi(Lambda, object = NULL, standardized = standardized)
  } else {
    Psi <- extract_psi(Lambda, object = object, standardized = standardized)
  }

  ############
  # 3. Subscale membership (computed once, reused across all indices)
  ############
  membership <- get_subscale_membership(Lambda, gen_col = gen_col, tol = tol)

  ############
  # 4. Compute all indices
  ############
  ecv    <- compute_ecv(Lambda, gen_col = gen_col, membership = membership)
  puc    <- compute_puc(Lambda, gen_col = gen_col,
                        general_only_items = general_only_items, tol = tol)
  omega  <- compute_omega(Lambda, Psi = Psi, gen_col = gen_col, membership = membership)
  H      <- compute_H(Lambda, gen_col = gen_col)
  fd     <- compute_fd(Lambda, Psi = Psi, gen_col = gen_col)
  alpha  <- compute_alpha(Lambda, Psi = Psi, gen_col = gen_col, membership = membership)

  ############
  # 5. ARPB if uni_object supplied
  ############
  arpb_results <- NULL
  if (!is.null(uni_object)) {
    if (is.numeric(uni_object) && !is.matrix(uni_object)) {
      uni_lam <- uni_object
    } else {
      uni_Lambda <- extract_loadings(uni_object, standardized = standardized)
      uni_lam    <- uni_Lambda[, 1]
    }
    arpb_results <- compute_arpb(Lambda, uni_lambda = uni_lam, gen_col = gen_col)
  }

  ############
  # 6. Assemble results
  ############
  results <- list(
    ECV_general    = ecv$ECV_general,
    ECV_subscale   = ecv$ECV_subscale,
    PUC            = puc,
    omega_h        = omega$omega_h,
    omega_s        = omega$omega_s,
    omega_t        = omega$omega_t,
    H_general      = H$H_general,
    H_subscale     = H$H_subscale,
    FD             = fd,
    alpha_total    = alpha$alpha_total,
    alpha_subscale = alpha$alpha_subscale
  )

  if (!is.null(arpb_results)) {
    results$ARPB       <- arpb_results$ARPB
    results$SRPB       <- arpb_results$SRPB
    results$ARPB_items <- arpb_results$ARPB_items
  }

  if (output == "list") {
    class(results) <- "bifactor_indices"
    return(results)
  }

  ############
  # 7. Flatten to dataframe if requested
  ############
  .to_dataframe(results)
}


#' Flatten bifactor indices list to a single-row data frame
#'
#' @param results Named list as returned by bifactor_indices() with output = "list".
#' @return A single-row data frame.
#' @keywords internal
.to_dataframe <- function(results) {
  flat <- list()
  for (nm in names(results)) {
    val <- results[[nm]]
    if (is.null(val)) next
    if (length(val) == 1 && is.null(names(val))) {
      flat[[nm]] <- val
    } else if (length(val) == 1 && !is.null(names(val))) {
      flat[[nm]] <- as.numeric(val)
    } else {
      # named vector or multi-element: expand to separate columns
      for (sub_nm in names(val)) {
        col_nm <- paste0(nm, "_", sub_nm)
        flat[[col_nm]] <- as.numeric(val[[sub_nm]])
      }
    }
  }
  as.data.frame(flat, check.names = FALSE)
}
