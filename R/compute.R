############
# Core index computation
# All functions accept Lambda (items x factors matrix, general factor first column)
# and Psi (unique variance vector) as primary inputs
# No model fitting occurs here
############


############
# ECV
############

#' Compute ECV (Explained Common Variance) indices
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param gen_col Integer. Column index of the general factor. Default 1.
#' @param membership Integer vector of subscale assignments from get_subscale_membership().
#'   If NULL, computed internally.
#'
#' @return A named list with elements:
#'   \item{ECV_general}{ECV for the general factor across all items.}
#'   \item{ECV_subscale}{Named numeric vector of ECV per subscale, computed
#'     using only items belonging to that subscale (Rodriguez et al., 2016).}
#'
#' @details
#' ECV_general is the proportion of total common variance attributable to the
#' general factor: sum of squared general factor loadings divided by sum of all
#' squared loadings across all factors.
#'
#' ECV per subscale uses only items belonging to that subscale in both numerator
#' and denominator, following Rodriguez, Reise, and Haviland (2016).
#'
#' ECV is most interpretable from the standardized solution. If unstandardized
#' loadings are passed, a warning is issued.
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @export
compute_ecv <- function(Lambda, gen_col = 1, membership = NULL) {
  if (any(rowSums(Lambda^2) > 1 + 1e-6)) {
    warning(
      "Some communalities exceed 1.0, suggesting unstandardized loadings. ",
      "ECV is most interpretable from the standardized solution."
    )
  }
  if (is.null(membership)) membership <- get_subscale_membership(Lambda, gen_col)
  spec_cols <- setdiff(seq_len(ncol(Lambda)), gen_col)
  lam_g <- Lambda[, gen_col]
  lam_s <- Lambda[, spec_cols, drop = FALSE]
  total_common <- sum(lam_g^2) + sum(lam_s^2)
  ecv_g <- sum(lam_g^2) / total_common
  # per-subscale ECV using subscale items only
  n_specific <- length(spec_cols)
  ecv_s <- setNames(numeric(n_specific), colnames(lam_s))
  for (s in seq_len(n_specific)) {
    items_s <- which(membership == s)
    if (length(items_s) == 0) {
      ecv_s[s] <- NA_real_
      next
    }
    num <- sum(Lambda[items_s, spec_cols[s]]^2)
    denom <- sum(Lambda[items_s, gen_col]^2) + sum(Lambda[items_s, spec_cols[s]]^2)
    ecv_s[s] <- num / denom
  }
  list(ECV_general = ecv_g, ECV_subscale = ecv_s)
}


############
# PUC
############

#' Compute PUC (Percent Uncontaminated Correlations)
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param gen_col Integer. Column index of the general factor. Default 1.
#' @param general_only_items Integer vector of item indices (row indices of Lambda)
#'   that load only on the general factor and are not assigned to any specific factor.
#'   Default NULL assumes all items belong to exactly one specific factor, consistent
#'   with Rodriguez et al. (2016) and Dueber (2017).
#' @param tol Numeric. Loadings below this value treated as zero for counting. Default 1e-6.
#'
#' @return A single numeric value between 0 and 1.
#'
#' @details
#' PUC is the proportion of inter-item correlations that are uncontaminated by
#' multidimensionality, i.e., item pairs that do not share a specific factor.
#' Follows the implementation of Dueber (2017), which is consistent with
#' Rodriguez et al. (2016).
#'
#' When general_only_items is specified, those items are excluded from specific
#' factor pair counts but contribute to total pair counts. Their correlations
#' with all other items are uncontaminated by definition.
#'
#' @references
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based
#' tool to calculate various indices relevant to bifactor CFA models.
#' doi:10.13023/edp.tool.01
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @export
compute_puc <- function(Lambda, gen_col = 1, general_only_items = NULL, tol = 1e-6) {
  n <- nrow(Lambda)
  spec_cols <- setdiff(seq_len(ncol(Lambda)), gen_col)
  total_pairs <- n * (n - 1) / 2
  # count items per specific factor, excluding general_only_items
  contaminated_pairs <- 0
  for (s in spec_cols) {
    items_on_s <- which(abs(Lambda[, s]) > tol)
    if (!is.null(general_only_items)) {
      items_on_s <- setdiff(items_on_s, general_only_items)
    }
    n_s <- length(items_on_s)
    contaminated_pairs <- contaminated_pairs + n_s * (n_s - 1) / 2
  }
  1 - contaminated_pairs / total_pairs
}


############
# Omega indices
############

#' Compute omega reliability indices
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param Psi Numeric vector of unique variances. If NULL, computed as
#'   1 - rowSums(Lambda^2).
#' @param gen_col Integer. Column index of the general factor. Default 1.
#' @param membership Integer vector of subscale assignments. If NULL, computed internally.
#'
#' @return A named list with elements:
#'   \item{omega_h}{Omega hierarchical for the general factor (total scale).}
#'   \item{omega_s}{Named numeric vector of omega subscale per specific factor.
#'     Denominator uses only items belonging to that subscale.}
#'   \item{omega_t}{Omega total for the full scale.}
#'
#' @details
#' Omega hierarchical (omega_h) reflects the proportion of total score variance
#' attributable to the general factor. Omega subscale (omega_s) reflects the
#' proportion of subscale score variance attributable to the specific factor,
#' after controlling for the general factor. Omega total (omega_t) reflects
#' the total reliability of the composite score.
#'
#' All formulas follow Rodriguez, Reise, and Haviland (2016). The denominator
#' for omega_s uses only items belonging to that subscale, consistent with
#' the subscale-level interpretation.
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @export
compute_omega <- function(Lambda, Psi = NULL, gen_col = 1, membership = NULL) {
  if (is.null(Psi)) Psi <- 1 - rowSums(Lambda^2)
  if (is.null(membership)) membership <- get_subscale_membership(Lambda, gen_col)
  spec_cols <- setdiff(seq_len(ncol(Lambda)), gen_col)
  lam_g <- Lambda[, gen_col]
  lam_s <- Lambda[, spec_cols, drop = FALSE]
  # omega_h: general factor reliability
  num_h   <- sum(lam_g)^2
  denom_h <- sum(lam_g)^2 + sum(colSums(lam_s)^2) + sum(Psi)
  omega_h <- num_h / denom_h
  # omega_t: total reliability
  all_loadings_sum <- sum(lam_g) + sum(lam_s)
  num_t   <- all_loadings_sum^2
  denom_t <- all_loadings_sum^2 + sum(Psi)
  omega_t <- num_t / denom_t
  # omega_s: per subscale, using subscale items only
  n_specific <- length(spec_cols)
  omega_s <- setNames(numeric(n_specific), colnames(lam_s))
  for (s in seq_len(n_specific)) {
    items_s <- which(membership == s)
    if (length(items_s) == 0) {
      omega_s[s] <- NA_real_
      next
    }
    g_s   <- lam_g[items_s]
    sp_s  <- Lambda[items_s, spec_cols[s]]
    psi_s <- Psi[items_s]
    num_s   <- sum(sp_s)^2
    denom_s <- sum(g_s)^2 + sum(sp_s)^2 + sum(psi_s)
    omega_s[s] <- num_s / denom_s
  }
  list(omega_h = omega_h, omega_s = omega_s, omega_t = omega_t)
}


############
# H index (construct replicability)
############

#' Compute H (construct replicability) indices
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param gen_col Integer. Column index of the general factor. Default 1.
#'
#' @return A named list with elements:
#'   \item{H_general}{H index for the general factor.}
#'   \item{H_subscale}{Named numeric vector of H per specific factor.}
#'
#' @details
#' H (construct replicability) reflects the degree to which a latent variable
#' is well-defined by its indicators. Values above .70 are generally considered
#' adequate.
#'
#' Unique variances are computed as total communality residuals:
#' psi_i = 1 - sum_f(lambda_fi^2), using all factors in Lambda. This differs
#' from Dueber (2017), who computes psi_i = 1 - lambda_fi^2 per factor column
#' independently. The total communality approach is consistent with the
#' standardized solution interpretation and Reise (2012).
#'
#' @references
#' Reise, S. P. (2012). The rediscovery of bifactor measurement models.
#' Multivariate Behavioral Research, 47(5), 667-696.
#'
#' @export
compute_H <- function(Lambda, gen_col = 1) {
  # total communality residual for each item across all factors
  total_comm <- rowSums(Lambda^2)
  psi_total  <- 1 - total_comm
  psi_total  <- pmax(psi_total, .Machine$double.eps)
  spec_cols  <- setdiff(seq_len(ncol(Lambda)), gen_col)
  # H for general factor
  lam_g <- Lambda[, gen_col]
  snr_g <- lam_g^2 / psi_total
  H_g   <- sum(snr_g) / (1 + sum(snr_g))
  # H per specific factor
  n_specific <- length(spec_cols)
  H_s <- setNames(numeric(n_specific), colnames(Lambda)[spec_cols])
  for (s in seq_len(n_specific)) {
    lam_f <- Lambda[, spec_cols[s]]
    # only items with non-negligible loadings on this specific factor
    items_f <- which(abs(lam_f) > 1e-6)
    if (length(items_f) == 0) {
      H_s[s] <- NA_real_
      next
    }
    snr_f <- lam_f[items_f]^2 / psi_total[items_f]
    H_s[s] <- sum(snr_f) / (1 + sum(snr_f))
  }
  list(H_general = H_g, H_subscale = H_s)
}


############
# Factor Determinacy (FD)
############

#' Compute Factor Determinacy (FD) indices
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param Psi Numeric vector of unique variances. If NULL, computed as
#'   1 - rowSums(Lambda^2).
#' @param Phi Factor correlation matrix. For standard bifactor models this is the
#'   identity matrix. Default NULL constructs identity matrix of appropriate size.
#' @param gen_col Integer. Column index of the general factor. Default 1.
#'
#' @return A named numeric vector of factor determinacy values, one per factor.
#'
#' @details
#' Factor determinacy reflects the degree to which factor scores can be accurately
#' estimated from the observed indicators. The model-implied correlation matrix
#' Sigma = Lambda %*% Phi %*% t(Lambda) + diag(Psi) is constructed from the full
#' loading matrix. Follows Grice (2001) as implemented in Dueber (2017).
#'
#' @references
#' Grice, J. W. (2001). Computing and evaluating factor scores. Psychological
#' Methods, 6(4), 430-450.
#'
#' Dueber, D. M. (2017). Bifactor Indices Calculator: A Microsoft Excel-based
#' tool to calculate various indices relevant to bifactor CFA models.
#' doi:10.13023/edp.tool.01
#'
#' @export
compute_fd <- function(Lambda, Psi = NULL, Phi = NULL, gen_col = 1) {
  if (is.null(Psi)) Psi <- 1 - rowSums(Lambda^2)
  Psi <- pmax(Psi, .Machine$double.eps)
  if (is.null(Phi)) Phi <- diag(ncol(Lambda))
  Sigma  <- Lambda %*% Phi %*% t(Lambda) + diag(Psi)
  FacDet <- sqrt(diag(Phi %*% t(Lambda) %*% solve(Sigma) %*% Lambda %*% Phi))
  names(FacDet) <- colnames(Lambda)
  FacDet
}


############
# Cronbach's alpha
############

#' Compute Cronbach's alpha from factor model parameters
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#'   General factor must be the first column.
#' @param Psi Numeric vector of unique variances. If NULL, computed as
#'   1 - rowSums(Lambda^2).
#' @param Phi Factor correlation matrix. Default NULL constructs identity matrix.
#' @param gen_col Integer. Column index of the general factor. Default 1.
#' @param membership Integer vector of subscale assignments. If NULL, computed internally.
#'
#' @return A named list with elements:
#'   \item{alpha_total}{Cronbach's alpha for the full item set.}
#'   \item{alpha_subscale}{Named numeric vector of alpha per subscale.}
#'
#' @details
#' Alpha is computed from the model-implied covariance matrix rather than raw data,
#' allowing computation from any fitted model without requiring item-level responses.
#' The formula used is the standardized alpha based on the implied correlation structure:
#' alpha = (n / (n-1)) * (1 - sum(psi) / (1'Sigma 1))
#' where Sigma = Lambda Phi Lambda' + diag(Psi).
#'
#' @export
compute_alpha <- function(Lambda, Psi = NULL, Phi = NULL, gen_col = 1, membership = NULL) {
  if (is.null(Psi)) Psi <- 1 - rowSums(Lambda^2)
  Psi <- pmax(Psi, .Machine$double.eps)
  if (is.null(Phi)) Phi <- diag(ncol(Lambda))
  if (is.null(membership)) membership <- get_subscale_membership(Lambda, gen_col)
  .alpha_from_model <- function(Lam, ps, Ph) {
    n <- nrow(Lam)
    if (n < 2) return(NA_real_)
    Sig  <- Lam %*% Ph %*% t(Lam) + diag(ps)
    ones <- rep(1, n)
    # Use tr(Sigma) in numerator - correct Cronbach formula for standardized items.
    # sum(psi) would give the Green & Yang (2009) omega-adjacent formula,
    # which can exceed 1.0 under high inter-item correlations.
    (n / (n - 1)) * (1 - sum(diag(Sig)) / as.numeric(t(ones) %*% Sig %*% ones))
  }
  alpha_total <- .alpha_from_model(Lambda, Psi, Phi)
  spec_cols   <- setdiff(seq_len(ncol(Lambda)), gen_col)
  n_specific  <- length(spec_cols)
  alpha_s     <- setNames(numeric(n_specific), colnames(Lambda)[spec_cols])
  for (s in seq_len(n_specific)) {
    items_s <- which(membership == s)
    if (length(items_s) < 2) {
      alpha_s[s] <- NA_real_
      next
    }
    Lam_s <- Lambda[items_s, , drop = FALSE]
    Psi_s <- Psi[items_s]
    alpha_s[s] <- .alpha_from_model(Lam_s, Psi_s, Phi)
  }
  list(alpha_total = alpha_total, alpha_subscale = alpha_s)
}


############
# ARPB
############

#' Compute ARPB (Average Relative Parameter Bias)
#'
#' @param Lambda Numeric matrix of standardized factor loadings from the bifactor model
#'   (items x factors). General factor must be the first column.
#' @param uni_lambda Numeric vector or single-column matrix of loadings from a fitted
#'   unidimensional model. Must have the same number of items as Lambda.
#' @param gen_col Integer. Column index of the general factor in Lambda. Default 1.
#'
#' @return A named list with elements:
#'   \item{ARPB}{Mean absolute relative parameter bias across all items.}
#'   \item{SRPB}{Signed relative parameter bias per item (unidimensional - bifactor general,
#'     divided by bifactor general). Named by item.}
#'   \item{ARPB_items}{Absolute relative parameter bias per item. Named by item.}
#'
#' @details
#' ARPB quantifies the degree to which general factor loadings from a bifactor model
#' differ from the corresponding loadings in a unidimensional model. Large ARPB
#' suggests that the unidimensional model loadings are substantially biased relative
#' to the bifactor solution.
#'
#' Both signed and absolute per-item bias are returned. The signed version reveals
#' whether the unidimensional model systematically over- or underestimates loadings.
#'
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @export
compute_arpb <- function(Lambda, uni_lambda, gen_col = 1) {
  if (is.matrix(uni_lambda) || is.data.frame(uni_lambda)) {
    uni_lambda <- as.numeric(uni_lambda[, 1])
  }
  gen_loadings <- Lambda[, gen_col]
  if (length(uni_lambda) != length(gen_loadings)) {
    stop(
      "uni_lambda must have the same number of items as Lambda. ",
      "Received ", length(uni_lambda), " vs ", length(gen_loadings), " items."
    )
  }
  item_names  <- rownames(Lambda)
  signed_bias <- (uni_lambda - gen_loadings) / gen_loadings
  abs_bias    <- abs(signed_bias)
  if (!is.null(item_names)) {
    names(signed_bias) <- item_names
    names(abs_bias)    <- item_names
  }
  list(
    ARPB       = mean(abs_bias),
    SRPB       = signed_bias,
    ARPB_items = abs_bias
  )
}
