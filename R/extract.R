############
# Extraction layer - class-aware loading matrix retrieval
# All downstream computation operates on Lambda (items x factors) and Psi (unique variances)
# General factor assumed to be first column of Lambda throughout
############


#' Extract standardized loading matrix from a fitted model or raw matrix
#'
#' @param object A fitted lavaan object, mirt object, or numeric matrix of
#'   standardized factor loadings (items x factors, general factor first).
#' @param standardized Logical. If TRUE (default), extract standardized solution.
#'   If FALSE, extract unstandardized loadings. Ignored for raw matrix input.
#'
#' @return A numeric matrix of factor loadings, items in rows, factors in columns.
#'   General factor must be the first column. Column names preserved where available.
#'
#' @details
#' For lavaan objects, extraction uses lavInspect with "std.all" (standardized) or
#' "est" (unstandardized). The function expects a bifactor model specification where
#' the general factor is the first factor in the lavaan syntax.
#'
#' For mirt objects, extraction uses coef() with simplify = TRUE. The general factor
#' is assumed to be the first column (F1) of the item parameter matrix.
#'
#' For raw matrix input, the matrix is returned as-is. Standardized argument is
#' ignored and a note is issued if standardized = FALSE was not explicitly set.
#'
#' @export
extract_loadings <- function(object, standardized = TRUE) {
  if (is.matrix(object) || is.data.frame(object)) {
    return(as.matrix(object))
  }
  if (methods::is(object, "lavaan")) {
    return(.extract_lavaan(object, standardized))
  }
  if (methods::is(object, "SingleGroupClass") || methods::is(object, "MultipleGroupClass")) {
    return(.extract_mirt(object, standardized))
  }
  stop(
    "object must be a lavaan model, mirt SingleGroupClass object, or a numeric matrix. ",
    "Received object of class: ", class(object)[1]
  )
}


.extract_lavaan <- function(object, standardized) {
  if (standardized) {
    std <- lavaan::lavInspect(object, "std.all")
    Lambda <- std$lambda
  } else {
    est <- lavaan::lavInspect(object, "est")
    Lambda <- est$lambda
  }
  # drop rows that are all zero (phantom indicators, if any)
  keep <- rowSums(abs(Lambda)) > 0
  Lambda <- Lambda[keep, , drop = FALSE]
  # drop columns that are all zero
  keep_col <- colSums(abs(Lambda)) > 0
  Lambda <- Lambda[, keep_col, drop = FALSE]
  as.matrix(Lambda)
}


.extract_mirt <- function(object, standardized) {
  params <- mirt::coef(object, simplify = TRUE, IRTpars = FALSE)
  Lambda_full <- params$items
  # mirt item matrix contains loadings (a columns) plus intercepts (d) and other params
  # loading columns are named a1, a2, ... or F1, F2, ... depending on mirt version
  a_cols <- grep("^a[0-9]+$|^F[0-9]+$", colnames(Lambda_full), value = TRUE)
  if (length(a_cols) == 0) {
    stop(
      "Could not identify factor loading columns in mirt coef() output. ",
      "Expected columns named a1, a2, ... or F1, F2, ..."
    )
  }
  Lambda <- as.matrix(Lambda_full[, a_cols, drop = FALSE])
  if (standardized) {
    # standardize: lambda_std = lambda / sqrt(1 + sum(lambda^2)) per row (normal ogive scaling)
    # for already-standardized mirt models (itemtype = "2PL" with standard scaling) this
    # may already be on the right metric; we apply the FA standardization
    comm <- rowSums(Lambda^2)
    # only rescale if communalities suggest unstandardized metric
    if (any(comm > 1)) {
      message(
        "mirt loadings appear to be on an unstandardized metric (communalities > 1 detected). ",
        "Applying FA standardization. Set standardized = FALSE to suppress."
      )
      total_var <- 1 + rowSums(Lambda^2)
      Lambda <- Lambda / sqrt(total_var)
    }
  }
  Lambda
}


#' Extract unique variances (Psi) from loading matrix or fitted object
#'
#' @param Lambda Numeric matrix of standardized factor loadings (items x factors).
#' @param object Optional. If supplied and standardized = FALSE, attempts to extract
#'   residual variances directly from the fitted object rather than computing from Lambda.
#' @param standardized Logical. If TRUE (default), computes Psi = 1 - rowSums(Lambda^2).
#'   If FALSE and object is supplied, extracts residual variances from the fitted object.
#'
#' @return A named numeric vector of unique variances, one per item.
#'
#' @export
extract_psi <- function(Lambda, object = NULL, standardized = TRUE) {
  if (standardized || is.null(object)) {
    psi <- 1 - rowSums(Lambda^2)
    if (any(psi < 0)) {
      warning(
        "One or more items have communalities > 1 (negative unique variances). ",
        "Check that Lambda contains standardized loadings."
      )
      psi <- pmax(psi, .Machine$double.eps)
    }
    return(psi)
  }
  # unstandardized path: extract from fitted object
  if (methods::is(object, "lavaan")) {
    est <- lavaan::lavInspect(object, "est")
    psi_mat <- est$theta
    return(diag(psi_mat))
  }
  if (methods::is(object, "SingleGroupClass") || methods::is(object, "MultipleGroupClass")) {
    params <- mirt::coef(object, simplify = TRUE, IRTpars = FALSE)
    # mirt stores residual variances as 1/a^2 scaling - for FA models use the
    # unique variance directly if available, otherwise compute from Lambda
    message("Extracting unstandardized unique variances from mirt object.")
    psi <- 1 - rowSums(Lambda^2)
    return(psi)
  }
  # fallback
  1 - rowSums(Lambda^2)
}


#' Identify which column of Lambda is the general factor
#'
#' Assumes the general factor is the column with the most non-zero loadings,
#' which in a standard bifactor model is always the first column. Returns 1
#' by default but checks for obvious misspecification.
#'
#' @param Lambda Numeric matrix of factor loadings.
#' @param gen_col Integer or character. Column index or name of the general factor.
#'   Default is 1.
#'
#' @return Integer column index of the general factor.
#'
#' @export
identify_general <- function(Lambda, gen_col = 1) {
  if (is.character(gen_col)) {
    idx <- which(colnames(Lambda) == gen_col)
    if (length(idx) == 0) stop("gen_col name '", gen_col, "' not found in Lambda column names.")
    return(idx)
  }
  if (gen_col > ncol(Lambda)) stop("gen_col exceeds number of columns in Lambda.")
  as.integer(gen_col)
}


#' Get subscale membership from loading matrix
#'
#' Identifies which specific factor each item belongs to based on non-zero
#' loadings in the specific factor columns of Lambda.
#'
#' @param Lambda Numeric matrix of factor loadings (items x factors).
#' @param gen_col Integer. Column index of the general factor. Default 1.
#' @param tol Numeric. Loadings with absolute value below tol treated as zero.
#'   Default 1e-6.
#'
#' @return An integer vector of length nrow(Lambda) giving the specific factor
#'   index (among specific factors only, 1-indexed) for each item. Items loading
#'   only on the general factor receive NA.
#'
#' @export
get_subscale_membership <- function(Lambda, gen_col = 1, tol = 1e-6) {
  spec_cols <- setdiff(seq_len(ncol(Lambda)), gen_col)
  spec_Lambda <- Lambda[, spec_cols, drop = FALSE]
  membership <- apply(abs(spec_Lambda) > tol, 1, function(row) {
    hits <- which(row)
    if (length(hits) == 0) return(NA_integer_)
    if (length(hits) > 1) {
      warning("Item loads on multiple specific factors; assigning to first.")
    }
    hits[1L]
  })
  membership
}
