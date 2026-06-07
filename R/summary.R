############
# Interpretation layer
# summary() reorganizes bifactor_indices output around the essential
# unidimensionality question and narrates the Reise-style conditional logic.
# No recomputation occurs here - this is a presentation layer only.
############


#' Print method for bifactor_indices objects
#'
#' @param x A bifactor_indices object returned by bifactor_indices() with output = "list".
#' @param digits Integer. Number of digits to display. Default 3.
#' @param ... Unused, for S3 consistency.
#'
#' @return Invisibly returns x. Called for its printing side effect.
#'
#' @export
print.bifactor_indices <- function(x, digits = 3, ...) {
  rnd <- function(v) round(v, digits)
  cat("Bifactor model indices\n")
  cat("----------------------\n")
  cat("ECV (general):  ", rnd(x$ECV_general), "\n")
  cat("PUC:            ", rnd(x$PUC), "\n")
  cat("omega_h:        ", rnd(x$omega_h), "\n")
  cat("omega_t:        ", rnd(x$omega_t), "\n")
  cat("H (general):    ", rnd(x$H_general), "\n")
  cat("\nSubscale indices\n")
  cat("ECV (subscale): ", paste(names(x$ECV_subscale), rnd(x$ECV_subscale), sep = "=", collapse = "  "), "\n")
  cat("omega_s:        ", paste(names(x$omega_s), rnd(x$omega_s), sep = "=", collapse = "  "), "\n")
  cat("H (subscale):   ", paste(names(x$H_subscale), rnd(x$H_subscale), sep = "=", collapse = "  "), "\n")
  cat("\nFactor determinacy: ", paste(names(x$FD), rnd(x$FD), sep = "=", collapse = "  "), "\n")
  cat("alpha (total):  ", rnd(x$alpha_total), "\n")
  if (!is.null(x$ARPB)) cat("ARPB:           ", rnd(x$ARPB), "\n")
  cat("\nUse summary() for an essential-unidimensionality interpretation.\n")
  invisible(x)
}


#' Summarize bifactor indices around essential unidimensionality
#'
#' Reorganizes a bifactor_indices object into the indices most relevant to
#' judging essential unidimensionality, displays each alongside a commonly cited
#' convention, narrates the conditional interpretive logic, and offers a tentative
#' lean. Conventions are shown as reference points, not decision rules.
#'
#' @param object A bifactor_indices object from bifactor_indices() with output = "list".
#' @param ecv_cut Numeric. Commonly cited ECV reference value. Default 0.70
#'   (Rodriguez, Reise, & Haviland, 2016; Reise et al., 2013).
#' @param omegah_cut Numeric. Commonly cited omega hierarchical reference value.
#'   Default 0.80 (Reise, 2012; Rodriguez et al., 2016).
#' @param puc_cut Numeric. PUC value above which general factor bias tends to be
#'   minimal and ECV/omega_h carry less weight. Default 0.80 (Reise et al., 2013).
#' @param omegas_cut Numeric. omega subscale value above which a specific factor
#'   is considered to carry substantial reliable variance beyond the general factor.
#'   Default 0.30 (Rodriguez et al., 2016).
#' @param H_cut Numeric. Construct replicability reference value. Default 0.70
#'   (Hancock & Mueller, 2001; Rodriguez et al., 2016).
#' @param digits Integer. Display digits. Default 3.
#' @param ... Unused, for S3 consistency.
#'
#' @return Invisibly returns a list with the organized indices, the cutoffs used,
#'   and the verdict text. Called primarily for its printed interpretation.
#'
#' @details
#' The narration follows the conditional logic described by Reise et al. (2013)
#' and Rodriguez et al. (2016): PUC moderates how much weight ECV and omega_h
#' should carry. When PUC is high, general factor parameter bias from fitting a
#' unidimensional model tends to be minimal even at moderate ECV, so omega_h
#' carries most of the interpretive weight. When PUC is lower, ECV and omega_h
#' become load-bearing for any claim of essential unidimensionality.
#'
#' The function does not return a categorical verdict. It states which direction
#' the evidence leans and names tensions when indices disagree, leaving the
#' judgment to the analyst.
#'
#' @references
#' Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct reliability
#' within latent variable systems. In Structural Equation Modeling: Present and
#' Future (pp. 195-216).
#'
#' Reise, S. P. (2012). The rediscovery of bifactor measurement models.
#' Multivariate Behavioral Research, 47(5), 667-696.
#'
#' Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
#' Multidimensionality and structural coefficient bias in structural equation
#' modeling: A bifactor perspective. Educational and Psychological Measurement,
#' 73(1), 5-26.
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices. Psychological
#' Methods, 21(2), 137-150.
#'
#' @export
summary.bifactor_indices <- function(object,
                                     ecv_cut    = 0.70,
                                     omegah_cut = 0.80,
                                     puc_cut    = 0.80,
                                     omegas_cut = 0.30,
                                     H_cut      = 0.70,
                                     digits     = 3,
                                     ...) {

  rnd <- function(v) round(v, digits)
  ecv    <- object$ECV_general
  puc    <- object$PUC
  omegah <- object$omega_h
  Hgen   <- object$H_general
  omegas <- object$omega_s

  ############
  # Primary indices table
  ############
  cat("Essential unidimensionality summary\n")
  cat("===================================\n\n")
  cat("Primary indices (general factor strength)\n")
  cat(sprintf("  %-14s %7s   convention: %s\n", "ECV (general)", rnd(ecv),
              paste0("> ", ecv_cut, " (Rodriguez et al., 2016)")))
  cat(sprintf("  %-14s %7s   convention: %s\n", "PUC", rnd(puc),
              paste0("> ", puc_cut, " (Reise et al., 2013)")))
  cat(sprintf("  %-14s %7s   convention: %s\n", "omega_h", rnd(omegah),
              paste0("> ", omegah_cut, " (Reise, 2012)")))
  cat(sprintf("  %-14s %7s   convention: %s\n", "H (general)", rnd(Hgen),
              paste0("> ", H_cut, " (Hancock & Mueller, 2001)")))
  cat("\n  Conventions are commonly cited reference points, not decision rules.\n\n")

  ############
  # Subscale strength flag
  ############
  strong_sub <- names(omegas)[!is.na(omegas) & omegas >= omegas_cut]
  cat("Specific factor reliable variance (omega_s)\n")
  for (nm in names(omegas)) {
    flag <- if (!is.na(omegas[nm]) && omegas[nm] >= omegas_cut) " *" else ""
    cat(sprintf("  %-8s %7s%s\n", nm, rnd(omegas[nm]), flag))
  }
  cat(sprintf("  (* omega_s >= %s; suggests reliable specific variance beyond the general factor.\n     Rodriguez et al., 2016. Note: * is a threshold flag, not a significance test or p-value.)\n\n",
              omegas_cut))

  ############
  # Conditional narration
  ############
  cat("Interpretation\n")
  cat("--------------\n")
  puc_high <- puc >= puc_cut
  if (puc_high) {
    cat(sprintf(
      "PUC = %s is at or above %s. General factor parameter bias from a\nunidimensional treatment tends to be minimal in this range, so omega_h\ncarries most of the interpretive weight (Reise et al., 2013).\n\n",
      rnd(puc), puc_cut))
  } else {
    cat(sprintf(
      "PUC = %s is below %s. ECV and omega_h are load-bearing here: a claim of\nessential unidimensionality rests on both reaching convention, since the\nproportion of contaminated correlations is non-trivial (Reise et al., 2013).\n\n",
      rnd(puc), puc_cut))
  }

  ecv_ok    <- ecv    >= ecv_cut
  omegah_ok <- omegah >= omegah_cut

  ############
  # Tentative lean
  ############
  # determine direction and surface tensions
  tensions <- character(0)
  if (length(strong_sub) > 0) {
    tensions <- c(tensions, sprintf(
      "subscale(s) %s show omega_s >= %s, indicating reliable specific variance that a purely unidimensional treatment would discard",
      paste(strong_sub, collapse = ", "), omegas_cut))
  }
  if (omegah_ok && !ecv_ok) {
    tensions <- c(tensions,
      "omega_h reaches convention but ECV does not, so common variance is reliable yet not strongly general-factor dominated")
  }
  if (ecv_ok && !omegah_ok) {
    tensions <- c(tensions,
      "ECV reaches convention but omega_h does not, so the general factor dominates common variance but total-score reliability for that factor is modest")
  }

  if (puc_high) {
    leaning <- if (omegah_ok) "toward" else "away from"
  } else {
    leaning <- if (ecv_ok && omegah_ok) "toward" else if (!ecv_ok && !omegah_ok) "away from" else "ambiguously on"
  }
  if (length(tensions) > 0 && leaning == "toward") leaning <- "tentatively toward (with caveats)"

  cat(sprintf("Tentative lean: evidence leans %s treating this scale as\nessentially unidimensional.\n", leaning))
  if (length(tensions) > 0) {
    cat("\nNoted tensions:\n")
    for (t in tensions) cat("  - ", t, "\n", sep = "")
  }
  cat("\nThis is a reading of the index constellation, not a categorical verdict.\nThe conventions above are choices; adjust the cutoff arguments to reflect\nyour own standards.\n")

  invisible(list(
    primary = list(ECV_general = ecv, PUC = puc, omega_h = omegah, H_general = Hgen),
    omega_s = omegas,
    cutoffs = list(ecv = ecv_cut, omegah = omegah_cut, puc = puc_cut,
                   omegas = omegas_cut, H = H_cut),
    leaning = leaning,
    tensions = tensions
  ))
}
