############
# Tests - math layer verification
# Reference values cross-checked against BifactorIndicesCalculator (Dueber, 2017)
# using the standard 12-item 3-specific-factor example Lambda
############

library(testthat)

# Dueber example Lambda - used as reference throughout
Lambda_ref <- matrix(c(
  .82, .10,  0,   0,
  .77, .35,  0,   0,
  .79, .32,  0,   0,
  .66, .39,  0,   0,
  .51,  0,  .71,  0,
  .56,  0,  .43,  0,
  .68,  0,  .13,  0,
  .60,  0,  .50,  0,
  .83,  0,   0,  .47,
  .60,  0,   0,  .27,
  .78,  0,   0,  .28,
  .55,  0,   0,  .75
), ncol = 4, byrow = TRUE)
colnames(Lambda_ref) <- c("G", "SF1", "SF2", "SF3")
rownames(Lambda_ref) <- paste0("Item", 1:12)


test_that("ECV_general matches Dueber reference", {
  ecv <- compute_ecv(Lambda_ref)
  # total common variance
  total <- sum(Lambda_ref^2)
  expected_g <- sum(Lambda_ref[, 1]^2) / total
  expect_equal(ecv$ECV_general, expected_g, tolerance = 1e-6)
  # rough range check consistent with Rodriguez et al. examples
  expect_gt(ecv$ECV_general, 0.5)
  expect_lt(ecv$ECV_general, 1.0)
})


test_that("ECV_subscale uses subscale items only", {
  ecv <- compute_ecv(Lambda_ref)
  # SF1: items 1-4
  items_sf1 <- 1:4
  num  <- sum(Lambda_ref[items_sf1, "SF1"]^2)
  denom <- sum(Lambda_ref[items_sf1, "G"]^2) + sum(Lambda_ref[items_sf1, "SF1"]^2)
  expect_equal(ecv$ECV_subscale["SF1"], num / denom, tolerance = 1e-6)
})


test_that("PUC matches Dueber formula", {
  puc <- compute_puc(Lambda_ref)
  # manual: 4 items SF1, 4 items SF2, 4 items SF3
  # contaminated = 4*3/2 * 3 = 18
  # total pairs = 12*11/2 = 66
  expected <- 1 - 18/66
  expect_equal(puc, expected, tolerance = 1e-6)
})


test_that("PUC with general_only_items reduces contaminated pairs", {
  # if item 1 is general only, SF1 drops to 3 items
  puc_standard <- compute_puc(Lambda_ref)
  puc_genonly  <- compute_puc(Lambda_ref, general_only_items = 1L)
  # contaminated: SF1 now 3 items = 3, SF2 = 6, SF3 = 6 -> 15 contaminated
  expected <- 1 - 15/66
  expect_equal(puc_genonly, expected, tolerance = 1e-6)
  expect_gt(puc_genonly, puc_standard)
})


test_that("omega_h is between 0 and 1 and less than omega_t", {
  omega <- compute_omega(Lambda_ref)
  expect_gt(omega$omega_h, 0)
  expect_lt(omega$omega_h, 1)
  expect_gt(omega$omega_t, omega$omega_h)
  expect_lt(omega$omega_t, 1)
})


test_that("omega_s values are between 0 and 1", {
  omega <- compute_omega(Lambda_ref)
  for (s in omega$omega_s) {
    expect_gt(s, 0)
    expect_lt(s, 1)
  }
})


test_that("H_general uses total communality psi", {
  H <- compute_H(Lambda_ref)
  # manually compute for general factor
  psi_total <- 1 - rowSums(Lambda_ref^2)
  lam_g <- Lambda_ref[, 1]
  snr <- lam_g^2 / psi_total
  expected_H_g <- sum(snr) / (1 + sum(snr))
  expect_equal(H$H_general, expected_H_g, tolerance = 1e-8)
})


test_that("H differs from Dueber per-factor psi approach", {
  H_total_comm <- compute_H(Lambda_ref)
  # Dueber approach: psi per column = 1 - lambda_col^2
  H_dueber_g <- 1 / (1 + 1 / sum(Lambda_ref[,1]^2 / (1 - Lambda_ref[,1]^2)))
  # should differ since total communality psi != per-factor psi for bifactor models
  expect_false(isTRUE(all.equal(H_total_comm$H_general, H_dueber_g, tolerance = 1e-4)))
})


test_that("FD values are between 0 and 1", {
  fd <- compute_fd(Lambda_ref)
  expect_true(all(fd >= 0 & fd <= 1))
  expect_length(fd, ncol(Lambda_ref))
})


test_that("FD matches Dueber formula", {
  Psi   <- 1 - rowSums(Lambda_ref^2)
  Phi   <- diag(4)
  Sigma <- Lambda_ref %*% Phi %*% t(Lambda_ref) + diag(Psi)
  expected <- sqrt(diag(Phi %*% t(Lambda_ref) %*% solve(Sigma) %*% Lambda_ref %*% Phi))
  fd <- compute_fd(Lambda_ref)
  expect_equal(fd, expected, tolerance = 1e-8)
})


test_that("alpha_total is between 0 and 1", {
  alpha <- compute_alpha(Lambda_ref)
  expect_gt(alpha$alpha_total, 0)
  expect_lt(alpha$alpha_total, 1)
})


test_that("alpha_subscale values are between 0 and 1", {
  alpha <- compute_alpha(Lambda_ref)
  for (a in alpha$alpha_subscale) {
    if (!is.na(a)) {
      expect_gt(a, 0)
      expect_lt(a, 1)
    }
  }
})


test_that("ARPB returns correct structure", {
  uni_lambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
  arpb <- compute_arpb(Lambda_ref, uni_lambda)
  expect_named(arpb, c("ARPB", "SRPB", "ARPB_items"))
  expect_length(arpb$ARPB, 1)
  expect_length(arpb$SRPB, 12)
  expect_length(arpb$ARPB_items, 12)
  expect_true(all(arpb$ARPB_items >= 0))
})


test_that("ARPB mean equals mean of absolute item-level bias", {
  uni_lambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
  arpb <- compute_arpb(Lambda_ref, uni_lambda)
  expect_equal(arpb$ARPB, mean(arpb$ARPB_items), tolerance = 1e-8)
})


test_that("bifactor_indices list output contains all expected elements", {
  res <- bifactor_indices(Lambda_ref)
  expected_names <- c("ECV_general", "ECV_subscale", "PUC", "omega_h",
                      "omega_s", "omega_t", "H_general", "H_subscale",
                      "FD", "alpha_total", "alpha_subscale")
  expect_true(all(expected_names %in% names(res)))
})


test_that("bifactor_indices dataframe output is a single-row dataframe", {
  res <- bifactor_indices(Lambda_ref, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  # subscale columns should be expanded
  expect_true(any(grepl("ECV_subscale_", colnames(res))))
  expect_true(any(grepl("omega_s_", colnames(res))))
})


test_that("bifactor_indices with uni_object returns ARPB", {
  uni_lambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
  res <- bifactor_indices(Lambda_ref, uni_object = uni_lambda)
  expect_true("ARPB" %in% names(res))
  expect_true("SRPB" %in% names(res))
})


test_that("bifactor_indices without uni_object does not return ARPB", {
  res <- bifactor_indices(Lambda_ref)
  expect_false("ARPB" %in% names(res))
})


test_that("dataframe output binds cleanly across replications", {
  # simulate SimDesign-like row-binding
  res1 <- bifactor_indices(Lambda_ref, output = "dataframe")
  res2 <- bifactor_indices(Lambda_ref, output = "dataframe")
  combined <- rbind(res1, res2)
  expect_equal(nrow(combined), 2)
  expect_equal(ncol(combined), ncol(res1))
})
