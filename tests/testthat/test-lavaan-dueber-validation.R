test_that("lavaan extraction agrees with Dueber where comparable", {
  
  skip_if_not_installed("lavaan")
  skip_if_not_installed("BifactorIndicesCalculator")
  
  bifactor_model <- "
    g  =~ x1 + x2 + x3 + x4 + x5 + x6
    s1 =~ x1 + x2 + x3
    s2 =~ x4 + x5 + x6

    g ~~ 0*s1
    g ~~ 0*s2
    s1 ~~ 0*s2
  "
  
  fit <- lavaan::cfa(
    bifactor_model,
    data = lavaan::HolzingerSwineford1939[, paste0("x", 1:6)],
    std.lv = TRUE,
    se = "none"
  )
  
  ours <- bifactor_indices(fit)
  theirs <- BifactorIndicesCalculator::bifactorIndices(fit)
  
  expect_equal(unname(ours$ECV_general), unname(theirs$ModelLevelIndices["ECV.g"]), tolerance = 1e-6)
  expect_equal(unname(ours$PUC), unname(theirs$ModelLevelIndices["PUC"]), tolerance = 1e-6)
  expect_equal(unname(ours$omega_h), unname(theirs$ModelLevelIndices["OmegaH.g"]), tolerance = 1e-6)
  expect_equal(unname(ours$omega_s), unname(theirs$FactorLevelIndices[c("s1", "s2"), "OmegaH"]), tolerance = 1e-6)
  expect_equal(unname(ours$FD), unname(theirs$FactorLevelIndices[c("g", "s1", "s2"), "FD"]), tolerance = 1e-6)
})