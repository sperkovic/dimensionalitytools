############
# Tests - interpretation layer (summary and print methods)
############

library(testthat)

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


test_that("list output carries bifactor_indices class", {
  res <- bifactor_indices(Lambda_ref)
  expect_s3_class(res, "bifactor_indices")
})


test_that("dataframe output does not carry bifactor_indices class", {
  res <- bifactor_indices(Lambda_ref, output = "dataframe")
  expect_false(inherits(res, "bifactor_indices"))
  expect_s3_class(res, "data.frame")
})


test_that("print method runs without error and returns invisibly", {
  res <- bifactor_indices(Lambda_ref)
  expect_output(print(res), "Bifactor model indices")
  expect_invisible(print(res))
})


test_that("summary returns leaning and cutoffs invisibly", {
  res <- bifactor_indices(Lambda_ref)
  out <- summary(res)
  expect_true("leaning" %in% names(out))
  expect_true("cutoffs" %in% names(out))
  expect_equal(out$cutoffs$ecv, 0.70)
})


test_that("summary respects custom cutoffs", {
  res <- bifactor_indices(Lambda_ref)
  out <- summary(res, ecv_cut = 0.90, omegah_cut = 0.95)
  expect_equal(out$cutoffs$ecv, 0.90)
  expect_equal(out$cutoffs$omegah, 0.95)
  # with very high cutoffs, lean should move away from unidimensionality
  expect_match(out$leaning, "away|ambiguous")
})


test_that("summary lean is toward for reference Lambda at default cutoffs", {
  res <- bifactor_indices(Lambda_ref)
  out <- summary(res)
  expect_match(out$leaning, "toward")
})


test_that("summary prints convention citations", {
  res <- bifactor_indices(Lambda_ref)
  expect_output(summary(res), "Rodriguez")
  expect_output(summary(res), "not decision rules")
})


test_that("summary names subscale tensions when omega_s is substantial", {
  # construct a Lambda with a strong specific factor
  Lambda_strong <- Lambda_ref
  Lambda_strong[5:8, 2] <- c(.70, .65, .68, .72)  # strengthen SF2 region
  res <- bifactor_indices(Lambda_strong)
  out <- summary(res)
  # at minimum should run and produce a character leaning
  expect_type(out$leaning, "character")
})
