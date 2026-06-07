test_that("mirt bfactor object can be processed", {
  
  skip_if_not_installed("mirt")
  
  set.seed(123)
  
  a <- matrix(
    c(
      1.0, 0.6, 0.0,
      1.0, 0.6, 0.0,
      1.0, 0.6, 0.0,
      1.0, 0.0, 0.6,
      1.0, 0.0, 0.6,
      1.0, 0.0, 0.6
    ),
    ncol = 3,
    byrow = TRUE
  )
  
  dat <- mirt::simdata(
    a = a,
    d = rep(0, 6),
    itemtype = rep("2PL", 6),
    N = 500
  )
  
  mod <- mirt::bfactor(
    dat,
    model = c(1, 1, 1, 2, 2, 2),
    itemtype = "2PL",
    technical = list(NCYCLES = 1000),
    verbose = FALSE
  )
  
  expect_message(
    res <- bifactor_indices(mod),
    "mirt loadings appear to be on an unstandardized metric"
  )
  
  expect_s3_class(res, "bifactor_indices")
  
  expect_true(is.numeric(res$ECV_general))
  expect_true(is.numeric(res$PUC))
  expect_true(is.numeric(res$omega_h))
  expect_true(is.numeric(res$omega_t))
  
  expect_true(res$ECV_general >= 0 && res$ECV_general <= 1)
  expect_true(res$PUC >= 0 && res$PUC <= 1)
  expect_true(res$omega_h >= 0 && res$omega_h <= 1)
  expect_true(res$omega_t >= 0 && res$omega_t <= 1)
  
  expect_equal(length(res$FD), 3)
  expect_equal(length(res$omega_s), 2)
})

