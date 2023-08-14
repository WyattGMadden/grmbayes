test_that("10 iter, GP, exponential covariance", {
              ctm_fit <- grm(Y = ctm_pm25$pm25,
                             X = ctm_pm25$ctm,
                             L = ctm_pm25[, c("elevation", "forestcover",
                                              "hwy_length", "lim_hwy_length", 
                                              "local_rd_length", "point_emi_any")],
                             M = ctm_pm25[, c("tmp", "wind")],
                             n.iter = 10,
                             burn = 2,
                             thin = 4,
                             nngp = F,
                             covariance = "exponential",
                             coords = ctm_pm25[, c("x", "y")],
                             space.id = ctm_pm25$space_id,
                             time.id = ctm_pm25$time_id,
                             spacetime.id = ctm_pm25$spacetime_id,
                             verbose.iter = 10)
              expect_equal(sum(is.na(ctm_fit$alpha.time)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.time)), 0)
              expect_equal(sum(is.na(ctm_fit$alpha.space)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.space)), 0)
              expect_equal(sum(is.na(ctm_fit$others)), 0)
})

test_that("10 iter, NNGP, exponential covariance", {
              ctm_fit <- grm(Y = ctm_pm25$pm25,
                             X = ctm_pm25$ctm,
                             L = ctm_pm25[, c("elevation", "forestcover",
                                              "hwy_length", "lim_hwy_length", 
                                              "local_rd_length", "point_emi_any")],
                             M = ctm_pm25[, c("tmp", "wind")],
                             n.iter = 10,
                             burn = 2,
                             thin = 4,
                             nngp = T,
                             covariance = "exponential",
                             coords = ctm_pm25[, c("x", "y")],
                             space.id = ctm_pm25$space_id,
                             time.id = ctm_pm25$time_id,
                             spacetime.id = ctm_pm25$spacetime_id,
                             verbose.iter = 10)
              expect_equal(sum(is.na(ctm_fit$alpha.time)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.time)), 0)
              expect_equal(sum(is.na(ctm_fit$alpha.space)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.space)), 0)
              expect_equal(sum(is.na(ctm_fit$others)), 0)
})
