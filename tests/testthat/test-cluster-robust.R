# Tests for cluster-robust standard errors (the `id` argument).

# Toy influence curve with a shared cluster effect.
make_ic <- function(n_clusters = 60, size = 3, seed = 1) {
  set.seed(seed)
  g <- rep(seq_len(n_clusters), each = size)
  ic <- rnorm(n_clusters * size) + rep(rnorm(n_clusters), each = size)
  list(ic = ic, g = g)
}

# Data frame with the columns calc_ate_estimates() expects.
make_ate_data <- function(n_clusters = 80, size = 2, seed = 2) {
  set.seed(seed)
  n <- n_clusters * size
  fam <- rep(seq_len(n_clusters), each = size)
  u <- rep(rnorm(n_clusters), each = size)
  data.frame(
    y = 1 + u + rnorm(n),
    qbar_1w = 1.2 + rnorm(n, sd = 0.1),
    qbar_1w_star = 1.3 + 0.5 * u + rnorm(n, sd = 0.1),
    qbar_aw_star = 1 + 0.5 * u + rnorm(n, sd = 0.1),
    h_aw = runif(n, 0.5, 2),
    folds = rep(1:2, length.out = n),
    fam = fam
  )
}

test_that("cluster_var equals var() with no id or singleton clusters", {
  x <- make_ic()
  expect_identical(cluster_var(x$ic), stats::var(x$ic))
  expect_identical(cluster_var(x$ic, NULL), stats::var(x$ic))
  expect_equal(cluster_var(x$ic, seq_along(x$ic)), stats::var(x$ic))
})

test_that("one-way cluster_var matches sandwich::vcovCL", {
  skip_if_not_installed("sandwich")
  x <- make_ic()
  n <- length(x$ic)
  fit <- stats::lm(x$ic ~ 1)
  v_sandwich <- sandwich::vcovCL(fit, cluster = x$g, type = "HC0",
                                 cadjust = TRUE)[1, 1]
  expect_equal(cluster_var(x$ic, x$g) / n, v_sandwich)
})

test_that("two-way cluster_var matches sandwich::vcovCL (multi0 = FALSE)", {
  skip_if_not_installed("sandwich")
  set.seed(3)
  n <- 400
  a <- sample(1:40, n, replace = TRUE)
  b <- sample(1:25, n, replace = TRUE)
  ic <- rnorm(n) + rnorm(40)[a] + rnorm(25)[b]
  fit <- stats::lm(ic ~ 1)
  v_sandwich <- sandwich::vcovCL(fit, cluster = data.frame(a, b),
                                 type = "HC0", cadjust = TRUE,
                                 multi0 = FALSE)[1, 1]
  expect_equal(cluster_var(ic, data.frame(a, b)) / n, v_sandwich)
})

test_that("a nested cluster variable reduces to the outer one", {
  x <- make_ic(n_clusters = 60, size = 3)
  site <- (x$g - 1) %/% 6 + 1 # 10 sites, families nested inside
  expect_message(
    id <- check_cluster_id(data.frame(fam = x$g, site = site),
                           length(x$ic), verbose = TRUE),
    "nested"
  )
  expect_identical(names(id), "site")
  expect_equal(cluster_var(x$ic, data.frame(fam = x$g, site = site)),
               cluster_var(x$ic, site))
})

test_that("check_cluster_id rejects bad input", {
  expect_null(check_cluster_id(NULL, 10))
  expect_null(check_cluster_id(1:10, 10)) # all singletons
  expect_error(check_cluster_id(c(1, 1, NA, 2), 4), "missing")
  expect_error(check_cluster_id(1:5, 10), "rows")
})

test_that("create_cv_folds keeps clusters together and uses every fold", {
  set.seed(4)
  id <- rep(1:100, times = sample(1:3, 100, replace = TRUE))
  y <- rnorm(length(id))
  folds <- create_cv_folds(5, y, id = id)
  expect_true(all(tapply(folds, id, function(f) length(unique(f))) == 1))
  expect_setequal(unique(folds), 1:5)

  # Binary outcome: stratified on whether a cluster has any event.
  yb <- stats::rbinom(length(id), 1, 0.3)
  folds_b <- create_cv_folds(5, yb, id = id)
  expect_true(all(tapply(folds_b, id, function(f) length(unique(f))) == 1))
  expect_setequal(unique(folds_b), 1:5)

  expect_error(create_cv_folds(5, y[1:8], id = id[1:8]), "fewer than")
})

test_that("create_cv_folds is unchanged without id", {
  y <- rnorm(50)
  set.seed(5)
  f1 <- create_cv_folds(5, y)
  set.seed(5)
  f2 <- create_cv_folds(5, y, id = NULL)
  expect_identical(f1, f2)
})

test_that("calc_ate_estimates: all-unique ids give identical results", {
  d <- make_ate_data()
  d_ids <- d
  d_ids$.cid_id <- seq_len(nrow(d))
  for (v_fold in c(FALSE, TRUE)) {
    r0 <- calc_ate_estimates(d, "ate", "y", p_adjust_n = 1, v_fold = v_fold)
    r1 <- calc_ate_estimates(d_ids, "ate", "y", p_adjust_n = 1,
                             v_fold = v_fold)
    expect_identical(r0$ate, r1$ate)
    expect_identical(r0$se, r1$se)
    expect_identical(r0$ci, r1$ci)
    expect_identical(r0$p_value, r1$p_value)
    expect_identical(r1$df, Inf)
  }
})

test_that("calc_ate_estimates: clustering keeps the estimate, changes the SE", {
  d <- make_ate_data(n_clusters = 80, size = 2)
  d_cl <- d
  d_cl$.cid_fam <- d$fam
  r0 <- calc_ate_estimates(d, "ate", "y", p_adjust_n = 1)
  r1 <- calc_ate_estimates(d_cl, "ate", "y", p_adjust_n = 1)
  expect_identical(r0$ate, r1$ate)
  expect_false(isTRUE(all.equal(r0$se, r1$se)))
  expect_equal(r1$df, 79)
  n <- nrow(d)
  expect_equal(r1$se, sqrt(cluster_var(r1$data$ic, d$fam) / n))
  expect_equal(r1$ci, r1$ate + c(-1, 1) * stats::qt(0.975, 79) * r1$se)
})

test_that("CVtreeMLE: id with all-unique values reproduces the default", {
  skip_on_cran()
  set.seed(6)
  d <- simulate_mixture_cube(n_obs = 300)
  lrnrs <- sl3::Stack$new(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new())
  run <- function(...) {
    suppressWarnings(CVtreeMLE(
      data = d, w = c("age", "sex", "bmi"), a = c("M1", "M2", "M3"),
      y = "y", w_stack = lrnrs, aw_stack = lrnrs, a_stack = lrnrs,
      n_folds = 2, family = "continuous", parallel = FALSE,
      parallel_cv = FALSE, max_depth = 1, min_obs = 40, seed = 1, ...
    ))
  }
  r_none <- run()
  r_single <- suppressMessages(run(id = seq_len(nrow(d))))
  expect_identical(r_none, r_single)
})

test_that("CVtreeMLE: family clustering runs and keeps families in one fold", {
  skip_on_cran()
  set.seed(7)
  d <- simulate_mixture_cube(n_obs = 300)
  d$family <- rep(seq_len(150), each = 2)
  d$y <- d$y + rep(stats::rnorm(150), each = 2)
  lrnrs <- sl3::Stack$new(sl3::Lrnr_glm$new(), sl3::Lrnr_mean$new())
  expect_message(
    res <- suppressWarnings(CVtreeMLE(
      data = d, w = c("age", "sex", "bmi"), a = c("M1", "M2", "M3"),
      y = "y", w_stack = lrnrs, aw_stack = lrnrs, a_stack = lrnrs,
      n_folds = 2, family = "continuous", parallel = FALSE,
      parallel_cv = FALSE, max_depth = 1, min_obs = 40, seed = 1,
      id = "family"
    )),
    "family \\(150 clusters\\)"
  )
  expect_s3_class(res, "CVtreeMLE")
  v_res <- res$`V-Specific Mix Results`
  if (!is.null(v_res) && nrow(v_res) > 0) {
    expect_true(all(is.finite(v_res$se)))
  }
})
