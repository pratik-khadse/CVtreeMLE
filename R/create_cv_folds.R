#' @title Stratified CV to insure balance (by one grouping variable, Y)
#' @description Creates a dummy variable that partitions the data into v
#' equal sized groups for v-fold CV.
#' @param v number of folds
#' @param y Outcome variable. If binary will be used for stratification.
#' @param verbose If T will display extra output.
#' @param id Optional cluster grouping, one value per observation. When given,
#' whole clusters are assigned to folds, so no cluster is split between
#' training and validation data. With a binary outcome, clusters are
#' stratified on whether they contain any event, which keeps events balanced
#' across folds as far as whole clusters allow.
#'
#' @return Vector of fold assignments.
#'
#' @importFrom cvTools cvFolds
#' @export
create_cv_folds <- function(v, y, verbose = FALSE, id = NULL) {
  if (!is.null(id)) {
    return(create_cluster_cv_folds(v, y, id, verbose = verbose))
  }
  ys <- unique(y)
  nys <- length(ys)
  nn <- length(y)
  # Binary outcome so we can do stratified fold generation.
  if (nys == 2) {
    out <- rep(NA, nn)
    for (i in 1:nys) {
      # Record how many observations have this Y value.
      n <- sum(y == ys[i])
      folds <- cvTools::cvFolds(n, K = v, R = 1, type = "random")$which

      out[y == ys[i]] <- folds
    }
    if (verbose) {
      cat("Cross-validation fold breakdown:\n")
      print(table(y, "Fold" = out, useNA = "ifany"))
    }
  } else {
    xx <- cvTools::cvFolds(nn, K = v, R = 1, type = "random")$which
    out <- xx
  }
  return(out)
}

#' Cluster-level CV folds
#'
#' Assigns whole clusters to folds. Used by create_cv_folds() when an `id` is
#' given.
#' @noRd
create_cluster_cv_folds <- function(v, y, id, verbose = FALSE) {
  if (length(id) != length(y)) {
    stop("id has ", length(id), " values but y has ", length(y), ".")
  }
  clusters <- unique(id)
  n_clusters <- length(clusters)
  if (n_clusters < v) {
    stop(
      "There are ", n_clusters, " clusters, fewer than n_folds = ", v,
      ". Reduce n_folds."
    )
  }
  cluster_folds <- rep(NA_integer_, n_clusters)

  ys <- unique(y)
  stratified <- FALSE
  if (length(ys) == 2) {
    # A cluster's stratum: whether it contains any event (the larger Y value).
    has_event <- tapply(y == max(ys), id, any)[as.character(clusters)]
    strata_sizes <- table(has_event)
    if (length(strata_sizes) == 2 && all(strata_sizes >= v)) {
      for (s in c(FALSE, TRUE)) {
        in_s <- has_event == s
        cluster_folds[in_s] <- random_folds(sum(in_s), v)
      }
      stratified <- TRUE
    }
  }
  if (!stratified) {
    cluster_folds <- random_folds(n_clusters, v)
  }

  out <- cluster_folds[match(id, clusters)]
  if (verbose) {
    cat("Cluster-grouped cross-validation folds:\n")
    print(table("Fold" = out, useNA = "ifany"))
  }
  out
}

#' Random fold assignment for n units
#'
#' cvTools::cvFolds() returns a random permutation in `subsets` and the fold
#' of each permuted position in `which`; `which` alone is just
#' 1, 2, ..., v, 1, 2, ... in input order. This maps the folds back to the
#' units, so the assignment is random.
#' @noRd
random_folds <- function(n, v) {
  cv <- cvTools::cvFolds(n, K = v, R = 1, type = "random")
  out <- integer(n)
  out[cv$subsets[, 1]] <- cv$which
  out
}
