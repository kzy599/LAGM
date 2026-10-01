# Compile the actual optimizer plus small adapters, with deterministic RNGs
# only in the test/validation process. The installed API is not instrumented.
load_pop_engine <- function(root, source = file.path(root, "src", "lagm_rcpp.cpp")) {
  code <- readLines(source)
  joint <- any(grepl("Rcpp::RObject pop_epsilon", code, fixed = TRUE))
  env <- new.env()
  Rcpp::sourceCpp(code = paste(c(
    if (joint) "#define LAGM_JOINT_POP",
    code, readLines(file.path(root, "tests", "fixtures", "pop-engine.cpp"))
  ), collapse = "\n"), env = env)
  env
}

pop_engine <- local({
  engine <- NULL
  function() {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
    if (!file.exists(file.path(root, "src", "lagm_rcpp.cpp")) &&
        dir.exists(file.path(root, "00_pkg_src", "lagm"))) {
      root <- file.path(root, "00_pkg_src", "lagm")
    }
    testthat::skip_if_not(file.exists(file.path(root, "src", "lagm_rcpp.cpp")),
                          "C++ engine probes require the source checkout")
    if (is.null(engine)) engine <<- load_pop_engine(root)
    engine
  }
})

pop_fixture <- function(metric = 1L) {
  fg <- matrix(c(0, 0, 1, 0, 1, 0, 1, 1, 0), 3, byrow = TRUE)
  mg <- matrix(c(1, 1, 0, 1, 0, 1, 0, 0, 1), 3, byrow = TRUE)
  K <- diag(6)
  K[1:3, 4:6] <- matrix(c(0.8, 0.1, 0.2, 0.2, 0.7, 0.1,
                          0.1, 0.2, 0.6), 3)
  K[4:6, 1:3] <- t(K[1:3, 4:6])
  div <- if (metric %in% c(0L, 1L)) {
    lagm::compute_expected_heterozygosity_cpp(fg, mg)
  } else {
    1 - K[1:3, 4:6] / 2
  }
  list(gain = outer(c(0.3, 0.8, 1.4), c(0.2, 0.6, 1.2), "+") / 2,
       div = div, metric = metric, fg = fg, mg = mg, K = K)
}

pop_api_args <- function(mode = "genomic") {
  fixture <- pop_fixture()
  ids <- paste0("p", 1:6)
  geno <- rbind(fixture$fg, fixture$mg)
  rownames(geno) <- ids
  K <- fixture$K
  dimnames(K) <- list(ids, ids)
  list(individual_ids = ids, female_ids = ids[1:3], male_ids = ids[4:6],
       ebv_vector = c(0.3, 0.8, 1.4, 0.2, 0.6, 1.2),
       n_crosses = 4L, lookahead_generations = 1L,
       female_max = rep(2L, 3), male_max = rep(2L, 3),
       diversity_mode = mode, diversity_level = "pop",
       geno_matrix = geno, relationship_matrix = K,
       n_iter = 40L, warmup_iter = 10L, n_pop = 2L, n_threads = 1L)
}
