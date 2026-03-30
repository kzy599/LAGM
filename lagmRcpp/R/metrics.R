calc_group_coancestry <- function(pop = NULL, relationship_matrix = NULL, geno_matrix = NULL) {
  if (!is.null(relationship_matrix)) {
    relationship_matrix <- as.matrix(relationship_matrix)
    return(mean(relationship_matrix[upper.tri(relationship_matrix, diag = TRUE)], na.rm = TRUE))
  }

  if (is.null(geno_matrix)) {
    if (is.null(pop)) {
      stop("Provide either relationship_matrix, geno_matrix, or pop.")
    }
    if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
      stop("AlphaSimR is required when pop is used. Otherwise provide geno_matrix directly.")
    }
    geno_matrix <- AlphaSimR::pullSnpGeno(pop)
  }

  geno_scaled <- scale(as.matrix(geno_matrix), center = TRUE, scale = FALSE)
  gmat <- tcrossprod(geno_scaled) / ncol(geno_scaled)
  mean(gmat[upper.tri(gmat, diag = TRUE)], na.rm = TRUE)
}

calc_inbreeding_proxy <- function(pop = NULL, geno_matrix = NULL) {
  if (is.null(geno_matrix)) {
    if (is.null(pop)) {
      stop("Provide pop or geno_matrix.")
    }
    if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
      stop("AlphaSimR is required when pop is used. Otherwise provide geno_matrix directly.")
    }
    geno_matrix <- AlphaSimR::pullSnpGeno(pop)
  }

  geno_matrix <- as.matrix(geno_matrix)
  mean(rowMeans(geno_matrix == 1, na.rm = TRUE), na.rm = TRUE)
}

calc_normalized_gain <- function(pop = NULL,
                                 founder_mean,
                                 founder_sd,
                                 trait_col = 1L,
                                 values = NULL) {
  if (is.na(founder_sd) || founder_sd <= 0) {
    stop("founder_sd must be positive.")
  }

  if (is.null(values)) {
    if (is.null(pop)) {
      stop("Provide pop or values.")
    }
    values <- pop@gv[, trait_col]
  }

  (mean(values, na.rm = TRUE) - founder_mean) / founder_sd
}
