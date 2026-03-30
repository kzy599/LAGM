`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Validate that IDs are unique and non-missing.
validate_unique_ids <- function(ids, arg = "ids") {
  ids <- as.character(ids)
  if (length(ids) == 0L) {
    stop(sprintf("%s must not be empty.", arg))
  }
  if (anyNA(ids) || any(ids == "")) {
    stop(sprintf("%s must not contain NA/empty values.", arg))
  }
  if (anyDuplicated(ids)) {
    stop(sprintf("%s must be unique.", arg))
  }
  ids
}

calc_expected_heterozygosity <- function(geno_matrix) {
  geno_matrix <- as.matrix(geno_matrix)
  if (!is.numeric(geno_matrix)) {
    stop("geno_matrix must be numeric.")
  }
  if (nrow(geno_matrix) == 0L || ncol(geno_matrix) == 0L) {
    stop("geno_matrix must have at least one row and one column.")
  }

  mean(colMeans(geno_matrix == 1, na.rm = TRUE))
}

calc_relationship_diversity <- function(relationship_matrix, female_ids = NULL, male_ids = NULL) {
  relationship_matrix <- as.matrix(relationship_matrix)
  if (!is.numeric(relationship_matrix)) {
    stop("relationship_matrix must be numeric.")
  }

  if (is.null(female_ids) || is.null(male_ids)) {
    rel_vals <- relationship_matrix[upper.tri(relationship_matrix, diag = TRUE)]
  } else {
    if (is.null(rownames(relationship_matrix)) || is.null(colnames(relationship_matrix))) {
      stop("relationship_matrix must have row/column names when female_ids/male_ids are provided.")
    }
    female_idx <- match(female_ids, rownames(relationship_matrix))
    male_idx <- match(male_ids, colnames(relationship_matrix))
    if (anyNA(female_idx) || anyNA(male_idx)) {
      stop("female_ids or male_ids were not found in relationship_matrix dimnames.")
    }
    rel_vals <- relationship_matrix[female_idx, male_idx, drop = FALSE]
  }

  mean(pmax(1e-12, 1 - rel_vals / 2), na.rm = TRUE)
}

make_lookahead_objective <- function(expected_gain, expected_diversity, base_diversity, lookahead_generations) {
  if (!all(dim(expected_gain) == dim(expected_diversity))) {
    stop("expected_gain and expected_diversity must have identical dimensions.")
  }

  if (is.null(base_diversity) || is.na(base_diversity) || base_diversity <= 0) {
    return(expected_gain * expected_diversity ^ lookahead_generations)
  }

  expected_gain * (expected_diversity / base_diversity) ^ lookahead_generations
}

build_lagm_input <- function(individual_ids,
                             female_ids,
                             male_ids,
                             ebv_vector,
                             gen,
                             diversity_mode = c("genomic", "relationship"),
                             base_diversity = NULL,
                             geno_matrix = NULL,
                             relationship_matrix = NULL) {
  diversity_mode <- match.arg(diversity_mode)

  individual_ids <- validate_unique_ids(individual_ids, "individual_ids")
  female_ids <- validate_unique_ids(female_ids, "female_ids")
  male_ids <- validate_unique_ids(male_ids, "male_ids")

  if (!all(female_ids %in% individual_ids)) {
    stop("All female_ids must be contained in individual_ids.")
  }
  if (!all(male_ids %in% individual_ids)) {
    stop("All male_ids must be contained in individual_ids.")
  }

  ebv_vector <- as.numeric(ebv_vector)
  if (length(ebv_vector) != length(individual_ids)) {
    stop("ebv_vector length must equal length(individual_ids).")
  }
  if (anyNA(ebv_vector)) {
    stop("ebv_vector must not contain missing values.")
  }

  female_idx <- match(female_ids, individual_ids)
  male_idx <- match(male_ids, individual_ids)

  out <- list(
    individual_ids = individual_ids,
    ebv_vector = ebv_vector,
    female_ids = female_ids,
    male_ids = male_ids,
    female_index = female_idx,
    male_index = male_idx,
    female_ebv = ebv_vector[female_idx],
    male_ebv = ebv_vector[male_idx],
    gen = gen,
    diversity_mode = diversity_mode
  )

  if (identical(diversity_mode, "genomic")) {
    if (is.null(geno_matrix)) {
      stop("geno_matrix must be supplied when diversity_mode = 'genomic'.")
    }

    geno_matrix <- as.matrix(geno_matrix)
    if (!is.numeric(geno_matrix)) {
      stop("geno_matrix must be numeric.")
    }
    if (nrow(geno_matrix) != length(individual_ids)) {
      stop("geno_matrix must have one row per individual in individual_ids.")
    }

    geno_scaled <- scale(geno_matrix, center = TRUE, scale = FALSE)
    kinship_matrix <- tcrossprod(geno_scaled) / ncol(geno_matrix)

    out$geno_matrix <- geno_matrix
    out$kinship_matrix <- kinship_matrix
    out$female_geno <- geno_matrix[female_idx, , drop = FALSE]
    out$male_geno <- geno_matrix[male_idx, , drop = FALSE]
    out$base_diversity <- base_diversity
  } else {
    if (is.null(relationship_matrix)) {
      stop("relationship_matrix must be supplied when diversity_mode = 'relationship'.")
    }

    relationship_matrix <- as.matrix(relationship_matrix)
    if (!is.numeric(relationship_matrix)) {
      stop("relationship_matrix must be numeric.")
    }
    if (nrow(relationship_matrix) != length(individual_ids) ||
        ncol(relationship_matrix) != length(individual_ids)) {
      stop("relationship_matrix must be square with dimension length(individual_ids).")
    }

    out$relationship_matrix <- relationship_matrix
    out$base_diversity <- base_diversity
  }

  out
}

safe_write_json <- function(x, path) {
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE)
  invisible(path)
}

read_python_mating_plan <- function(path) {
  out <- jsonlite::fromJSON(path)
  pairs <- out$breeding_pairs
  data.table::as.data.table(pairs)
}
