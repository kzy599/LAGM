# AlphaSimR convenience helpers (not required by the core optimizer).

create_candidate_pool <- function(pop, n_females = 4L, n_males = 2L, use = "ebv") {
  females <- AlphaSimR::selectWithinFam(pop, nInd = n_females, use = use, sex = "F", trait = 1)
  males <- AlphaSimR::selectWithinFam(pop, nInd = n_males, use = use, sex = "M", trait = 1)

  list(
    females = females,
    males = males,
    candidate = c(females, males)
  )
}

create_candidate <- function(pop, n_females = 4L, n_males = 2L, use = "ebv") {
  pool <- create_candidate_pool(pop, n_females = n_females, n_males = n_males, use = use)
  pool$candidate
}

create_ped <- function(pop) {
  ped <- data.table::data.table(id = pop@id, sire = pop@father, dam = pop@mother)
  ped
}

make_bv_ebv <- function(pop, trait = 1L) {
  pop@ebv <- matrix(pop@gv[, trait], ncol = 1)
  pop
}

calc_a_matrix <- function(ped, keep) {
  if (!requireNamespace("optiSel", quietly = TRUE)) {
    stop("Package 'optiSel' is required for pedigree relationship matrices.")
  }

  optiSel::prePed(as.data.frame(ped), keep = keep)
  amat <- as.matrix(optiSel::makeA(as.data.frame(ped)))
  row_idx <- match(keep, rownames(amat))
  col_idx <- match(keep, colnames(amat))
  amat[row_idx, col_idx, drop = FALSE]
}

inject_gebv_into_pop <- function(pop, gebv_dt, id_col = "ID", value_col = "GEBV") {
  if (!data.table::is.data.table(gebv_dt)) {
    gebv_dt <- data.table::as.data.table(gebv_dt)
  }

  if (!id_col %in% names(gebv_dt) || !value_col %in% names(gebv_dt)) {
    stop("gebv_dt must contain the requested ID and value columns.")
  }

  idx <- match(pop@id, gebv_dt[[id_col]])
  if (anyNA(idx)) {
    missing_ids <- pop@id[is.na(idx)]
    stop(sprintf("Missing GEBV values for %d IDs: %s", length(missing_ids), paste(utils::head(missing_ids, 10), collapse = ", ")))
  }

  gebv <- gebv_dt[[value_col]][idx]
  if (length(gebv) != pop@nInd) {
    stop("Injected GEBV length does not match population size.")
  }
  if (anyNA(gebv)) {
    stop("Injected GEBV contains missing values after ID matching.")
  }

  pop@ebv <- matrix(gebv, ncol = 1)
  pop
}

compute_candidate_gebv <- function(pheno_file,
                                   trait_pos,
                                   addG,
                                   domG = "",
                                   out_prefix = tempfile(pattern = "hiblup_"),
                                   hiblup_runner = NULL) {
  if (is.null(hiblup_runner)) {
    if (!exists("run_hiblup", mode = "function")) {
      stop("run_hiblup() is not available. Source legacy/utils_addition.r or provide hiblup_runner.")
    }
    hiblup_runner <- get("run_hiblup", mode = "function")
  }

  res <- hiblup_runner(
    phename = pheno_file,
    trait_pos = trait_pos,
    addG = addG,
    domG = domG,
    out_prefix = out_prefix
  )

  model_idx <- which(vapply(res, function(x) identical(x$comp, "gv"), logical(1)))
  if (length(model_idx) == 0L) {
    stop("run_hiblup() did not return a gv component.")
  }

  gebv_dt <- data.table::copy(res[[model_idx[1]]]$dt)
  data.table::setnames(gebv_dt, old = names(gebv_dt)[names(gebv_dt) == "Prediction"], new = "GEBV")
  gebv_dt
}

# VanRaden Method 2 genomic relationship matrix.
#
# G = Z W Z' / m_eff, where Z = M - 2p, columns of Z are scaled by
# 1 / sqrt(2 p (1 - p)), and m_eff is the number of polymorphic loci kept
# after MAF filtering.  Loci with heterozygosity below 2*min_maf*(1-min_maf)
# are dropped to avoid division by ~0.
#
# `geno_matrix` should have rows = individuals (rownames are individual IDs)
# and columns = SNP markers, coded as 0/1/2 dosages of the alternate allele.
# Returns a square matrix with row/column names equal to rownames(geno_matrix).
compute_vr2_grm <- function(geno_matrix, min_maf = 1e-3) {
  geno <- as.matrix(geno_matrix)
  if (!is.numeric(geno)) {
    stop("geno_matrix must be numeric (0/1/2 dosages).")
  }
  if (nrow(geno) == 0L || ncol(geno) == 0L) {
    stop("geno_matrix must have at least one row and one column.")
  }

  p <- colMeans(geno, na.rm = TRUE) / 2
  het <- 2 * p * (1 - p)
  keep <- het > 2 * min_maf * (1 - min_maf)
  if (sum(keep) < 10L) {
    stop("Too few polymorphic loci for VR2 GRM (need at least 10 after MAF filter).")
  }

  geno_sub <- geno[, keep, drop = FALSE]
  p_sub    <- p[keep]
  het_sub  <- 2 * p_sub * (1 - p_sub)

  Z <- sweep(geno_sub, 2, 2 * p_sub, FUN = "-")
  Z <- sweep(Z,        2, sqrt(het_sub), FUN = "/")
  G <- tcrossprod(Z) / sum(keep)

  ids <- rownames(geno_matrix)
  if (!is.null(ids)) {
    rownames(G) <- ids
    colnames(G) <- ids
  }
  G
}
