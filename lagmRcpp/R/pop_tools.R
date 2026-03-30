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
