default_python_binary <- function() {
  configured <- Sys.getenv("LAGM_PYTHON", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }

  for (candidate in c("python3", "python")) {
    if (nzchar(Sys.which(candidate))) {
      return(candidate)
    }
  }

  stop("No Python interpreter was found. Set LAGM_PYTHON or install python/python3.")
}

default_python_script <- function(script_name = "optMatingP.py") {
  installed_script <- system.file("python", script_name, package = "lagm")
  if (nzchar(installed_script)) {
    return(installed_script)
  }

  local_script <- file.path(getwd(), "inst", "python", script_name)
  if (file.exists(local_script)) {
    return(normalizePath(local_script, winslash = "/", mustWork = TRUE))
  }

  stop(sprintf("Python backend script '%s' was not found.", script_name))
}

run_python_lagm <- function(input_json,
                            output_json,
                            python = default_python_binary(),
                            script = default_python_script("optMatingP.py")) {
  python_path <- Sys.which(python)
  if (!nzchar(python_path) && !file.exists(python)) {
    stop(sprintf("Python interpreter '%s' was not found.", python))
  }

  script_path <- normalizePath(script, winslash = "/", mustWork = TRUE)
  cmd <- c(shQuote(script_path), shQuote(input_json), "--output_pairs", shQuote(output_json))
  status <- system2(command = python, args = cmd)

  if (!identical(status, 0L)) {
    stop("Python LAGM backend failed.")
  }

  invisible(output_json)
}
