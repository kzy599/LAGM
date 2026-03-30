default_python_binary <- function() {
  stop("Python backend has been removed. Use lagm_plan() or lagm_mating() instead.")
}

default_python_script <- function(script_name = "optMatingP.py") {
  stop("Python backend has been removed. Use lagm_plan() or lagm_mating() instead.")
}

run_python_lagm <- function(input_json, output_json, python = NULL, script = NULL) {
  stop("Python backend has been removed. Use lagm_plan() or lagm_mating() instead.")
}

safe_write_json <- function(x, path) {
  stop("JSON/Python backend helpers are no longer used in the refactored package.")
}

read_python_mating_plan <- function(path) {
  stop("JSON/Python backend helpers are no longer used in the refactored package.")
}
