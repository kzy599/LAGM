#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

// Test-only observers: no RNG draws or writes to the optimizer's state.
std::vector<List> sa_trials;
std::vector<List> sa_temperatures;

void record_sa_trial(const arma::uvec& fp, const arma::uvec& mp,
                     const arma::uvec& trial_fp, const arma::uvec& trial_mp,
                     const arma::rowvec& sum_p, const arma::rowvec& trial_sum_p,
                     double score, double trial_score, bool warmup) {
  sa_trials.push_back(List::create(
    _["fp"] = fp, _["mp"] = mp, _["trial_fp"] = trial_fp, _["trial_mp"] = trial_mp,
    _["sum_p"] = sum_p, _["trial_sum_p"] = trial_sum_p,
    _["score"] = score, _["trial_score"] = trial_score,
    _["warmup"] = warmup, _["accepted"] = false));
}

void record_sa_temperature(double sum, int count, double temperature) {
  sa_temperatures.push_back(List::create(
    _["sum"] = sum, _["count"] = count, _["temperature"] = temperature));
}

// [[Rcpp::export]]
void reset_sa_trace() {
  sa_trials.clear();
  sa_temperatures.clear();
}

// [[Rcpp::export]]
List inspect_sa_trace() {
  return List::create(_["trials"] = wrap(sa_trials),
                      _["temperatures"] = wrap(sa_temperatures));
}

// The loader inserts the actual production calibration block here.
// [[Rcpp::export]]
double inspect_warmup(NumericVector deltas, double init_prob) {
  double sum_worse_delta = 0.0;
  int count_worse = 0;
  for (double delta : deltas) {
    if (delta < 0.0) {
      sum_worse_delta += delta;
      ++count_worse;
    }
  }
  // WARMUP_CALIBRATION
  return current_temp;
}
