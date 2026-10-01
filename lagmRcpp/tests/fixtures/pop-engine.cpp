// Test/validation adapters, appended to the production source by sourceCpp.
// No changes to production seeding, neighborhoods, or annealing are needed.

// [[Rcpp::export]]
List inspect_plan(const arma::mat& gain, const arma::mat& div,
                  const arma::uvec& fp, const arma::uvec& mp,
                  int metric, const arma::mat& fg, const arma::mat& mg,
                  const arma::mat& K, double epsilon = 0.005,
                  int mode = 3, double Gmin = 0.0, double Gmax = 2.0,
                  double Dmin = -2.0, double Dmax = 1.0,
                  double base = 1.0, double horizon = 1.0) {
  arma::vec x = arma::conv_to<arma::vec>::from(
    arma::join_cols(count_plan_cpp(fp, gain.n_rows),
                    count_plan_cpp(mp, gain.n_cols)));
  double avg_gain, avg_div;
  double S = evaluate_plan_cpp(fp, mp, gain, div, mode,
    Gmin, Gmax, Dmin, Dmax, base, horizon, metric,
    &fg, &mg, nullptr, &K, &x, &avg_gain, &avg_div
#ifdef LAGM_JOINT_POP
    , epsilon, div.min(), div.max()
#endif
  );
  double q = 0.0;
  for (unsigned int k = 0; k < fp.n_elem; ++k) q += div(fp[k], mp[k]);
  q /= fp.n_elem;
  double Q = div.max() > div.min()
    ? std::max(0.0, std::min(1.0, (q - div.min()) / (div.max() - div.min())))
    : 0.0;
  double J = evaluate_pair_cpp(avg_gain, avg_div, mode,
    Gmin, Gmax, Dmin, Dmax, base, horizon);
  return List::create(_["J"] = J, _["q"] = q, _["Q"] = Q, _["S"] = S,
                      _["gain"] = avg_gain, _["diversity"] = avg_div);
}

// [[Rcpp::export]]
List seeded_search(const arma::mat& gain, const arma::mat& div,
                   int metric, const arma::mat& fg, const arma::mat& mg,
                   const arma::mat& K, int seed = 123, double epsilon = 0.005,
                   int mode = 3, double Gmin = 0.0, double Gmax = 2.0,
                   double Dmin = -2.0, double Dmax = 1.0,
                   int iterations = 200, int warmup = 40,
                   double swap = 0.2, int restarts = 3) {
  arma::ivec fmin(gain.n_rows, arma::fill::zeros);
  arma::ivec mmin(gain.n_cols, arma::fill::zeros);
  arma::ivec fmax(gain.n_rows); fmax.fill(2);
  arma::ivec mmax(gain.n_cols); mmax.fill(2);
  SAResult best;
  best.score = -std::numeric_limits<double>::infinity();
  for (int r = 0; r < restarts; ++r) {
    std::mt19937 rng(seed + r);
    SAResult run = sa_single_run_cpp(gain, div, fmin, fmax, mmin, mmax,
      4, mode, Gmin, Gmax, Dmin, Dmax, 1.0, 1.0,
      iterations, swap, 0.5, 0.8, 0.995, 1000, 1e-8, warmup, metric,
      &fg, &mg, &K, rng
#ifdef LAGM_JOINT_POP
      , epsilon, div.min(), div.max()
#endif
    );
    if (run.score > best.score) best = run;
  }
  return List::create(_["fp"] = best.female_plan, _["mp"] = best.male_plan,
                      _["score"] = best.score, _["gain"] = best.avg_gain,
                      _["diversity"] = best.avg_div);
}
