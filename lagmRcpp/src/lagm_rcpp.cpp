#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <chrono>
#include <limits>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

struct SAResult {
  arma::uvec female_plan;
  arma::uvec male_plan;
  double score;
  double avg_gain;
  double avg_div;
};

inline uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

arma::ivec count_plan_cpp(const arma::uvec& plan, const int n_parent) {
  arma::ivec counts(n_parent, arma::fill::zeros);
  for (unsigned int i = 0; i < plan.n_elem; ++i) {
    counts[plan[i]] += 1;
  }
  return counts;
}

template <typename RNG>
arma::uvec make_feasible_parent_plan_rng(const arma::ivec& min_contrib,
                                         const arma::ivec& max_contrib,
                                         const int n_crosses,
                                         RNG& rng) {
  if (min_contrib.n_elem != max_contrib.n_elem) {
    stop("min_contrib and max_contrib must have the same length.");
  }

  const unsigned int n_parent = min_contrib.n_elem;
  arma::ivec counts(n_parent, arma::fill::zeros);

  int total_cap = 0;
  for (unsigned int i = 0; i < n_parent; ++i) {
    if (min_contrib[i] < 0 || max_contrib[i] < 0) {
      stop("Contribution bounds must be non-negative.");
    }
    if (min_contrib[i] > max_contrib[i]) {
      stop("Each minimum contribution must be <= the matching maximum contribution.");
    }
    total_cap += max_contrib[i];
  }

  if (total_cap < n_crosses) {
    stop("Sum of maximum contributions is smaller than n_crosses.");
  }

  bool built = false;
  for (int attempt = 0; attempt < 400 && !built; ++attempt) {
    counts.zeros();
    int remaining = n_crosses;

    while (remaining > 0) {
      std::vector<unsigned int> feasible;
      feasible.reserve(n_parent);

      for (unsigned int i = 0; i < n_parent; ++i) {
        if (counts[i] >= max_contrib[i]) {
          continue;
        }

        const int block = (counts[i] == 0)
          ? std::max(1, static_cast<int>(min_contrib[i]))
          : 1;

        if (block > remaining) {
          continue;
        }

        int cap_after = 0;
        for (unsigned int j = 0; j < n_parent; ++j) {
          int c = counts[j];
          if (j == i) {
            c += block;
          }
          cap_after += (max_contrib[j] - c);
        }

        if (cap_after >= (remaining - block)) {
          feasible.push_back(i);
        }
      }

      if (feasible.empty()) {
        break;
      }

      std::uniform_int_distribution<int> pick(0, static_cast<int>(feasible.size()) - 1);
      const unsigned int chosen = feasible[pick(rng)];
      const int block = (counts[chosen] == 0)
        ? std::max(1, static_cast<int>(min_contrib[chosen]))
        : 1;

      counts[chosen] += block;
      remaining -= block;
    }

    if (remaining == 0) {
      built = true;
    }
  }

  if (!built) {
    stop("Unable to create a feasible plan under 0-or-[min,max] constraints.");
  }

  arma::uvec plan(n_crosses, arma::fill::zeros);
  int used = 0;
  for (unsigned int i = 0; i < n_parent; ++i) {
    for (int k = 0; k < counts[i]; ++k) {
      plan[used++] = i;
    }
  }

  for (int i = n_crosses - 1; i > 0; --i) {
    std::uniform_int_distribution<int> pick_swap(0, i);
    int j = pick_swap(rng);
    std::swap(plan[i], plan[j]);
  }

  return plan;
}

double evaluate_pair_cpp(const double gain,
                         const double div,
                         const int opt_mode,
                         const double Gmin,
                         const double Gmax,
                         const double Dmin,
                         const double Dmax,
                         const double base_div,
                         const double lookahead_t) {
  if (opt_mode == 1) {
    return gain;
  }
  if (opt_mode == 2) {
    return div;
  }

  const double eps = 1e-12;
  const double div_t = std::pow(div/base_div, lookahead_t);
  const double Dmax_t = std::pow(Dmax/base_div, lookahead_t);
  const double Dmin_t = std::pow(Dmin/base_div, lookahead_t);

  const double ratio_D = std::max((div_t - Dmin_t) / (Dmax_t - Dmin_t + eps), eps);
  const double ratio_G = std::max((gain - Gmin) / (Gmax - Gmin + eps), eps);

  return std::log(ratio_G) + lookahead_t * std::log(ratio_D);

}

// Compute population-level expected heterozygosity (He) from a mating plan.
// For each selected pair k, the offspring allele frequency at locus l is
// p_off_{k,l} = (p_{f_k,l} + p_{m_k,l}) / 2.  The population-mean frequency
// is p_bar_l = mean_k(p_off_{k,l}), and He = mean_l(2 * p_bar_l * (1 - p_bar_l)).
// This captures the between-family variance component (Wahlund: H_T = H_S + 2*Var(p)).
double compute_population_He_from_plan(const arma::uvec& female_plan,
                                       const arma::uvec& male_plan,
                                       const arma::mat& female_geno,
                                       const arma::mat& male_geno) {
  const int n = static_cast<int>(female_plan.n_elem);
  const unsigned int n_loci = female_geno.n_cols;

  arma::rowvec sum_p(n_loci, arma::fill::zeros);
  for (int k = 0; k < n; ++k) {
    sum_p += (female_geno.row(female_plan[k]) + male_geno.row(male_plan[k]));
  }
  // p_bar_l = (sum of (geno_f + geno_m) / 2) / n  =  sum / (2 * n)
  arma::rowvec p_bar = sum_p / (2.0 * n);

  return arma::mean(2.0 * p_bar % (1.0 - p_bar));
}

double evaluate_plan_cpp(const arma::uvec& female_plan,
                         const arma::uvec& male_plan,
                         const arma::mat& gain_mat,
                         const arma::mat& div_mat,
                         int opt_mode,
                         double Gmin,
                         double Gmax,
                         double Dmin,
                         double Dmax,
                         double base_div,
                         double lookahead_t,
                         int diversity_metric = 0,
                         const arma::mat* female_geno_ptr = nullptr,
                         const arma::mat* male_geno_ptr = nullptr,
                         double* avg_gain_out = nullptr,
                         double* avg_div_out = nullptr) {
  double sum_gain = 0.0;
  double sum_div = 0.0;
  const int n = female_plan.n_elem;

  for (int k = 0; k < n; ++k) {
    sum_gain += gain_mat(female_plan[k], male_plan[k]);
    sum_div += div_mat(female_plan[k], male_plan[k]);
  }

  const double avg_gain = sum_gain / n;
  double avg_div;

  if (diversity_metric == 1 && female_geno_ptr != nullptr && male_geno_ptr != nullptr) {
    avg_div = compute_population_He_from_plan(female_plan, male_plan,
                                              *female_geno_ptr, *male_geno_ptr);
  } else {
    avg_div = sum_div / n;
  }

  if (avg_gain_out != nullptr) {
    *avg_gain_out = avg_gain;
  }
  if (avg_div_out != nullptr) {
    *avg_div_out = avg_div;
  }

  return evaluate_pair_cpp(
    avg_gain,
    avg_div,
    opt_mode,
    Gmin,
    Gmax,
    Dmin,
    Dmax,
    base_div,
    lookahead_t
  );
}

template <typename RNG>
SAResult sa_single_run_cpp(const arma::mat& gain_mat,
                           const arma::mat& div_mat,
                           const arma::ivec& female_min,
                           const arma::ivec& female_max,
                           const arma::ivec& male_min,
                           const arma::ivec& male_max,
                           const int n_crosses,
                           const int opt_mode,
                           const double Gmin,
                           const double Gmax,
                           const double Dmin,
                           const double Dmax,
                           const double base_div,
                           const double lookahead_t,
                           const int n_iter,
                           const double swap_prob,
                           const double mutate_female_prob,
                           const double init_prob,
                           const double cooling_rate,
                           const int stop_window,
                           const double stop_eps,
                           const int warmup_iter,
                           const int diversity_metric,
                           const arma::mat* female_geno_ptr,
                           const arma::mat* male_geno_ptr,
                           RNG& rng) {
  const int n_f = gain_mat.n_rows;
  const int n_m = gain_mat.n_cols;

  arma::uvec female_plan = make_feasible_parent_plan_rng(female_min, female_max, n_crosses, rng);
  arma::uvec male_plan = make_feasible_parent_plan_rng(male_min, male_max, n_crosses, rng);
  arma::ivec female_counts = count_plan_cpp(female_plan, n_f);
  arma::ivec male_counts = count_plan_cpp(male_plan, n_m);

  auto propose_mutation = [&](const arma::uvec& in_female_plan,
                              const arma::uvec& in_male_plan,
                              const arma::ivec& in_female_counts,
                              const arma::ivec& in_male_counts,
                              arma::uvec& out_female_plan,
                              arma::uvec& out_male_plan,
                              arma::ivec& out_female_counts,
                              arma::ivec& out_male_counts) -> bool {
    out_female_plan = in_female_plan;
    out_male_plan = in_male_plan;
    out_female_counts = in_female_counts;
    out_male_counts = in_male_counts;

    std::bernoulli_distribution coin_swap(swap_prob);
    std::bernoulli_distribution coin_sex(mutate_female_prob);
    bool use_swap = coin_swap(rng);
    bool mutate_female = coin_sex(rng);
    bool valid_move = true;

    if (use_swap && n_crosses > 1) {
      std::uniform_int_distribution<int> pick_slot(0, n_crosses - 1);
      int i = pick_slot(rng);
      int j = pick_slot(rng);
      while (j == i) {
        j = pick_slot(rng);
      }

      if (mutate_female) {
        std::swap(out_female_plan[i], out_female_plan[j]);
      } else {
        std::swap(out_male_plan[i], out_male_plan[j]);
      }
      return true;
    }

    if (mutate_female) {
      arma::uvec active = arma::find(out_female_counts > 0);
      if (active.n_elem == 0) {
        return false;
      }
      std::uniform_int_distribution<int> pick_active(0, static_cast<int>(active.n_elem) - 1);
      unsigned int A = active[pick_active(rng)];

      std::uniform_int_distribution<int> pick_parent(0, n_f - 1);
      unsigned int B = A;
      for (int t = 0; t < 30 && B == A; ++t) {
        B = static_cast<unsigned int>(pick_parent(rng));
      }
      if (B == A) {
        return false;
      }

      const int cA = out_female_counts[A];
      const int cB = out_female_counts[B];
      const int newB = cB + cA;
      const int minIfSelectedB = std::max(1, static_cast<int>(female_min[B]));
      const bool legalB = newB <= female_max[B] && (newB == 0 || newB >= minIfSelectedB);
      if (!legalB) {
        return false;
      }

      for (unsigned int k = 0; k < out_female_plan.n_elem; ++k) {
        if (out_female_plan[k] == A) {
          out_female_plan[k] = B;
        }
      }
      out_female_counts[A] = 0;
      out_female_counts[B] = newB;
    } else {
      arma::uvec active = arma::find(out_male_counts > 0);
      if (active.n_elem == 0) {
        return false;
      }
      std::uniform_int_distribution<int> pick_active(0, static_cast<int>(active.n_elem) - 1);
      unsigned int A = active[pick_active(rng)];

      std::uniform_int_distribution<int> pick_parent(0, n_m - 1);
      unsigned int B = A;
      for (int t = 0; t < 30 && B == A; ++t) {
        B = static_cast<unsigned int>(pick_parent(rng));
      }
      if (B == A) {
        return false;
      }

      const int cA = out_male_counts[A];
      const int cB = out_male_counts[B];
      const int newB = cB + cA;
      const int minIfSelectedB = std::max(1, static_cast<int>(male_min[B]));
      const bool legalB = newB <= male_max[B] && (newB == 0 || newB >= minIfSelectedB);
      if (!legalB) {
        return false;
      }

      for (unsigned int k = 0; k < out_male_plan.n_elem; ++k) {
        if (out_male_plan[k] == A) {
          out_male_plan[k] = B;
        }
      }
      out_male_counts[A] = 0;
      out_male_counts[B] = newB;
    }

    return valid_move;
  };

  double current_avg_gain = 0.0;
  double current_avg_div = 0.0;
  double current_score = evaluate_plan_cpp(
    female_plan,
    male_plan,
    gain_mat,
    div_mat,
    opt_mode,
    Gmin,
    Gmax,
    Dmin,
    Dmax,
    base_div,
    lookahead_t,
    diversity_metric,
    female_geno_ptr,
    male_geno_ptr,
    &current_avg_gain,
    &current_avg_div
  );

  double best_score = current_score;
  double best_avg_gain = current_avg_gain;
  double best_avg_div = current_avg_div;
  arma::uvec best_female_plan = female_plan;
  arma::uvec best_male_plan = male_plan;

  // --- Auto Warm-up Phase ---
  double sum_worse_delta = 0.0;
  int count_worse = 0;

  for (int w = 0; w < warmup_iter; ++w) {
    arma::uvec trial_female_plan;
    arma::uvec trial_male_plan;
    arma::ivec trial_female_counts;
    arma::ivec trial_male_counts;

    bool valid_move = propose_mutation(
      female_plan,
      male_plan,
      female_counts,
      male_counts,
      trial_female_plan,
      trial_male_plan,
      trial_female_counts,
      trial_male_counts
    );

    if (!valid_move) {
      continue;
    }

    double trial_score = evaluate_plan_cpp(
      trial_female_plan,
      trial_male_plan,
      gain_mat,
      div_mat,
      opt_mode,
      Gmin,
      Gmax,
      Dmin,
      Dmax,
      base_div,
      lookahead_t,
      diversity_metric,
      female_geno_ptr,
      male_geno_ptr
    );

    double delta = trial_score - current_score;
    if (delta < 0.0) {
      sum_worse_delta += delta;
      count_worse++;
    }
  }

  double current_temp = 0.01;
  if (count_worse > 0 && init_prob > 0.0 && init_prob < 1.0) {
    double avg_worse_delta = sum_worse_delta / static_cast<double>(count_worse);
    current_temp = -avg_worse_delta / std::log(init_prob);
    if (!std::isfinite(current_temp) || current_temp <= 0.0) {
      current_temp = 0.01;
    }
  }

  int iter_without_improvement = 0;
  // --- End Warm-up ---

  for (int iter = 0; iter < n_iter; ++iter) {
    arma::uvec trial_female_plan;
    arma::uvec trial_male_plan;
    arma::ivec trial_female_counts;
    arma::ivec trial_male_counts;

    bool valid_move = propose_mutation(
      female_plan,
      male_plan,
      female_counts,
      male_counts,
      trial_female_plan,
      trial_male_plan,
      trial_female_counts,
      trial_male_counts
    );

    if (valid_move) {
      double trial_avg_gain = 0.0;
      double trial_avg_div = 0.0;
      double trial_score = evaluate_plan_cpp(
        trial_female_plan,
        trial_male_plan,
        gain_mat,
        div_mat,
        opt_mode,
        Gmin,
        Gmax,
        Dmin,
        Dmax,
        base_div,
        lookahead_t,
        diversity_metric,
        female_geno_ptr,
        male_geno_ptr,
        &trial_avg_gain,
        &trial_avg_div
      );

      double delta = trial_score - current_score;
      bool accept = delta >= 0.0;

      if (!accept && current_temp > 0.0) {
        std::uniform_real_distribution<double> pick_prob(0.0, 1.0);
        accept = pick_prob(rng) < std::exp(delta / current_temp);
      }

      if (accept) {
        female_plan = trial_female_plan;
        male_plan = trial_male_plan;
        female_counts = trial_female_counts;
        male_counts = trial_male_counts;
        current_score = trial_score;
        current_avg_gain = trial_avg_gain;
        current_avg_div = trial_avg_div;

        if (current_score > best_score) {
          if ((current_score - best_score) > stop_eps) {
            iter_without_improvement = 0;
          }
          best_score = current_score;
          best_avg_gain = current_avg_gain;
          best_avg_div = current_avg_div;
          best_female_plan = female_plan;
          best_male_plan = male_plan;
        }
      }
    }

    iter_without_improvement++;
    if (iter_without_improvement >= stop_window) {
      break;
    }
    current_temp *= cooling_rate;
  }

  return SAResult{best_female_plan, best_male_plan, best_score, best_avg_gain, best_avg_div};
}

// NOTE: Despite the name, this function returns the per-pair observed
// heterozygosity (Ho), computed as Ho = p_f + p_m - 2*p_f*p_m.
// The name is retained for backward API compatibility.
// For the population-level expected heterozygosity (He) under the full
// mating plan, see compute_population_He_from_plan().
// [[Rcpp::export]]
arma::mat compute_expected_heterozygosity_cpp(const arma::mat& female_geno,
                                              const arma::mat& male_geno) {
  const unsigned int n_females = female_geno.n_rows;
  const unsigned int n_males = male_geno.n_rows;
  const unsigned int n_markers = female_geno.n_cols;

  if (male_geno.n_cols != n_markers) {
    stop("Female and male genotype matrices must have the same number of columns.");
  }

  arma::mat out(n_females, n_males, arma::fill::zeros);
  for (unsigned int i = 0; i < n_females; ++i) {
    arma::rowvec pf = female_geno.row(i) / 2.0;
    for (unsigned int j = 0; j < n_males; ++j) {
      arma::rowvec pm = male_geno.row(j) / 2.0;
      arma::rowvec he = pf + pm - 2.0 * (pf % pm);
      out(i, j) = arma::mean(he);
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::mat compute_pair_gain_cpp(const arma::vec& female_ebv,
                                const arma::vec& male_ebv) {
  arma::mat out(female_ebv.n_elem, male_ebv.n_elem, arma::fill::zeros);

  for (unsigned int i = 0; i < female_ebv.n_elem; ++i) {
    for (unsigned int j = 0; j < male_ebv.n_elem; ++j) {
      out(i, j) = 0.5 * (female_ebv[i] + male_ebv[j]);
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::mat compute_pair_relationship_diversity_cpp(const arma::mat& relationship_matrix,
                                                  const arma::uvec& female_index,
                                                  const arma::uvec& male_index) {
  const unsigned int n_females = female_index.n_elem;
  const unsigned int n_males = male_index.n_elem;

  arma::mat out(n_females, n_males, arma::fill::zeros);

  for (unsigned int i = 0; i < n_females; ++i) {
    unsigned int fi = female_index[i];
    if (fi >= relationship_matrix.n_rows) {
      stop("female_index is out of bounds for relationship_matrix.");
    }
    for (unsigned int j = 0; j < n_males; ++j) {
      unsigned int mj = male_index[j];
      if (mj >= relationship_matrix.n_cols) {
        stop("male_index is out of bounds for relationship_matrix.");
      }
      const double rel = relationship_matrix(fi, mj);
      out(i, j) = 1.0 - rel / 2.0;
    }
  }

  return out;
}

// Returns raw gain/diversity matrices only.
// [[Rcpp::export]]
List lagm_score_grid_cpp(const arma::mat& female_geno,
                         const arma::mat& male_geno,
                         const arma::vec& female_ebv,
                         const arma::vec& male_ebv) {
  arma::mat diversity = compute_expected_heterozygosity_cpp(female_geno, male_geno);
  arma::mat gain = compute_pair_gain_cpp(female_ebv, male_ebv);

  return List::create(
    Named("expected_diversity") = diversity,
    Named("expected_gain") = gain
  );
}

// Returns raw gain/diversity matrices only.
// [[Rcpp::export]]
List lagm_relationship_score_grid_cpp(const arma::mat& relationship_matrix,
                                      const arma::uvec& female_index,
                                      const arma::uvec& male_index,
                                      const arma::vec& female_ebv,
                                      const arma::vec& male_ebv) {
  arma::mat diversity = compute_pair_relationship_diversity_cpp(
    relationship_matrix,
    female_index,
    male_index
  );
  arma::mat gain = compute_pair_gain_cpp(female_ebv, male_ebv);

  return List::create(
    Named("expected_diversity") = diversity,
    Named("expected_gain") = gain
  );
}

// [[Rcpp::export]]
List optimize_mating_plan_cpp(const arma::mat& gain_mat,
                              const arma::mat& div_mat,
                              const IntegerVector& female_min,
                              const IntegerVector& female_max,
                              const IntegerVector& male_min,
                              const IntegerVector& male_max,
                              const int n_crosses,
                              const int opt_mode = 3,
                              const double Gmin = 0.0,
                              const double Gmax = 1.0,
                              const double Dmin = 0.0,
                              const double Dmax = 1.0,
                              const double base_div = 1.0,
                              const double lookahead_t = 1.0,
                              const int n_iter = 2000,
                              const double swap_prob = 0.2,
                              const double mutate_female_prob = 0.5,
                              const double init_prob = 0.8,
                              const double cooling_rate = 0.995,
                              const int stop_window = 1000,
                              const double stop_eps = 1e-8,
                              const int warmup_iter = 100,
                              const int n_pop = 50,
                              const int n_threads = 4,
                              const int diversity_metric = 1,
                              Rcpp::Nullable<Rcpp::NumericMatrix> female_geno = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericMatrix> male_geno = R_NilValue) {
  const int n_f = gain_mat.n_rows;
  const int n_m = gain_mat.n_cols;

  if (div_mat.n_rows != n_f || div_mat.n_cols != n_m) {
    stop("gain_mat and div_mat must have identical dimensions.");
  }

  if (female_min.size() != n_f || female_max.size() != n_f ||
      male_min.size() != n_m || male_max.size() != n_m) {
    stop("Contribution bound vectors must match gain/div matrix dimensions.");
  }
  if (n_crosses <= 0) {
    stop("n_crosses must be positive.");
  }
  if (n_pop <= 0) {
    stop("n_pop must be positive.");
  }
  if (n_threads <= 0) {
    stop("n_threads must be positive.");
  }
  if (swap_prob < 0.0 || swap_prob > 1.0) {
    stop("swap_prob must be in [0,1].");
  }
  if (mutate_female_prob < 0.0 || mutate_female_prob > 1.0) {
    stop("mutate_female_prob must be in [0,1].");
  }

  // Resolve optional genotype matrices for pop_He mode
  arma::mat female_geno_arma;
  arma::mat male_geno_arma;
  const arma::mat* female_geno_ptr = nullptr;
  const arma::mat* male_geno_ptr = nullptr;

  if (diversity_metric == 1) {
    if (female_geno.isNull() || male_geno.isNull()) {
      stop("diversity_metric = 1 (pop_He) requires female_geno and male_geno. "
           "Only genomic mode supports pop_He.");
    }
    female_geno_arma = as<arma::mat>(Rcpp::NumericMatrix(female_geno));
    male_geno_arma   = as<arma::mat>(Rcpp::NumericMatrix(male_geno));
    if (static_cast<int>(female_geno_arma.n_rows) != n_f) {
      stop("female_geno must have the same number of rows as gain_mat.");
    }
    if (static_cast<int>(male_geno_arma.n_rows) != n_m) {
      stop("male_geno must have the same number of rows as gain_mat (columns).");
    }
    female_geno_ptr = &female_geno_arma;
    male_geno_ptr   = &male_geno_arma;
  }

  arma::ivec female_min_arma = as<arma::ivec>(female_min);
  arma::ivec female_max_arma = as<arma::ivec>(female_max);
  arma::ivec male_min_arma = as<arma::ivec>(male_min);
  arma::ivec male_max_arma = as<arma::ivec>(male_max);

  std::vector<arma::uvec> female_plans(n_pop);
  std::vector<arma::uvec> male_plans(n_pop);
  std::vector<double> scores(n_pop, -std::numeric_limits<double>::infinity());
  std::vector<double> avg_gains(n_pop, NA_REAL);
  std::vector<double> avg_divs(n_pop, NA_REAL);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(static)
#endif
  for (int p = 0; p < n_pop; ++p) {
    uint64_t seed64 = static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );
#ifdef _OPENMP
    seed64 ^= static_cast<uint64_t>(omp_get_thread_num() + 1) * 0x9e3779b97f4a7c15ULL;
#endif
    seed64 ^= static_cast<uint64_t>(p + 1) * 0xbf58476d1ce4e5b9ULL;
    seed64 = splitmix64(seed64);

    std::mt19937 rng(static_cast<uint32_t>(seed64 & 0xffffffffULL));

    SAResult run = sa_single_run_cpp(
      gain_mat,
      div_mat,
      female_min_arma,
      female_max_arma,
      male_min_arma,
      male_max_arma,
      n_crosses,
      opt_mode,
      Gmin,
      Gmax,
      Dmin,
      Dmax,
      base_div,
      lookahead_t,
      n_iter,
      swap_prob,
      mutate_female_prob,
      init_prob,
      cooling_rate,
      stop_window,
      stop_eps,
      warmup_iter,
      diversity_metric,
      female_geno_ptr,
      male_geno_ptr,
      rng
    );

    female_plans[p] = run.female_plan;
    male_plans[p] = run.male_plan;
    scores[p] = run.score;
    avg_gains[p] = run.avg_gain;
    avg_divs[p] = run.avg_div;
  }

  int best_idx = 0;
  for (int p = 1; p < n_pop; ++p) {
    if (scores[p] > scores[best_idx]) {
      best_idx = p;
    }
  }

  arma::uvec best_female_plan = female_plans[best_idx];
  arma::uvec best_male_plan = male_plans[best_idx];
  double best_score = scores[best_idx];

  IntegerVector female_index(best_female_plan.begin(), best_female_plan.end());
  IntegerVector male_index(best_male_plan.begin(), best_male_plan.end());
  female_index = female_index + 1;
  male_index = male_index + 1;

  NumericVector score(best_female_plan.n_elem);
  NumericVector pair_gain(best_female_plan.n_elem);
  NumericVector pair_div(best_female_plan.n_elem);

  for (unsigned int k = 0; k < best_female_plan.n_elem; ++k) {
    pair_gain[k] = gain_mat(best_female_plan[k], best_male_plan[k]);
    pair_div[k] = div_mat(best_female_plan[k], best_male_plan[k]);
    score[k] = evaluate_pair_cpp(
      pair_gain[k],
      pair_div[k],
      opt_mode,
      Gmin,
      Gmax,
      Dmin,
      Dmax,
      base_div,
      lookahead_t
    );
  }

  return List::create(
    Named("female_index") = female_index,
    Named("male_index") = male_index,
    Named("score") = score,
    Named("pair_gain") = pair_gain,
    Named("pair_diversity") = pair_div,
    Named("objective_sum") = best_score,
    Named("avg_gain") = avg_gains[best_idx],
    Named("avg_diversity") = avg_divs[best_idx],
    Named("female_counts") = wrap(count_plan_cpp(best_female_plan, n_f)),
    Named("male_counts") = wrap(count_plan_cpp(best_male_plan, n_m))
  );
}