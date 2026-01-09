#' =============================================================================
#' BAYESIAN V DE AIKEN - BayesAiken Package
#' =============================================================================
#'
#' Functions for Bayesian estimation of Aiken's V coefficient using two models:
#' - Beta-Binomial (standard approach)
#' - Dirichlet-Multinomial (conservative approach, recommended for small panels)
#'
#' Based on: Ventura-León, J. (2025). A Dirichlet-Multinomial Framework for
#' Bayesian Aiken's V in Content Validity Studies with Small Expert Panels.
#'
#' =============================================================================

# =============================================================================
# MAIN FUNCTION: V_aiken()
# =============================================================================

#' Bayesian Aiken's V Coefficient
#'
#' Calculates Aiken's V coefficient with Bayesian credible intervals using
#' either the Beta-Binomial or Dirichlet-Multinomial model.
#'
#' @param ratings Vector of judge ratings or frequency vector (see details)
#' @param l Minimum scale value (default = 0)
#' @param s Maximum scale value (default = 3, so scale is 0-3)
#' @param model Model to use: "dirichlet" (default, conservative) or "beta"
#' @param prior_alpha Prior concentration parameter (default = 1 for uniform)
#' @param cred_level Credibility level for intervals (default = 0.95)
#' @param n_samples Number of Monte Carlo samples for Dirichlet model (default = 10000)
#' @param threshold Decision threshold for probability calculation (default = 0.70)
#' @param seed Random seed for reproducibility (default = 2025)
#' @param input_type Type of input: "ratings" (individual ratings) or "counts" (frequency vector)
#' @param verbose Print detailed results (default = TRUE)
#'
#' @return Object of class 'bayes_V' containing:
#' \itemize{
#'   \item model: Model used ("Beta-Binomial" or "Dirichlet-Multinomial")
#'   \item n_judges: Number of judges
#'   \item V_classic: Classical V estimate
#'   \item V_mean: Posterior mean
#'   \item V_median: Posterior median
#'   \item V_sd: Posterior standard deviation
#'   \item HDI: Highest Density Interval
#'   \item ETI: Equal-Tailed Interval
#'   \item prob_threshold: P(V > threshold | data)
#'   \item counts: Category frequency vector
#'   \item posterior_params: Posterior distribution parameters
#' }
#'
#' @details
#' The function accepts two input formats:
#'
#' **Ratings format** (input_type = "ratings"):
#' Vector of individual judge ratings, e.g., c(3, 3, 2, 3, 3) for 5 judges.
#'
#' **Counts format** (input_type = "counts"):
#' Frequency vector where position i represents count for category i-1.
#' E.g., c(0, 0, 1, 4) means 0 judges chose 0, 0 chose 1, 1 chose 2, 4 chose 3.
#'
#' **Model comparison:**
#' - Beta-Binomial: Treats each rating point as a binary trial. Can be
#'   optimistic with small panels as it inflates effective sample size.
#' - Dirichlet-Multinomial: Models category frequencies directly. More
#'   conservative and appropriate for small panels (n < 10).
#'
#' @examples
#' # Example from the article: Item 1 (Relevance)
#' # 5 judges, scale 0-3, ratings: one judge chose 2, four chose 3
#'
#' # Using ratings vector
#' ratings <- c(2, 3, 3, 3, 3)
#' result <- V_aiken(ratings, l = 0, s = 3)
#'
#' # Using frequency counts
#' counts <- c(0, 0, 1, 4)  # n0=0, n1=0, n2=1, n3=4
#' result <- V_aiken(counts, l = 0, s = 3, input_type = "counts")
#'
#' # Compare models
#' result_dir <- V_aiken(ratings, model = "dirichlet")
#' result_beta <- V_aiken(ratings, model = "beta")
#'
#' @references
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability
#' and validity of ratings. Educational and Psychological Measurement, 45(1), 131-142.
#'
#' Ventura-León, J. (2025). A Dirichlet-Multinomial Framework for Bayesian
#' Aiken's V in Content Validity Studies with Small Expert Panels.
#'
#' @export
V_aiken <- function(ratings,
                    l = 0,
                    s = 3,
                    model = c("dirichlet", "beta"),
                    prior_alpha = 1,
                    cred_level = 0.95,
                    n_samples = 10000,
                    threshold = 0.70,
                    seed = 2025,
                    input_type = c("ratings", "counts"),
                    verbose = TRUE) {

  # Match arguments

  model <- match.arg(model)
  input_type <- match.arg(input_type)

  # Scale parameters
  k <- s - l
  n_categories <- k + 1

  # Process input
  if (input_type == "ratings") {
    ratings <- ratings[!is.na(ratings)]
    n <- length(ratings)

    if (n == 0) stop("No valid ratings provided")
    if (any(ratings < l) || any(ratings > s)) {
      stop(sprintf("Ratings must be between %d and %d", l, s))
    }

    # Convert to counts
    ratings_idx <- ratings - l
    counts <- tabulate(ratings_idx + 1, nbins = n_categories)

  } else {
    # Input is already counts
    counts <- ratings
    if (length(counts) != n_categories) {
      stop(sprintf("Counts vector must have %d elements (categories 0 to %d)", n_categories, k))
    }
    n <- sum(counts)
  }

  # Classical V
  category_values <- 0:k
  V_classic <- sum(category_values * counts) / (n * k)

  # Calculate based on model
  if (model == "dirichlet") {
    result <- .V_dirichlet(counts, k, prior_alpha, cred_level, n_samples,
                           threshold, seed, n, V_classic)
  } else {
    result <- .V_beta_binomial(counts, k, prior_alpha, cred_level,
                               threshold, n, V_classic)
  }

  # Add common elements
  result$n_judges <- n
  result$V_classic <- round(V_classic, 4)
  result$counts <- counts
  result$scale <- c(l = l, s = s, k = k)
  result$threshold <- threshold
  result$cred_level <- cred_level

  class(result) <- c("bayes_V", "list")

  # Print if verbose

  if (verbose) {
    print(result)
  }

  invisible(result)
}


# =============================================================================
# DIRICHLET-MULTINOMIAL MODEL
# =============================================================================

#' @noRd
.V_dirichlet <- function(counts, k, prior_alpha, cred_level, n_samples,
                         threshold, seed, n, V_classic) {

  n_categories <- k + 1
  category_values <- 0:k

  # Prior and posterior
  if (length(prior_alpha) == 1) {
    prior_alpha_vec <- rep(prior_alpha, n_categories)
  } else {
    prior_alpha_vec <- prior_alpha
  }
  post_alpha_vec <- prior_alpha_vec + counts

  # Monte Carlo sampling
  if (!is.null(seed)) set.seed(seed)

  V_samples <- numeric(n_samples)
  for (i in 1:n_samples) {
    gamma_samples <- stats::rgamma(n_categories, shape = post_alpha_vec, rate = 1)
    pi_samples <- gamma_samples / sum(gamma_samples)
    V_samples[i] <- sum(pi_samples * category_values) / k
  }

  # Posterior statistics
  V_mean <- mean(V_samples)
  V_median <- stats::median(V_samples)
  V_sd <- stats::sd(V_samples)

  # HDI (Highest Density Interval)
  sorted_samples <- sort(V_samples)
  ci_mass <- floor(cred_level * n_samples)
  n_cis <- n_samples - ci_mass
  ci_widths <- sorted_samples[(ci_mass + 1):n_samples] - sorted_samples[1:n_cis]
  best_ci_idx <- which.min(ci_widths)
  HDI <- c(sorted_samples[best_ci_idx], sorted_samples[best_ci_idx + ci_mass])

  # ETI (Equal-Tailed Interval)
  alpha_ci <- 1 - cred_level
  ETI <- stats::quantile(V_samples, probs = c(alpha_ci/2, 1 - alpha_ci/2))

  # Probability above threshold
  prob_threshold <- mean(V_samples >= threshold)

  # Posterior mean (analytical formula for verification)
  V_mean_analytical <- sum((category_values / k) * post_alpha_vec) / sum(post_alpha_vec)

  list(
    model = "Dirichlet-Multinomial",
    V_mean = round(V_mean, 4),
    V_median = round(V_median, 4),
    V_sd = round(V_sd, 4),
    HDI = round(HDI, 4),
    ETI = round(as.numeric(ETI), 4),
    prob_threshold = round(prob_threshold, 4),
    V_mean_analytical = round(V_mean_analytical, 4),
    posterior_params = list(
      alpha = post_alpha_vec,
      prior_alpha = prior_alpha
    ),
    samples = V_samples
  )
}


# =============================================================================
# BETA-BINOMIAL MODEL
# =============================================================================

#' @noRd
.V_beta_binomial <- function(counts, k, prior_alpha, cred_level,
                             threshold, n, V_classic) {

  category_values <- 0:k

  # Total successes and trials
  sum_success <- sum(category_values * counts)
  total_trials <- n * k

  # Prior (Beta)
  prior_beta <- prior_alpha

  # Posterior parameters
  post_alpha <- prior_alpha + sum_success
  post_beta <- prior_beta + (total_trials - sum_success)

  # Posterior statistics (analytical)
  V_mean <- post_alpha / (post_alpha + post_beta)
  V_median <- stats::qbeta(0.5, post_alpha, post_beta)
  V_mode <- if (post_alpha > 1 && post_beta > 1) {
    (post_alpha - 1) / (post_alpha + post_beta - 2)
  } else NA
  V_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervals
  alpha_ci <- 1 - cred_level
  ETI <- c(stats::qbeta(alpha_ci/2, post_alpha, post_beta),
           stats::qbeta(1 - alpha_ci/2, post_alpha, post_beta))

  # HDI (numerical search)
  HDI <- .compute_hdi_beta(post_alpha, post_beta, cred_level)

  # Probability above threshold
  prob_threshold <- 1 - stats::pbeta(threshold, post_alpha, post_beta)

  list(
    model = "Beta-Binomial",
    V_mean = round(V_mean, 4),
    V_median = round(V_median, 4),
    V_mode = round(V_mode, 4),
    V_sd = round(V_sd, 4),
    HDI = round(HDI, 4),
    ETI = round(ETI, 4),
    prob_threshold = round(prob_threshold, 4),
    posterior_params = list(
      alpha = post_alpha,
      beta = post_beta,
      prior_alpha = prior_alpha
    )
  )
}


# =============================================================================
# HDI COMPUTATION
# =============================================================================

#' @noRd
.compute_hdi_beta <- function(alpha, beta, cred_level, tol = 1e-8) {
  # Find HDI for Beta distribution using optimization
  target_mass <- cred_level

  # Objective: minimize interval width
  optimize_fn <- function(lower) {
    upper <- stats::qbeta(stats::pbeta(lower, alpha, beta) + target_mass, alpha, beta)
    if (is.na(upper) || upper > 1) return(Inf)
    upper - lower
  }

  # Search range
  search_upper <- stats::qbeta(1 - target_mass, alpha, beta)
  if (is.na(search_upper)) search_upper <- 0.001

  result <- stats::optimize(optimize_fn, interval = c(0, search_upper), tol = tol)

  lower <- result$minimum
  upper <- stats::qbeta(stats::pbeta(lower, alpha, beta) + target_mass, alpha, beta)

  c(lower, upper)
}


# =============================================================================
# PRINT METHOD
# =============================================================================

#' @export
print.bayes_V <- function(x, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("         BAYESIAN AIKEN'S V - ", x$model, "\n")
  cat("=======================================================\n\n")

  cat("DATA SUMMARY:\n")
  cat(sprintf("  Number of judges:     %d\n", x$n_judges))
  cat(sprintf("  Scale range:          [%d, %d]\n", x$scale["l"], x$scale["s"]))
  cat(sprintf("  Category frequencies: %s\n", paste(x$counts, collapse = ", ")))
  cat("\n")

  cat("ESTIMATES:\n")
  cat(sprintf("  Classical V:          %.4f\n", x$V_classic))
  cat(sprintf("  Posterior mean:       %.4f\n", x$V_mean))
  cat(sprintf("  Posterior median:     %.4f\n", x$V_median))
  cat(sprintf("  Posterior SD:         %.4f\n", x$V_sd))
  cat("\n")

  cat(sprintf("CREDIBLE INTERVALS (%.0f%%):\n", x$cred_level * 100))
  cat(sprintf("  HDI: [%.4f, %.4f]  (width = %.4f)\n",
              x$HDI[1], x$HDI[2], x$HDI[2] - x$HDI[1]))
  cat(sprintf("  ETI: [%.4f, %.4f]  (width = %.4f)\n",
              x$ETI[1], x$ETI[2], x$ETI[2] - x$ETI[1]))
  cat("\n")

  cat("DECISION SUPPORT:\n")
  cat(sprintf("  P(V >= %.2f | data) = %.4f\n", x$threshold, x$prob_threshold))

  # Interpretation
  if (x$prob_threshold >= 0.95) {
    decision <- "Strong evidence for adequate validity"
  } else if (x$prob_threshold >= 0.80) {
    decision <- "Moderate evidence for adequate validity"
  } else if (x$prob_threshold >= 0.50) {
    decision <- "Weak evidence for adequate validity"
  } else {
    decision <- "Insufficient evidence for adequate validity"
  }
  cat(sprintf("  Interpretation:       %s\n", decision))

  cat("\n=======================================================\n")

  invisible(x)
}
