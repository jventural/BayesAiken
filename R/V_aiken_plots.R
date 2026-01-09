#' =============================================================================
#' VISUALIZATION FUNCTIONS - BayesAiken Package
#' =============================================================================
#'
#' Plot functions for Bayesian Aiken's V analysis.
#'
#' =============================================================================

# =============================================================================
# PLOT METHOD FOR bayes_V
# =============================================================================

#' Plot Bayesian V Posterior Distribution
#'
#' @param x Object of class 'bayes_V' from V_aiken()
#' @param show_prior Show prior distribution (default = TRUE for Beta, FALSE for Dirichlet)
#' @param show_hdi Show HDI shaded region (default = TRUE)
#' @param ... Additional arguments (ignored)
#'
#' @export
plot.bayes_V <- function(x, show_prior = NULL, show_hdi = TRUE, ...) {

  if (x$model == "Dirichlet-Multinomial") {
    if (is.null(show_prior)) show_prior <- FALSE
    .plot_V_dirichlet_internal(x, show_prior, show_hdi)
  } else {
    if (is.null(show_prior)) show_prior <- TRUE
    .plot_V_beta_internal(x, show_prior, show_hdi)
  }

  invisible(x)
}


# =============================================================================
# COMPARISON PLOT: Both Models
# =============================================================================

#' Compare Beta-Binomial and Dirichlet-Multinomial Posteriors
#'
#' Creates a comparison plot showing posterior distributions from both models,
#' similar to Figure 1 in Ventura-León (2025).
#'
#' @param ratings Vector of judge ratings
#' @param l Minimum scale value (default = 0)
#' @param s Maximum scale value (default = 3)
#' @param prior_alpha Prior concentration parameter (default = 1)
#' @param cred_level Credibility level (default = 0.95)
#' @param n_samples Number of Monte Carlo samples (default = 10000)
#' @param threshold Decision threshold (default = 0.70)
#' @param seed Random seed (default = 2025)
#' @param main Plot title (optional)
#' @param col_beta Color for Beta-Binomial (default = "steelblue")
#' @param col_dir Color for Dirichlet (default = "coral")
#'
#' @return Invisible list with both model results
#'
#' @examples
#' ratings <- c(2, 3, 3, 3, 3)
#' plot_V_comparison(ratings, l = 0, s = 3)
#'
#' @export
plot_V_comparison <- function(ratings,
                               l = 0,
                               s = 3,
                               prior_alpha = 1,
                               cred_level = 0.95,
                               n_samples = 10000,
                               threshold = 0.70,
                               seed = 2025,
                               main = NULL,
                               col_beta = "steelblue",
                               col_dir = "coral") {

  # Calculate both models
  result_beta <- V_aiken(ratings, l = l, s = s, model = "beta",
                         prior_alpha = prior_alpha, cred_level = cred_level,
                         threshold = threshold, verbose = FALSE)

  result_dir <- V_aiken(ratings, l = l, s = s, model = "dirichlet",
                        prior_alpha = prior_alpha, cred_level = cred_level,
                        n_samples = n_samples, threshold = threshold,
                        seed = seed, verbose = FALSE)

  # Title
  if (is.null(main)) {
    main <- sprintf("Posterior Comparison (n = %d judges)", result_beta$n_judges)
  }

  # Beta-Binomial density
  x_seq <- seq(0.001, 0.999, length.out = 500)
  y_beta <- stats::dbeta(x_seq,
                         result_beta$posterior_params$alpha,
                         result_beta$posterior_params$beta)

  # Dirichlet density (kernel estimate from samples)
  dens_dir <- stats::density(result_dir$samples, from = 0, to = 1, n = 512)

  # Setup plot
  y_max <- max(c(y_beta, dens_dir$y)) * 1.15
  old_par <- graphics::par(mar = c(5, 5, 4, 2))
  on.exit(graphics::par(old_par))

  # Base plot
  plot(x_seq, y_beta, type = "l", col = col_beta, lwd = 3,
       xlim = c(0.3, 1), ylim = c(0, y_max),
       xlab = "Aiken's V", ylab = "Posterior Density",
       main = main, cex.lab = 1.2, cex.axis = 1.1)

  # Grid
  graphics::grid(col = "gray90", lty = 1)

  # Redraw Beta line
  graphics::lines(x_seq, y_beta, col = col_beta, lwd = 3)

  # Dirichlet density
  graphics::lines(dens_dir$x, dens_dir$y, col = col_dir, lwd = 3)

  # Reference lines
  graphics::abline(v = result_beta$V_classic, lty = 2, col = "gray40", lwd = 2)
  graphics::abline(v = threshold, lty = 3, col = "gray60", lwd = 1.5)

  # Legend
  graphics::legend("topleft",
                   legend = c("Beta-Binomial", "Dirichlet-Multinomial",
                              bquote(hat(V) == .(sprintf("%.2f", result_beta$V_classic))),
                              sprintf("Threshold = %.2f", threshold)),
                   col = c(col_beta, col_dir, "gray40", "gray60"),
                   lty = c(1, 1, 2, 3),
                   lwd = c(3, 3, 2, 1.5),
                   bty = "n", cex = 0.9)

  # Statistics annotations (bottom right, where curves are low)
  y_pos_beta <- y_max * 0.38
  y_pos_dir <- y_max * 0.18

  # Beta-Binomial stats
  .text_box(0.99, y_pos_beta,
            sprintf("Beta-Binomial:\n95%% HDI = [%.2f, %.2f]\nWidth = %.2f",
                    result_beta$HDI[1], result_beta$HDI[2],
                    result_beta$HDI[2] - result_beta$HDI[1]),
            col = col_beta, adj = 1)

  # Dirichlet stats
  .text_box(0.99, y_pos_dir,
            sprintf("Dirichlet:\n95%% HDI = [%.2f, %.2f]\nWidth = %.2f",
                    result_dir$HDI[1], result_dir$HDI[2],
                    result_dir$HDI[2] - result_dir$HDI[1]),
            col = col_dir, adj = 1)

  invisible(list(beta = result_beta, dirichlet = result_dir))
}


# =============================================================================
# SINGLE MODEL PLOTS (Internal)
# =============================================================================

#' @noRd
.plot_V_beta_internal <- function(result, show_prior, show_hdi) {

  post_alpha <- result$posterior_params$alpha
  post_beta <- result$posterior_params$beta
  prior_alpha <- result$posterior_params$prior_alpha

  x_seq <- seq(0.001, 0.999, length.out = 500)
  y_post <- stats::dbeta(x_seq, post_alpha, post_beta)

  y_max <- max(y_post) * 1.15
  old_par <- graphics::par(mar = c(5, 5, 4, 2))
  on.exit(graphics::par(old_par))

  # Base plot
  plot(x_seq, y_post, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "Aiken's V", ylab = "Posterior Density",
       main = sprintf("Beta-Binomial Posterior (n = %d judges)", result$n_judges),
       cex.lab = 1.1, cex.axis = 1)

  graphics::grid(col = "gray90", lty = 1)

  # Prior
  if (show_prior) {
    y_prior <- stats::dbeta(x_seq, prior_alpha, prior_alpha)
    if (max(y_prior) > 0 && is.finite(max(y_prior))) {
      y_prior_scaled <- y_prior * (y_max * 0.25) / max(y_prior)
      graphics::lines(x_seq, y_prior_scaled, col = "gray60", lwd = 2, lty = 2)
    }
  }

  # HDI shading
  if (show_hdi) {
    x_hdi <- x_seq[x_seq >= result$HDI[1] & x_seq <= result$HDI[2]]
    y_hdi <- stats::dbeta(x_hdi, post_alpha, post_beta)
    if (length(x_hdi) > 0) {
      graphics::polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
                        c(0, y_hdi, 0),
                        col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
                        border = NA)
    }
  }

  # Posterior
  graphics::lines(x_seq, y_post, col = "steelblue", lwd = 3)

  # Reference lines
  graphics::abline(v = result$threshold, col = "darkorange", lty = 2, lwd = 2)
  graphics::abline(v = result$V_mean, col = "steelblue", lty = 3, lwd = 2)

  # Stats text
  graphics::text(0.02, y_max * 0.95,
                 sprintf("V classic = %.3f", result$V_classic),
                 adj = c(0, 1), cex = 0.9, col = "gray40")
  graphics::text(0.02, y_max * 0.87,
                 sprintf("V posterior = %.3f", result$V_mean),
                 adj = c(0, 1), cex = 0.9, font = 2, col = "steelblue")
  graphics::text(0.02, y_max * 0.79,
                 sprintf("95%% HDI: [%.3f, %.3f]", result$HDI[1], result$HDI[2]),
                 adj = c(0, 1), cex = 0.85, col = "gray30")
  graphics::text(0.02, y_max * 0.71,
                 sprintf("P(V > %.2f) = %.3f", result$threshold, result$prob_threshold),
                 adj = c(0, 1), cex = 0.9, font = 2, col = "darkorange")

  # Legend
  legend_items <- c(
    sprintf("Posterior Beta(%.1f, %.1f)", post_alpha, post_beta),
    "95% HDI",
    sprintf("Threshold %.2f", result$threshold)
  )
  legend_cols <- c("steelblue", grDevices::adjustcolor("steelblue", 0.5), "darkorange")
  legend_lty <- c(1, NA, 2)
  legend_pch <- c(NA, 15, NA)

  if (show_prior) {
    legend_items <- c(legend_items, sprintf("Prior Beta(%.0f, %.0f)", prior_alpha, prior_alpha))
    legend_cols <- c(legend_cols, "gray60")
    legend_lty <- c(legend_lty, 2)
    legend_pch <- c(legend_pch, NA)
  }

  graphics::legend("topright",
                   legend = legend_items,
                   col = legend_cols,
                   lty = legend_lty,
                   lwd = c(3, NA, 2, if(show_prior) 2 else NULL),
                   pch = legend_pch,
                   pt.cex = c(NA, 2, NA, NA),
                   bty = "n", cex = 0.8)
}


#' @noRd
.plot_V_dirichlet_internal <- function(result, show_prior, show_hdi) {

  samples <- result$samples
  dens <- stats::density(samples, from = 0, to = 1, n = 512)

  y_max <- max(dens$y) * 1.15
  old_par <- graphics::par(mar = c(5, 5, 4, 2))
  on.exit(graphics::par(old_par))

  # Base plot
  plot(dens$x, dens$y, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "Aiken's V", ylab = "Posterior Density",
       main = sprintf("Dirichlet-Multinomial Posterior (n = %d judges)", result$n_judges),
       cex.lab = 1.1, cex.axis = 1)

  graphics::grid(col = "gray90", lty = 1)

  # HDI shading
  if (show_hdi) {
    x_hdi <- dens$x[dens$x >= result$HDI[1] & dens$x <= result$HDI[2]]
    y_hdi <- dens$y[dens$x >= result$HDI[1] & dens$x <= result$HDI[2]]
    if (length(x_hdi) > 0) {
      graphics::polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
                        c(0, y_hdi, 0),
                        col = grDevices::adjustcolor("coral", alpha.f = 0.3),
                        border = NA)
    }
  }

  # Posterior
  graphics::lines(dens$x, dens$y, col = "coral", lwd = 3)

  # Reference lines
  graphics::abline(v = result$threshold, col = "darkorange", lty = 2, lwd = 2)
  graphics::abline(v = result$V_mean, col = "coral", lty = 3, lwd = 2)

  # Stats text
  graphics::text(0.02, y_max * 0.95,
                 sprintf("V classic = %.3f", result$V_classic),
                 adj = c(0, 1), cex = 0.9, col = "gray40")
  graphics::text(0.02, y_max * 0.87,
                 sprintf("V posterior = %.3f", result$V_mean),
                 adj = c(0, 1), cex = 0.9, font = 2, col = "coral")
  graphics::text(0.02, y_max * 0.79,
                 sprintf("95%% HDI: [%.3f, %.3f]", result$HDI[1], result$HDI[2]),
                 adj = c(0, 1), cex = 0.85, col = "gray30")
  graphics::text(0.02, y_max * 0.71,
                 sprintf("P(V > %.2f) = %.3f", result$threshold, result$prob_threshold),
                 adj = c(0, 1), cex = 0.9, font = 2, col = "darkorange")

  # Legend
  graphics::legend("topright",
                   legend = c(
                     sprintf("Posterior Dirichlet(%s)",
                             paste(result$posterior_params$alpha, collapse = ",")),
                     "95% HDI",
                     sprintf("Threshold %.2f", result$threshold)
                   ),
                   col = c("coral", grDevices::adjustcolor("coral", 0.5), "darkorange"),
                   lty = c(1, NA, 2),
                   lwd = c(3, NA, 2),
                   pch = c(NA, 15, NA),
                   pt.cex = c(NA, 2, NA),
                   bty = "n", cex = 0.8)
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' @noRd
.text_box <- function(x, y, label, col, adj = 0, cex = 0.85) {
  # Draw text with white background box
  w <- graphics::strwidth(label, cex = cex) * 1.1
  h <- graphics::strheight(label, cex = cex) * 1.3

  if (adj == 1) {
    x_left <- x - w
    x_right <- x + w * 0.05
  } else {
    x_left <- x - w * 0.05
    x_right <- x + w
  }

  graphics::rect(x_left, y - h * 0.6, x_right, y + h * 0.5,
                 col = "white", border = NA)
  graphics::text(x, y, label, cex = cex, col = col, adj = adj)
}


# =============================================================================
# FOREST PLOT FOR MULTIPLE ITEMS
# =============================================================================

#' Forest Plot for Multiple Items
#'
#' Creates a forest plot displaying Bayesian V estimates with credible
#' intervals for multiple items.
#'
#' @param results_list Named list of bayes_V objects from V_aiken()
#' @param order_by How to order items: "none", "value", "prob" (default = "none")
#' @param main Plot title (optional)
#'
#' @examples
#' # Analyze multiple items
#' items <- list(
#'   "Item 1" = V_aiken(c(3,3,3,3,2), verbose = FALSE),
#'   "Item 2" = V_aiken(c(3,3,2,2,2), verbose = FALSE),
#'   "Item 3" = V_aiken(c(3,3,3,2,3), verbose = FALSE)
#' )
#' plot_V_forest(items)
#'
#' @export
plot_V_forest <- function(results_list,
                           order_by = c("none", "value", "prob"),
                           main = "Bayesian Content Validity Analysis") {

  order_by <- match.arg(order_by)
  n_items <- length(results_list)

  # Extract data
  labels <- names(results_list)
  if (is.null(labels)) labels <- paste("Item", 1:n_items)

  values <- sapply(results_list, function(x) x$V_mean)
  lower <- sapply(results_list, function(x) x$HDI[1])
  upper <- sapply(results_list, function(x) x$HDI[2])
  probs <- sapply(results_list, function(x) x$prob_threshold)
  threshold <- results_list[[1]]$threshold

  # Order if requested
  if (order_by == "value") {
    ord <- order(values, decreasing = TRUE)
  } else if (order_by == "prob") {
    ord <- order(probs, decreasing = TRUE)
  } else {
    ord <- 1:n_items
  }

  labels <- labels[ord]
  values <- values[ord]
  lower <- lower[ord]
  upper <- upper[ord]
  probs <- probs[ord]

  # Colors based on probability
  colors <- ifelse(probs >= 0.95, "forestgreen",
                   ifelse(probs >= 0.80, "darkorange", "firebrick"))

  # Setup plot
  old_par <- graphics::par(mar = c(5, 10, 4, 8))
  on.exit(graphics::par(old_par))

  y_pos <- n_items:1

  plot(values, y_pos,
       xlim = c(0, 1),
       ylim = c(0.5, n_items + 0.5),
       xlab = "Aiken's V",
       ylab = "",
       yaxt = "n",
       main = main,
       pch = 19, col = colors, cex = 1.5)

  graphics::axis(2, at = y_pos, labels = labels, las = 1, cex.axis = 0.8)
  graphics::abline(h = y_pos, col = "gray90", lty = 1)
  graphics::segments(lower, y_pos, upper, y_pos, col = colors, lwd = 2)
  graphics::abline(v = threshold, col = "darkorange", lty = 2, lwd = 2)

  # Legend
  graphics::legend("bottomright",
                   legend = c(sprintf("P(V>%.2f) >= 95%%", threshold),
                              sprintf("P(V>%.2f) 80-95%%", threshold),
                              sprintf("P(V>%.2f) < 80%%", threshold),
                              sprintf("Threshold = %.2f", threshold)),
                   col = c("forestgreen", "darkorange", "firebrick", "darkorange"),
                   pch = c(19, 19, 19, NA),
                   lty = c(NA, NA, NA, 2),
                   lwd = c(NA, NA, NA, 2),
                   bty = "n", cex = 0.8)

  invisible(NULL)
}
