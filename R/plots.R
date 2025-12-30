#' =============================================================================
#' FUNCIONES DE VISUALIZACION - BayesAiken
#' =============================================================================

#' Plot del objeto bayes_aiken
#'
#' @param x Objeto de clase bayes_aiken
#' @param show_prior Mostrar la distribucion prior (default = TRUE)
#' @param show_hdi Mostrar HDI (default = TRUE)
#' @param ... Argumentos adicionales
#'
#' @export
plot.bayes_aiken <- function(x, show_prior = TRUE, show_hdi = TRUE, ...) {

  if (x$input_type == "data.frame" && !is.null(x$summary_table)) {
    # Para data.frame, mostrar comparacion de items
    plot_comparison(x)
  } else {
    # Para vector, mostrar posterior
    for (coef in names(x$results)) {
      result <- x$results[[coef]]
      plot_posterior(result$post_alpha, result$post_beta,
                     coef = coef,
                     prior_alpha = x$prior$alpha,
                     prior_beta = x$prior$beta,
                     cred_level = x$cred_level,
                     show_prior = show_prior,
                     show_hdi = show_hdi)
    }
  }

  invisible(x)
}


#' Graficar distribucion posterior
#'
#' @param post_alpha Parametro alpha de la posterior
#' @param post_beta Parametro beta de la posterior
#' @param coef Nombre del coeficiente
#' @param prior_alpha Parametro alpha del prior (default = 1)
#' @param prior_beta Parametro beta del prior (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param show_prior Mostrar prior (default = TRUE)
#' @param show_hdi Mostrar HDI (default = TRUE)
#' @param main Titulo del grafico (opcional)
#'
#' @return Invisible NULL
#'
#' @examples
#' plot_posterior(post_alpha = 25, post_beta = 5, coef = "V")
#'
#' @export
plot_posterior <- function(post_alpha, post_beta, coef = "V",
                           prior_alpha = 1, prior_beta = 1,
                           cred_level = 0.95,
                           show_prior = TRUE, show_hdi = TRUE,
                           main = NULL) {

  # Secuencia de valores
  x <- seq(0, 1, length.out = 500)

  # Densidades
  posterior <- dbeta(x, post_alpha, post_beta)

  # Calcular HDI
  hdi <- compute_hdi(post_alpha, post_beta, cred_level)

  # Media posterior
  mean_post <- post_alpha / (post_alpha + post_beta)

  # Configurar titulo
  if (is.null(main)) {
    coef_names <- c(
      V = "Aiken's V (Content Validity)",
      H = "H (Homogeneity)",
      R = "R (Reproducibility)",
      C = "C (Consensus)",
      A = "A (Judge Agreement)",
      I = "I (Inter-item Consistency)"
    )
    main <- paste("Posterior Distribution -", coef_names[coef])
  }

  # Configurar margenes
  old_par <- par(mar = c(5, 4, 4, 2) + 0.1)
  on.exit(par(old_par))

  # Plot base
  y_max <- max(posterior) * 1.15

  plot(x, posterior, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = paste("Coefficient", coef),
       ylab = "Density",
       main = main,
       las = 1)

  # Prior (si se solicita)
  if (show_prior) {
    prior <- dbeta(x, prior_alpha, prior_beta)
    lines(x, prior, col = "gray60", lwd = 2, lty = 2)
  }

  # HDI shading
  if (show_hdi) {
    x_hdi <- x[x >= hdi[1] & x <= hdi[2]]
    y_hdi <- dbeta(x_hdi, post_alpha, post_beta)
    polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
            c(0, y_hdi, 0),
            col = rgb(0.2, 0.5, 0.8, 0.3),
            border = NA)
  }

  # Posterior
  lines(x, posterior, col = "steelblue", lwd = 3)

  # Lineas de referencia
  abline(v = 0.70, col = "orange", lty = 2, lwd = 1.5)
  abline(v = 0.80, col = "red", lty = 2, lwd = 1.5)

  # Media posterior
  abline(v = mean_post, col = "steelblue", lty = 3, lwd = 2)

  # HDI lines
  if (show_hdi) {
    abline(v = hdi[1], col = "steelblue", lty = 2, lwd = 1.5)
    abline(v = hdi[2], col = "steelblue", lty = 2, lwd = 1.5)
  }

  # Leyenda
  legend_items <- c(
    paste0("Posterior: Beta(", round(post_alpha, 1), ", ", round(post_beta, 1), ")"),
    paste0("Mean: ", round(mean_post, 3)),
    paste0(round(cred_level * 100), "% HDI: [", round(hdi[1], 3), ", ", round(hdi[2], 3), "]"),
    "Threshold 0.70",
    "Threshold 0.80"
  )
  legend_cols <- c("steelblue", "steelblue", "steelblue", "orange", "red")
  legend_lty <- c(1, 3, 2, 2, 2)

  if (show_prior) {
    legend_items <- c(paste0("Prior: Beta(", prior_alpha, ", ", prior_beta, ")"), legend_items)
    legend_cols <- c("gray60", legend_cols)
    legend_lty <- c(2, legend_lty)
  }

  legend("topright",
         legend = legend_items,
         col = legend_cols,
         lty = legend_lty,
         lwd = c(2, rep(2, length(legend_items) - 1)),
         bty = "n",
         cex = 0.8)

  # Probabilidades en el grafico
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  mtext(sprintf("P(%s > 0.70) = %.3f  |  P(%s > 0.80) = %.3f",
                coef, prob_70, coef, prob_80),
        side = 1, line = 3.5, cex = 0.9)

  invisible(NULL)
}


#' Graficar todos los coeficientes
#'
#' @param results Lista de resultados de coeficientes
#' @param show_prior Mostrar prior (default = FALSE)
#'
#' @export
plot_all_coefficients <- function(results, show_prior = FALSE) {

  n_coefs <- length(results$results)

  if (n_coefs == 0) {
    message("No hay resultados para graficar")
    return(invisible(NULL))
  }

  # Configurar layout
  if (n_coefs > 1) {
    n_cols <- min(2, n_coefs)
    n_rows <- ceiling(n_coefs / n_cols)
    old_par <- par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    on.exit(par(old_par))
  }

  # Graficar cada coeficiente
  for (coef in names(results$results)) {
    result <- results$results[[coef]]
    plot_posterior(result$post_alpha, result$post_beta,
                   coef = coef,
                   prior_alpha = results$prior$alpha,
                   prior_beta = results$prior$beta,
                   cred_level = results$cred_level,
                   show_prior = show_prior,
                   show_hdi = TRUE)
  }

  invisible(NULL)
}


#' Grafico de comparacion de items
#'
#' @param x Objeto bayes_aiken con multiples items
#' @param coef Coeficiente a graficar (si hay varios)
#' @param order_by Ordenar por: "none", "value", "prob" (default = "none")
#'
#' @export
plot_comparison <- function(x, coef = NULL, order_by = "none") {

  # Obtener datos
  if (is.data.frame(x$summary_table)) {
    data <- x$summary_table
  } else if (!is.null(x$results) && is.data.frame(x$results[[1]])) {
    data <- x$results[[1]]
    if (is.null(coef)) coef <- names(x$results)[1]
  } else {
    message("No hay datos de multiples items para comparar")
    return(invisible(NULL))
  }

  # Determinar columnas
  if ("V_bayesiano" %in% names(data)) {
    value_col <- "V_bayesiano"
    ci_lower <- "IC_lower"
    ci_upper <- "IC_upper"
    prob_col <- "P_70"
  } else if ("bayesiano" %in% names(data)) {
    value_col <- "bayesiano"
    ci_lower <- "IC_lower"
    ci_upper <- "IC_upper"
    prob_col <- "P_70"
  } else {
    message("Formato de datos no reconocido")
    return(invisible(NULL))
  }

  n_items <- nrow(data)

  # Ordenar si se solicita
  if (order_by == "value") {
    data <- data[order(data[[value_col]], decreasing = TRUE), ]
  } else if (order_by == "prob") {
    data <- data[order(data[[prob_col]], decreasing = TRUE), ]
  }

  # Configurar grafico
  old_par <- par(mar = c(5, 8, 4, 2))
  on.exit(par(old_par))

  # Colores segun probabilidad > 0.70
  colors <- ifelse(data[[prob_col]] > 0.95, "forestgreen",
                   ifelse(data[[prob_col]] > 0.80, "orange", "red"))

  # Posiciones
  y_pos <- n_items:1

  # Labels
  if ("item" %in% names(data)) {
    labels <- as.character(data$item)
  } else {
    labels <- paste("Item", 1:n_items)
  }

  # Grafico base
  plot(data[[value_col]], y_pos,
       xlim = c(0, 1),
       ylim = c(0.5, n_items + 0.5),
       xlab = "Coefficient Value",
       ylab = "",
       yaxt = "n",
       main = paste("Bayesian Content Validity -", ifelse(!is.null(coef), coef, "V")),
       pch = 19, col = colors, cex = 1.5)

  # Eje Y
  axis(2, at = y_pos, labels = labels, las = 1, cex.axis = 0.8)

  # Barras de error (IC)
  segments(data[[ci_lower]], y_pos, data[[ci_upper]], y_pos,
           col = colors, lwd = 2)

  # Lineas de referencia
  abline(v = 0.70, col = "orange", lty = 2, lwd = 1.5)
  abline(v = 0.80, col = "red", lty = 2, lwd = 1.5)

  # Grid horizontal
  abline(h = y_pos, col = "gray90", lty = 1)

  # Leyenda
  legend("bottomright",
         legend = c("P(>0.70) > 95%", "P(>0.70) 80-95%", "P(>0.70) < 80%",
                    "Threshold 0.70", "Threshold 0.80"),
         col = c("forestgreen", "orange", "red", "orange", "red"),
         pch = c(19, 19, 19, NA, NA),
         lty = c(NA, NA, NA, 2, 2),
         lwd = c(NA, NA, NA, 1.5, 1.5),
         bty = "n",
         cex = 0.8)

  invisible(NULL)
}


#' Grafico interno para un solo coeficiente
#' @noRd
.plot_single <- function(result, coef, cred_level, prior_alpha, prior_beta) {
  plot_posterior(result$post_alpha, result$post_beta,
                 coef = coef,
                 prior_alpha = prior_alpha,
                 prior_beta = prior_beta,
                 cred_level = cred_level,
                 show_prior = TRUE,
                 show_hdi = TRUE)
}


#' Generar plot y retornarlo como objeto
#' @noRd
.generate_plot <- function(result, coef, cred_level) {
  # Retorna los datos necesarios para replotear
  list(
    post_alpha = result$post_alpha,
    post_beta = result$post_beta,
    coef = coef,
    cred_level = cred_level
  )
}
