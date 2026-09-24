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
#' @param compact Modo compacto para paneles multiples (default = FALSE)
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
                           main = NULL, compact = FALSE) {

  # Secuencia de valores
  x <- seq(0.001, 0.999, length.out = 500)

  # Densidades
  posterior <- dbeta(x, post_alpha, post_beta)

  # Calcular HDI
  hdi <- compute_hdi(post_alpha, post_beta, cred_level)

  # Media posterior
  mean_post <- post_alpha / (post_alpha + post_beta)

  # Probabilidades

  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

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

  # Configurar margenes segun modo
  if (compact) {
    old_par <- par(mar = c(4, 4, 3, 1))
  } else {
    old_par <- par(mar = c(6, 4.5, 4, 8), xpd = TRUE)
  }
  on.exit(par(old_par))

  # Plot base
  y_max <- max(posterior) * 1.1

  plot(x, posterior, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "",
       ylab = "Density",
       main = main,
       las = 1,
       cex.main = ifelse(compact, 0.9, 1.1),
       cex.lab = ifelse(compact, 0.8, 1),
       cex.axis = ifelse(compact, 0.7, 0.9))

  # Etiqueta del eje X con espacio adecuado
  mtext(paste("Coefficient", coef), side = 1, line = 2.5,
        cex = ifelse(compact, 0.8, 1))

  # Grid de fondo
  grid(col = "gray90", lty = 1, lwd = 0.5)

  # Prior (si se solicita)
  if (show_prior && !compact) {
    prior <- dbeta(x, prior_alpha, prior_beta)
    # Escalar el prior si es muy diferente
    if (max(prior) > 0 && is.finite(max(prior))) {
      prior_scaled <- prior * (y_max * 0.3) / max(prior)
      lines(x, prior_scaled, col = "gray50", lwd = 2, lty = 2)
    }
  }

  # HDI shading
  if (show_hdi) {
    x_hdi <- x[x >= hdi[1] & x <= hdi[2]]
    y_hdi <- dbeta(x_hdi, post_alpha, post_beta)
    if (length(x_hdi) > 0) {
      polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
              c(0, y_hdi, 0),
              col = rgb(0.27, 0.51, 0.71, 0.35),
              border = NA)
    }
  }

  # Posterior
  lines(x, posterior, col = "steelblue", lwd = 3)

  # Linea de referencia (threshold)
  abline(v = 0.70, col = "darkorange", lty = 2, lwd = 2)

  # Media posterior
  abline(v = mean_post, col = "steelblue", lty = 3, lwd = 2)

  # HDI lines
  if (show_hdi) {
    segments(hdi[1], 0, hdi[1], dbeta(hdi[1], post_alpha, post_beta),
             col = "steelblue", lty = 2, lwd = 1.5)
    segments(hdi[2], 0, hdi[2], dbeta(hdi[2], post_alpha, post_beta),
             col = "steelblue", lty = 2, lwd = 1.5)
  }

  # Texto de estadisticos en el grafico (parte superior izquierda)
  if (!compact) {
    # Caja de texto con estadisticos
    text_x <- 0.02
    text_y <- y_max * 0.95

    text(text_x, text_y,
         sprintf("%s = %.3f", coef, mean_post),
         adj = c(0, 1), cex = 0.9, font = 2, col = "steelblue")
    text(text_x, text_y - y_max * 0.08,
         sprintf("%.0f%% HDI: [%.3f, %.3f]", cred_level * 100, hdi[1], hdi[2]),
         adj = c(0, 1), cex = 0.85, col = "gray30")
    text(text_x, text_y - y_max * 0.16,
         sprintf("P(%s > 0.70) = %.3f", coef, prob_70),
         adj = c(0, 1), cex = 0.85, col = "darkorange")

    # Leyenda FUERA del grafico (a la derecha)
    legend_x <- 1.02
    legend_y <- y_max * 0.95

    legend(legend_x, legend_y,
           legend = c(
             paste0("Posterior Beta(", round(post_alpha, 1), ", ", round(post_beta, 1), ")"),
             paste0("Mean = ", round(mean_post, 3)),
             paste0(round(cred_level * 100), "% HDI"),
             "Threshold 0.70"
           ),
           col = c("steelblue", "steelblue", rgb(0.27, 0.51, 0.71, 0.5), "darkorange"),
           lty = c(1, 3, NA, 2),
           lwd = c(3, 2, NA, 2),
           pch = c(NA, NA, 15, NA),
           pt.cex = c(NA, NA, 2, NA),
           bty = "n",
           cex = 0.75,
           xpd = TRUE)

    # Si se muestra prior, agregar a la leyenda
    if (show_prior) {
      legend(legend_x, legend_y - y_max * 0.45,
             legend = paste0("Prior Beta(", prior_alpha, ", ", prior_beta, ")"),
             col = "gray50",
             lty = 2,
             lwd = 2,
             bty = "n",
             cex = 0.75,
             xpd = TRUE)
    }
  } else {
    # Modo compacto: solo mostrar valor y HDI dentro del grafico
    text(0.05, y_max * 0.9,
         sprintf("%s=%.2f [%.2f,%.2f]", coef, mean_post, hdi[1], hdi[2]),
         adj = c(0, 1), cex = 0.7, col = "steelblue")
  }

  invisible(NULL)
}


#' Graficar multiples items en un panel
#'
#' Genera un panel con multiples graficos de distribucion posterior
#'
#' @param results_list Lista de resultados o objeto bayes_aiken con multiples items
#' @param n_cols Numero de columnas en el panel (default = 4)
#' @param show_prior Mostrar prior (default = FALSE)
#' @param main_title Titulo principal del panel (opcional
#' @param coef Tipo de coeficiente (default = "V")
#' @param prior_alpha Alpha del prior (default = 1)
#' @param prior_beta Beta del prior (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#'
#' @return Invisible NULL
#'
#' @examples
#' # Crear datos de ejemplo
#' items_data <- lapply(1:20, function(i) {
#'   list(post_alpha = 20 + sample(-5:10, 1),
#'        post_beta = 3 + sample(0:3, 1),
#'        item_name = paste0("Item ", i))
#' })
#' plot_multiple_posteriors(items_data, n_cols = 5)
#'
#' @export
plot_multiple_posteriors <- function(results_list,
                                     n_cols = 4,
                                     show_prior = FALSE,
                                     main_title = "Posterior Distributions by Item",
                                     coef = "V",
                                     prior_alpha = 1,
                                     prior_beta = 1,
                                     cred_level = 0.95) {

  # Si es un objeto bayes_aiken, extraer los resultados
 if (inherits(results_list, "bayes_aiken")) {
    if (!is.null(results_list$results) && is.data.frame(results_list$results[[1]])) {
      df <- results_list$results[[1]]
      # Reconstruir lista de resultados desde el data.frame
      items <- unique(df$item)
      results_list <- lapply(seq_along(items), function(i) {
        row <- df[i, ]
        # Recalcular post_alpha y post_beta desde los datos
        # Aproximacion basada en media y varianza
        mean_val <- row$V_bayesiano
        n_judges <- row$n_jueces
        # Estimacion aproximada
        list(
          post_alpha = mean_val * n_judges * 3 + 1,
          post_beta = (1 - mean_val) * n_judges * 3 + 1,
          item_name = as.character(row$item)
        )
      })
    }
  }

  n_items <- length(results_list)

  if (n_items == 0) {
    message("No hay items para graficar")
    return(invisible(NULL))
  }

  # Calcular filas necesarias
  n_rows <- ceiling(n_items / n_cols)

  # Configurar layout
  old_par <- par(mfrow = c(n_rows, n_cols),
                 mar = c(3, 3, 2.5, 1),
                 oma = c(2, 2, 3, 1),
                 mgp = c(2, 0.7, 0))
  on.exit(par(old_par))

  # Graficar cada item
  for (i in seq_along(results_list)) {
    item <- results_list[[i]]

    # Obtener nombre del item
    item_name <- item$item_name
    if (is.null(item_name)) {
      item_name <- paste("Item", i)
    }

    # Obtener parametros
    post_alpha <- item$post_alpha
    post_beta <- item$post_beta

    if (is.null(post_alpha) || is.null(post_beta)) {
      next
    }

    # Graficar
    .plot_posterior_compact(post_alpha, post_beta,
                            item_name = item_name,
                            coef = coef,
                            prior_alpha = prior_alpha,
                            prior_beta = prior_beta,
                            cred_level = cred_level,
                            show_prior = show_prior)
  }

  # Titulo principal
  mtext(main_title, outer = TRUE, line = 1, cex = 1.2, font = 2)

  # Leyenda global en el margen inferior
  mtext("Orange line = 0.70 threshold | Blue area = HDI",
        outer = TRUE, side = 1, line = 0.5, cex = 0.8, col = "gray40")

  invisible(NULL)
}


#' Grafico compacto para paneles
#' @noRd
.plot_posterior_compact <- function(post_alpha, post_beta,
                                    item_name = "Item",
                                    coef = "V",
                                    prior_alpha = 1, prior_beta = 1,
                                    cred_level = 0.95,
                                    show_prior = FALSE) {

  # Secuencia de valores
  x <- seq(0.001, 0.999, length.out = 300)

  # Densidades
  posterior <- dbeta(x, post_alpha, post_beta)

  # Calcular HDI
  hdi <- compute_hdi(post_alpha, post_beta, cred_level)

  # Media posterior
  mean_post <- post_alpha / (post_alpha + post_beta)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)

  # Color del titulo segun probabilidad
  title_col <- if (prob_70 >= 0.95) "forestgreen" else if (prob_70 >= 0.80) "darkorange" else "firebrick"

  # Plot base
  y_max <- max(posterior) * 1.1

  plot(x, posterior, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "", ylab = "",
       main = item_name,
       col.main = title_col,
       cex.main = 0.9,
       font.main = 2,
       las = 1,
       cex.axis = 0.7,
       xaxt = "n")

  # Eje X simplificado
  axis(1, at = c(0, 0.5, 0.7, 1), labels = c("0", ".5", ".7", "1"), cex.axis = 0.7)

  # Grid
  abline(h = seq(0, y_max, length.out = 5), col = "gray95", lty = 1)
  abline(v = c(0.5, 0.7), col = "gray90", lty = 1)

  # HDI shading
  x_hdi <- x[x >= hdi[1] & x <= hdi[2]]
  y_hdi <- dbeta(x_hdi, post_alpha, post_beta)
  if (length(x_hdi) > 0) {
    polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
            c(0, y_hdi, 0),
            col = rgb(0.27, 0.51, 0.71, 0.4),
            border = NA)
  }

  # Posterior
  lines(x, posterior, col = "steelblue", lwd = 2)

  # Linea de referencia
  abline(v = 0.70, col = "darkorange", lty = 2, lwd = 1.5)

  # Media
  abline(v = mean_post, col = "steelblue", lty = 3, lwd = 1.5)

  # Texto con valor
  text(0.03, y_max * 0.92,
       sprintf("%s=%.2f", coef, mean_post),
       adj = c(0, 1), cex = 0.65, font = 2, col = "steelblue")

  text(0.03, y_max * 0.75,
       sprintf("P>.70=%.2f", prob_70),
       adj = c(0, 1), cex = 0.6, col = title_col)
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
    old_par <- par(mfrow = c(n_rows, n_cols), mar = c(5, 4, 3, 6), oma = c(0, 0, 2, 0))
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
                   show_hdi = TRUE,
                   compact = (n_coefs > 1))
  }

  if (n_coefs > 1) {
    mtext("Bayesian Content Validity Analysis", outer = TRUE, line = 0.5, cex = 1.2, font = 2)
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

  # Configurar grafico con margen derecho para leyenda
  old_par <- par(mar = c(5, 10, 4, 10), xpd = TRUE)
  on.exit(par(old_par))

  # Colores segun probabilidad > 0.70
  colors <- ifelse(data[[prob_col]] > 0.95, "forestgreen",
                   ifelse(data[[prob_col]] > 0.80, "darkorange", "firebrick"))

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

  # Eje Y con labels
  axis(2, at = y_pos, labels = labels, las = 1, cex.axis = 0.8)

  # Grid horizontal
  abline(h = y_pos, col = "gray90", lty = 1)

  # Barras de error (IC)
  segments(data[[ci_lower]], y_pos, data[[ci_upper]], y_pos,
           col = colors, lwd = 2)

  # Linea de referencia
  abline(v = 0.70, col = "darkorange", lty = 2, lwd = 2)

  # Leyenda FUERA del grafico (a la derecha)
  legend(1.02, n_items * 0.8,
         legend = c("P(>0.70) > 95%", "P(>0.70) 80-95%", "P(>0.70) < 80%",
                    "Threshold 0.70"),
         col = c("forestgreen", "darkorange", "firebrick", "darkorange"),
         pch = c(19, 19, 19, NA),
         lty = c(NA, NA, NA, 2),
         lwd = c(NA, NA, NA, 2),
         bty = "n",
         cex = 0.8,
         xpd = TRUE)

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


# =============================================================================
# FUNCIONES DE PLOT SEPARADAS PARA CADA MODELO
# =============================================================================

#' Graficar densidad posterior del modelo Beta-Binomial
#'
#' Genera un grafico de la distribucion posterior de V usando el modelo
#' Beta-Binomial estandar.
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param threshold Umbral de referencia (default = 0.70)
#' @param main Titulo del grafico (opcional)
#' @param col Color de la densidad (default = "steelblue")
#' @param show_prior Mostrar distribucion prior (default = TRUE)
#' @param show_hdi Mostrar region HDI sombreada (default = TRUE)
#'
#' @return Invisible: lista con parametros de la posterior
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3)
#' plot_V_beta_binomial(ratings, l = 0, s = 3)
#'
#' @export
plot_V_beta_binomial <- function(ratings, l = 0, s = 3,
                                  prior_alpha = 1, prior_beta = 1,
                                  cred_level = 0.95,
                                  threshold = 0.70,
                                  main = NULL,
                                  col = "steelblue",
                                  show_prior = TRUE,
                                  show_hdi = TRUE) {

  # Calcular parametros
  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)
  k <- s - l

  # Exitos y fracasos
  sum_success <- sum(ratings - l)
  total_trials <- n * k

  # Posterior Beta
  post_alpha <- prior_alpha + sum_success
  post_beta <- prior_beta + (total_trials - sum_success)

  # Estadisticos
  V_mean <- post_alpha / (post_alpha + post_beta)
  V_classic <- (mean(ratings) - l) / k
  hdi <- compute_hdi(post_alpha, post_beta, cred_level)
  prob_threshold <- 1 - pbeta(threshold, post_alpha, post_beta)

  # Titulo

  if (is.null(main)) {
    main <- sprintf("Posterior Beta-Binomial (n = %d jueces)", n)
  }

  # Secuencia de valores
  x <- seq(0.001, 0.999, length.out = 500)
  y <- dbeta(x, post_alpha, post_beta)

  # Configurar grafico
  old_par <- par(mar = c(5, 5, 4, 2))
  on.exit(par(old_par))

  y_max <- max(y) * 1.15

  # Plot base
  plot(x, y, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "V de Aiken",
       ylab = "Densidad posterior",
       main = main,
       cex.lab = 1.1, cex.axis = 1, cex.main = 1.2,
       las = 1)

  # Grid
  grid(col = "gray90", lty = 1)

  # Prior (escalado)
  if (show_prior) {
    y_prior <- dbeta(x, prior_alpha, prior_beta)
    if (max(y_prior) > 0 && is.finite(max(y_prior))) {
      y_prior_scaled <- y_prior * (y_max * 0.25) / max(y_prior)
      lines(x, y_prior_scaled, col = "gray60", lwd = 2, lty = 2)
    }
  }

  # HDI shading
  if (show_hdi) {
    x_hdi <- x[x >= hdi[1] & x <= hdi[2]]
    y_hdi <- dbeta(x_hdi, post_alpha, post_beta)
    if (length(x_hdi) > 0) {
      polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
              c(0, y_hdi, 0),
              col = adjustcolor(col, alpha.f = 0.3),
              border = NA)
    }
  }

  # Densidad posterior
  lines(x, y, col = col, lwd = 3)

  # Lineas de referencia
  abline(v = threshold, col = "darkorange", lty = 2, lwd = 2)
  abline(v = V_mean, col = col, lty = 3, lwd = 2)
  abline(v = V_classic, col = "gray40", lty = 4, lwd = 1.5)

  # Texto con estadisticos
  text(0.02, y_max * 0.95,
       sprintf("V clasica = %.3f", V_classic),
       adj = c(0, 1), cex = 0.9, col = "gray40")
  text(0.02, y_max * 0.87,
       sprintf("V posterior = %.3f", V_mean),
       adj = c(0, 1), cex = 0.9, font = 2, col = col)
  text(0.02, y_max * 0.79,
       sprintf("%.0f%% HDI: [%.3f, %.3f]", cred_level * 100, hdi[1], hdi[2]),
       adj = c(0, 1), cex = 0.85, col = "gray30")
  text(0.02, y_max * 0.71,
       sprintf("P(V > %.2f) = %.3f", threshold, prob_threshold),
       adj = c(0, 1), cex = 0.9, font = 2, col = "darkorange")

  # Leyenda
  legend("topright",
         legend = c(
           sprintf("Posterior Beta(%.1f, %.1f)", post_alpha, post_beta),
           sprintf("%.0f%% HDI", cred_level * 100),
           sprintf("Umbral %.2f", threshold),
           if (show_prior) sprintf("Prior Beta(%.0f, %.0f)", prior_alpha, prior_beta) else NULL
         ),
         col = c(col, adjustcolor(col, 0.5), "darkorange", if (show_prior) "gray60" else NULL),
         lty = c(1, NA, 2, if (show_prior) 2 else NULL),
         lwd = c(3, NA, 2, if (show_prior) 2 else NULL),
         pch = c(NA, 15, NA, NA),
         pt.cex = c(NA, 2, NA, NA),
         bty = "n", cex = 0.8)

  # Retornar parametros
  invisible(list(
    modelo = "Beta-Binomial",
    n_jueces = n,
    V_clasica = V_classic,
    V_posterior = V_mean,
    post_alpha = post_alpha,
    post_beta = post_beta,
    HDI = hdi,
    prob_threshold = prob_threshold
  ))
}


#' Graficar densidad posterior del modelo Dirichlet-Multinomial
#'
#' Genera un grafico de la distribucion posterior de V usando el modelo
#' Dirichlet-Multinomial (mas conservador que Beta-Binomial).
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Parametro alpha del prior Dirichlet (default = 1, uniforme)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param n_samples Numero de muestras Monte Carlo (default = 10000)
#' @param threshold Umbral de referencia (default = 0.70)
#' @param main Titulo del grafico (opcional)
#' @param col Color de la densidad (default = "coral")
#' @param show_prior Mostrar distribucion prior (default = FALSE)
#' @param show_hdi Mostrar region HDI sombreada (default = TRUE)
#' @param seed Semilla para reproducibilidad (default = NULL)
#'
#' @return Invisible: lista con parametros y muestras de la posterior
#'
#' @details
#' El modelo Dirichlet-Multinomial modela directamente las frecuencias
#' de cada categoria, produciendo intervalos mas amplios y conservadores
#' que el modelo Beta-Binomial, especialmente con muestras pequenas.
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3)
#' plot_V_dirichlet(ratings, l = 0, s = 3)
#'
#' @export
plot_V_dirichlet <- function(ratings, l = 0, s = 3,
                              prior_alpha = 1,
                              cred_level = 0.95,
                              n_samples = 10000,
                              threshold = 0.70,
                              main = NULL,
                              col = "coral",
                              show_prior = FALSE,
                              show_hdi = TRUE,
                              seed = 42) {

  # Calcular parametros
  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)
  k <- s - l
  n_categories <- k + 1

  # V clasica
  V_classic <- (mean(ratings) - l) / k

  # Contar frecuencias
  ratings_idx <- ratings - l
  counts <- tabulate(ratings_idx + 1, nbins = n_categories)

  # Prior y posterior Dirichlet
  if (length(prior_alpha) == 1) {
    prior_alpha_vec <- rep(prior_alpha, n_categories)
  } else {
    prior_alpha_vec <- prior_alpha
  }
  post_alpha_vec <- prior_alpha_vec + counts

  # Simulacion Monte Carlo
  if (!is.null(seed)) set.seed(seed)

  V_samples <- numeric(n_samples)
  category_values <- 0:k

  for (i in 1:n_samples) {
    gamma_samples <- rgamma(n_categories, shape = post_alpha_vec, rate = 1)
    pi_samples <- gamma_samples / sum(gamma_samples)
    V_samples[i] <- sum(pi_samples * category_values) / k
  }

  # Estadisticos
  V_mean <- mean(V_samples)
  V_median <- median(V_samples)
  V_sd <- sd(V_samples)

  # HDI
  sorted_samples <- sort(V_samples)
  ci_mass <- floor(cred_level * n_samples)
  n_cis <- n_samples - ci_mass
  ci_widths <- sorted_samples[(ci_mass + 1):n_samples] - sorted_samples[1:n_cis]
  best_ci_idx <- which.min(ci_widths)
  hdi <- c(sorted_samples[best_ci_idx], sorted_samples[best_ci_idx + ci_mass])

  # P(V > threshold)
  prob_threshold <- mean(V_samples > threshold)

  # Titulo
  if (is.null(main)) {
    main <- sprintf("Posterior Dirichlet-Multinomial (n = %d jueces)", n)
  }

  # Densidad via kernel
  dens <- density(V_samples, from = 0, to = 1, n = 512)

  # Configurar grafico
  old_par <- par(mar = c(5, 5, 4, 2))
  on.exit(par(old_par))

  y_max <- max(dens$y) * 1.15

  # Plot base
  plot(dens$x, dens$y, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "V de Aiken",
       ylab = "Densidad posterior",
       main = main,
       cex.lab = 1.1, cex.axis = 1, cex.main = 1.2,
       las = 1)

  # Grid
  grid(col = "gray90", lty = 1)

  # HDI shading
  if (show_hdi) {
    x_hdi <- dens$x[dens$x >= hdi[1] & dens$x <= hdi[2]]
    y_hdi <- dens$y[dens$x >= hdi[1] & dens$x <= hdi[2]]
    if (length(x_hdi) > 0) {
      polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
              c(0, y_hdi, 0),
              col = adjustcolor(col, alpha.f = 0.3),
              border = NA)
    }
  }

  # Densidad posterior
  lines(dens$x, dens$y, col = col, lwd = 3)

  # Lineas de referencia
  abline(v = threshold, col = "darkorange", lty = 2, lwd = 2)
  abline(v = V_mean, col = col, lty = 3, lwd = 2)
  abline(v = V_classic, col = "gray40", lty = 4, lwd = 1.5)

  # Texto con estadisticos
  text(0.02, y_max * 0.95,
       sprintf("V clasica = %.3f", V_classic),
       adj = c(0, 1), cex = 0.9, col = "gray40")
  text(0.02, y_max * 0.87,
       sprintf("V posterior = %.3f (SD = %.3f)", V_mean, V_sd),
       adj = c(0, 1), cex = 0.9, font = 2, col = col)
  text(0.02, y_max * 0.79,
       sprintf("%.0f%% HDI: [%.3f, %.3f]", cred_level * 100, hdi[1], hdi[2]),
       adj = c(0, 1), cex = 0.85, col = "gray30")
  text(0.02, y_max * 0.71,
       sprintf("P(V > %.2f) = %.3f", threshold, prob_threshold),
       adj = c(0, 1), cex = 0.9, font = 2, col = "darkorange")

  # Frecuencias observadas
  freq_text <- paste0("n = (", paste(counts, collapse = ", "), ")")
  text(0.02, y_max * 0.60,
       sprintf("Frecuencias: %s", freq_text),
       adj = c(0, 1), cex = 0.8, col = "gray50")

  # Leyenda
  legend("topright",
         legend = c(
           "Posterior Dirichlet",
           sprintf("%.0f%% HDI", cred_level * 100),
           sprintf("Umbral %.2f", threshold)
         ),
         col = c(col, adjustcolor(col, 0.5), "darkorange"),
         lty = c(1, NA, 2),
         lwd = c(3, NA, 2),
         pch = c(NA, 15, NA),
         pt.cex = c(NA, 2, NA),
         bty = "n", cex = 0.8)

  # Retornar parametros
  invisible(list(
    modelo = "Dirichlet-Multinomial",
    n_jueces = n,
    n_categories = n_categories,
    counts = counts,
    V_clasica = V_classic,
    V_posterior = V_mean,
    V_sd = V_sd,
    post_alpha = post_alpha_vec,
    HDI = hdi,
    prob_threshold = prob_threshold,
    V_samples = V_samples
  ))
}
