#' =============================================================================
#' FUNCIONES AVANZADAS - BayesAiken
#' Analisis bayesianos avanzados para validez de contenido
#'
#' @author Jose Ventura-Leon
#' @version 1.0.0
#' =============================================================================

# =============================================================================
# 1. REGLAS DE DECISION PROBABILISTICAS
# =============================================================================

#' Regla de decision bayesiana para validez de contenido
#'
#' Evalua si un item cumple con criterios de validez usando probabilidades
#' posteriores en lugar de intervalos de confianza.
#'
#' @param result Objeto bayes_aiken o resultado de coef_V
#' @param threshold Umbral de validez (default = 0.70)
#' @param min_prob Probabilidad minima requerida (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con decision, probabilidad y justificacion
#'
#' @details
#' La regla de decision es:
#' - ACEPTAR: si P(V >= threshold) >= min_prob
#' - RECHAZAR: si P(V < threshold) >= min_prob
#' - INDECISO: en otro caso
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 3, 3, 2, 3, 3)
#' result <- coef_V(ratings, l = 0, s = 3)
#' decision <- decision_rule(result, threshold = 0.70, min_prob = 0.95)
#'
#' @export
decision_rule <- function(result,
                          threshold = 0.70,
                          min_prob = 0.95,
                          verbose = TRUE) {

  # Extraer parametros de la posterior
  if (inherits(result, "bayes_aiken")) {
    # Buscar el primer coeficiente disponible
    coef_name <- names(result$results)[1]
    post_alpha <- result$results[[coef_name]]$post_alpha
    post_beta <- result$results[[coef_name]]$post_beta
  } else if (is.list(result) && !is.null(result$post_alpha)) {
    post_alpha <- result$post_alpha
    post_beta <- result$post_beta
    coef_name <- "V"
  } else {
    stop("Formato de resultado no reconocido")
  }

  # Calcular probabilidades
  prob_above <- 1 - pbeta(threshold, post_alpha, post_beta)
  prob_below <- pbeta(threshold, post_alpha, post_beta)

  # Media posterior
  mean_post <- post_alpha / (post_alpha + post_beta)

  # Tomar decision
  if (prob_above >= min_prob) {
    decision <- "ACEPTAR"
    decision_color <- "green"
    justification <- sprintf(
      "P(%s >= %.2f) = %.3f >= %.2f",
      coef_name, threshold, prob_above, min_prob
    )
  } else if (prob_below >= min_prob) {
    decision <- "RECHAZAR"
    decision_color <- "red"
    justification <- sprintf(
      "P(%s < %.2f) = %.3f >= %.2f",
      coef_name, threshold, prob_below, min_prob
    )
  } else {
    decision <- "INDECISO"
    decision_color <- "orange"
    justification <- sprintf(
      "P(%s >= %.2f) = %.3f y P(%s < %.2f) = %.3f, ambas < %.2f",
      coef_name, threshold, prob_above,
      coef_name, threshold, prob_below, min_prob
    )
  }

  # Resultado
  output <- list(
    decision = decision,
    threshold = threshold,
    min_prob = min_prob,
    prob_above_threshold = prob_above,
    prob_below_threshold = prob_below,
    mean_posterior = mean_post,
    justification = justification,
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("        REGLA DE DECISION BAYESIANA\n")
    cat("================================================================\n\n")
    cat(sprintf("  Umbral de validez:     %.2f\n", threshold))
    cat(sprintf("  Probabilidad minima:   %.2f\n", min_prob))
    cat("\n")
    cat(sprintf("  Media posterior:       %.3f\n", mean_post))
    cat(sprintf("  P(%s >= %.2f):         %.3f\n", coef_name, threshold, prob_above))
    cat(sprintf("  P(%s <  %.2f):         %.3f\n", coef_name, threshold, prob_below))
    cat("\n")
    cat(sprintf("  DECISION: %s\n", decision))
    cat(sprintf("  Justificacion: %s\n", justification))
    cat("\n")
    cat("================================================================\n")
  }

  class(output) <- "bayes_decision"
  return(output)
}


#' Aplicar regla de decision a multiples items
#'
#' @param results Objeto bayes_aiken con multiples items o lista de resultados
#' @param threshold Umbral de validez (default = 0.70)
#' @param min_prob Probabilidad minima requerida (default = 0.95)
#'
#' @return Data.frame con decisiones para cada item
#'
#' @export
decision_rule_multi <- function(results,
                                threshold = 0.70,
                                min_prob = 0.95) {

  # Si es objeto bayes_aiken con data.frame
  if (inherits(results, "bayes_aiken") && results$input_type == "data.frame") {
    # Obtener tabla de resultados
    coef_name <- names(results$results)[1]
    df <- results$results[[coef_name]]

    decisions <- lapply(1:nrow(df), function(i) {
      # Reconstruir parametros
      mean_val <- df[i, paste0(coef_name, "_bayesiano")]
      if (is.null(mean_val) || is.na(mean_val)) {
        mean_val <- df[i, "V_bayesiano"]
      }
      n <- df[i, "n_jueces"]

      # Aproximar alpha y beta
      post_alpha <- mean_val * n * 3 + 1
      post_beta <- (1 - mean_val) * n * 3 + 1

      prob_above <- 1 - pbeta(threshold, post_alpha, post_beta)
      prob_below <- pbeta(threshold, post_alpha, post_beta)

      if (prob_above >= min_prob) {
        decision <- "ACEPTAR"
      } else if (prob_below >= min_prob) {
        decision <- "RECHAZAR"
      } else {
        decision <- "INDECISO"
      }

      data.frame(
        item = df[i, "item"],
        criterio = df[i, "criterio"],
        mean = mean_val,
        prob_above = prob_above,
        decision = decision,
        stringsAsFactors = FALSE
      )
    })

    return(do.call(rbind, decisions))
  }

  # Si es lista de resultados individuales
  if (is.list(results) && !inherits(results, "bayes_aiken")) {
    decisions <- lapply(seq_along(results), function(i) {
      res <- results[[i]]
      dec <- decision_rule(res, threshold, min_prob, verbose = FALSE)

      data.frame(
        item = names(results)[i],
        mean = dec$mean_posterior,
        prob_above = dec$prob_above_threshold,
        decision = dec$decision,
        stringsAsFactors = FALSE
      )
    })

    return(do.call(rbind, decisions))
  }

  stop("Formato de resultados no reconocido")
}


# =============================================================================
# 2. COMPARACION DE ITEMS/DIMENSIONES
# =============================================================================

#' Comparar dos items usando diferencias posteriores
#'
#' Calcula la distribucion posterior de la diferencia entre dos coeficientes
#' y reporta P(Delta > 0) o P(Delta > delta).
#'
#' @param result1 Resultado del primer item (objeto de coef_V o lista)
#' @param result2 Resultado del segundo item
#' @param delta Diferencia minima relevante (default = 0)
#' @param n_samples Numero de muestras Monte Carlo (default = 10000)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con distribucion de diferencias y probabilidades
#'
#' @examples
#' ratings1 <- c(3, 3, 3, 3, 2, 3, 3, 3, 3, 3)
#' ratings2 <- c(2, 3, 2, 3, 2, 2, 3, 2, 3, 2)
#' res1 <- coef_V(ratings1, l = 0, s = 3, verbose = FALSE)
#' res2 <- coef_V(ratings2, l = 0, s = 3, verbose = FALSE)
#' compare_items(res1, res2)
#'
#' @export
compare_items <- function(result1, result2,
                          delta = 0,
                          n_samples = 10000,
                          verbose = TRUE) {

  # Extraer parametros
  params1 <- .extract_posterior_params(result1)
  params2 <- .extract_posterior_params(result2)

  # Generar muestras de las posteriors
  set.seed(123)  # Para reproducibilidad
  samples1 <- rbeta(n_samples, params1$alpha, params1$beta)
  samples2 <- rbeta(n_samples, params2$alpha, params2$beta)

  # Calcular diferencias
  diff_samples <- samples1 - samples2

  # Estadisticos
  mean_diff <- mean(diff_samples)
  sd_diff <- sd(diff_samples)
  hdi_diff <- quantile(diff_samples, c(0.025, 0.975))

  # Probabilidades
  prob_positive <- mean(diff_samples > 0)
  prob_greater_delta <- mean(diff_samples > delta)
  prob_less_neg_delta <- mean(diff_samples < -delta)

  # Resultado
  output <- list(
    item1 = list(
      mean = params1$alpha / (params1$alpha + params1$beta),
      alpha = params1$alpha,
      beta = params1$beta
    ),
    item2 = list(
      mean = params2$alpha / (params2$alpha + params2$beta),
      alpha = params2$alpha,
      beta = params2$beta
    ),
    difference = list(
      mean = mean_diff,
      sd = sd_diff,
      hdi_95 = hdi_diff,
      samples = diff_samples
    ),
    probabilities = list(
      P_diff_positive = prob_positive,
      P_diff_greater_delta = prob_greater_delta,
      P_diff_less_neg_delta = prob_less_neg_delta,
      delta = delta
    ),
    n_samples = n_samples
  )

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("        COMPARACION DE ITEMS (Diferencias Posteriores)\n")
    cat("================================================================\n\n")
    cat(sprintf("  Item 1 - Media posterior:    %.3f\n", output$item1$mean))
    cat(sprintf("  Item 2 - Media posterior:    %.3f\n", output$item2$mean))
    cat("\n")
    cat("  Diferencia (Item1 - Item2):\n")
    cat(sprintf("    Media:                     %.3f\n", mean_diff))
    cat(sprintf("    Desv. Estandar:            %.3f\n", sd_diff))
    cat(sprintf("    95%% HDI:                   [%.3f, %.3f]\n", hdi_diff[1], hdi_diff[2]))
    cat("\n")
    cat("  Probabilidades:\n")
    cat(sprintf("    P(Item1 > Item2):          %.3f\n", prob_positive))
    if (delta > 0) {
      cat(sprintf("    P(Diff > %.2f):             %.3f\n", delta, prob_greater_delta))
      cat(sprintf("    P(Diff < -%.2f):            %.3f\n", delta, prob_less_neg_delta))
    }
    cat("\n")

    # Conclusion
    if (prob_positive >= 0.95) {
      cat("  CONCLUSION: Item 1 es superior a Item 2 (probabilidad >= 95%)\n")
    } else if (prob_positive <= 0.05) {
      cat("  CONCLUSION: Item 2 es superior a Item 1 (probabilidad >= 95%)\n")
    } else {
      cat("  CONCLUSION: No hay diferencia clara entre los items\n")
    }
    cat("\n")
    cat("================================================================\n")
  }

  class(output) <- "bayes_comparison"
  return(output)
}


#' Plot de comparacion de items
#'
#' @param comparison Objeto de compare_items
#' @param show_rope Mostrar ROPE (default = FALSE)
#' @param rope Limites del ROPE si show_rope = TRUE
#'
#' @export
plot.bayes_comparison <- function(x, show_rope = FALSE, rope = c(-0.05, 0.05), ...) {

  diff_samples <- x$difference$samples

  # Histograma de diferencias
  old_par <- par(mar = c(5, 4, 4, 2))
  on.exit(par(old_par))

  hist(diff_samples,
       breaks = 50,
       col = rgb(0.27, 0.51, 0.71, 0.5),
       border = "white",
       main = "Distribucion Posterior de la Diferencia",
       xlab = expression(Delta ~ "= Item1 - Item2"),
       ylab = "Frecuencia",
       las = 1)

  # Linea en cero
  abline(v = 0, col = "red", lwd = 2, lty = 2)

  # Linea en media
  abline(v = x$difference$mean, col = "steelblue", lwd = 2)

  # HDI
  abline(v = x$difference$hdi_95[1], col = "steelblue", lwd = 1.5, lty = 2)
  abline(v = x$difference$hdi_95[2], col = "steelblue", lwd = 1.5, lty = 2)

  # ROPE si se solicita
  if (show_rope) {
    rect(rope[1], 0, rope[2], par("usr")[4], col = rgb(1, 0.5, 0, 0.2), border = NA)
    abline(v = rope, col = "darkorange", lwd = 1.5, lty = 3)
  }

  # Leyenda
  legend("topright",
         legend = c(
           sprintf("Media = %.3f", x$difference$mean),
           "Cero",
           "95% HDI",
           if(show_rope) "ROPE" else NULL
         ),
         col = c("steelblue", "red", "steelblue", if(show_rope) "darkorange" else NULL),
         lty = c(1, 2, 2, if(show_rope) 3 else NULL),
         lwd = c(2, 2, 1.5, if(show_rope) 1.5 else NULL),
         bty = "n")

  # Texto con probabilidad
  text(x$difference$mean, par("usr")[4] * 0.9,
       sprintf("P(Item1 > Item2) = %.3f", x$probabilities$P_diff_positive),
       pos = 3, cex = 0.9)
}


# =============================================================================
# 3. COMPARACION DE GRUPOS DE JUECES
# =============================================================================

#' Comparar grupos de jueces
#'
#' Compara las evaluaciones de dos grupos independientes de jueces
#' (ej. clinicos vs investigadores, Lima vs provincias).
#'
#' @param ratings_group1 Vector o matriz de calificaciones del grupo 1
#' @param ratings_group2 Vector o matriz de calificaciones del grupo 2
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param delta Diferencia minima relevante para equivalencia (default = 0.05)
#' @param prior_alpha Alpha del prior (default = 1)
#' @param prior_beta Beta del prior (default = 1)
#' @param n_samples Numero de muestras Monte Carlo (default = 10000)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Lista con comparacion de grupos
#'
#' @examples
#' # Grupo 1: Clinicos
#' grupo1 <- c(3, 3, 2, 3, 3)
#' # Grupo 2: Investigadores
#' grupo2 <- c(3, 2, 3, 3, 2)
#' compare_judge_groups(grupo1, grupo2, l = 0, s = 3)
#'
#' @export
compare_judge_groups <- function(ratings_group1,
                                  ratings_group2,
                                  l = 0,
                                  s = 3,
                                  delta = 0.05,
                                  prior_alpha = 1,
                                  prior_beta = 1,
                                  n_samples = 10000,
                                  verbose = TRUE) {

  # Calcular V para cada grupo
  n1 <- length(ratings_group1)
  n2 <- length(ratings_group2)

  # Exitos = suma de (rating - l)
  successes1 <- sum(ratings_group1 - l)
  successes2 <- sum(ratings_group2 - l)

  # Total posible
  total1 <- n1 * (s - l)
  total2 <- n2 * (s - l)

  # Parametros posterior
  alpha1 <- prior_alpha + successes1
  beta1 <- prior_beta + (total1 - successes1)
  alpha2 <- prior_alpha + successes2
  beta2 <- prior_beta + (total2 - successes2)

  # Medias
  mean1 <- alpha1 / (alpha1 + beta1)
  mean2 <- alpha2 / (alpha2 + beta2)

  # Muestras Monte Carlo
  set.seed(123)
  samples1 <- rbeta(n_samples, alpha1, beta1)
  samples2 <- rbeta(n_samples, alpha2, beta2)
  diff_samples <- samples1 - samples2

  # Probabilidades
  prob_g1_greater <- mean(diff_samples > 0)
  prob_g2_greater <- mean(diff_samples < 0)
  prob_equivalent <- mean(abs(diff_samples) < delta)

  # HDI de la diferencia
  hdi_diff <- quantile(diff_samples, c(0.025, 0.975))

  output <- list(
    group1 = list(
      n = n1,
      ratings = ratings_group1,
      mean = mean1,
      alpha = alpha1,
      beta = beta1
    ),
    group2 = list(
      n = n2,
      ratings = ratings_group2,
      mean = mean2,
      alpha = alpha2,
      beta = beta2
    ),
    difference = list(
      mean = mean(diff_samples),
      sd = sd(diff_samples),
      hdi_95 = hdi_diff,
      samples = diff_samples
    ),
    probabilities = list(
      P_group1_greater = prob_g1_greater,
      P_group2_greater = prob_g2_greater,
      P_equivalent = prob_equivalent,
      delta = delta
    )
  )

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("        COMPARACION DE GRUPOS DE JUECES\n")
    cat("================================================================\n\n")
    cat(sprintf("  Grupo 1: n = %d jueces, V = %.3f\n", n1, mean1))
    cat(sprintf("  Grupo 2: n = %d jueces, V = %.3f\n", n2, mean2))
    cat("\n")
    cat("  Diferencia (Grupo1 - Grupo2):\n")
    cat(sprintf("    Media:                     %.3f\n", mean(diff_samples)))
    cat(sprintf("    95%% HDI:                   [%.3f, %.3f]\n", hdi_diff[1], hdi_diff[2]))
    cat("\n")
    cat("  Probabilidades:\n")
    cat(sprintf("    P(Grupo1 > Grupo2):        %.3f\n", prob_g1_greater))
    cat(sprintf("    P(Grupo2 > Grupo1):        %.3f\n", prob_g2_greater))
    cat(sprintf("    P(|Diff| < %.2f):          %.3f  (equivalencia)\n", delta, prob_equivalent))
    cat("\n")

    # Conclusion
    if (prob_equivalent >= 0.95) {
      cat("  CONCLUSION: Los grupos son practicamente equivalentes\n")
    } else if (prob_g1_greater >= 0.95) {
      cat("  CONCLUSION: Grupo 1 evalua mas alto que Grupo 2\n")
    } else if (prob_g2_greater >= 0.95) {
      cat("  CONCLUSION: Grupo 2 evalua mas alto que Grupo 1\n")
    } else {
      cat("  CONCLUSION: Evidencia inconclusa\n")
    }
    cat("\n")
    cat("================================================================\n")
  }

  class(output) <- c("bayes_group_comparison", "bayes_comparison")
  return(output)
}


# =============================================================================
# 4. ANALISIS ROPE (Region of Practical Equivalence)
# =============================================================================

#' Analisis ROPE para equivalencia practica
#'
#' Evalua si un coeficiente cae dentro de una region de equivalencia practica
#' usando el enfoque HDI + ROPE.
#'
#' @param result Objeto bayes_aiken o resultado de coef_V
#' @param rope_lower Limite inferior del ROPE (default = 0.70)
#' @param rope_upper Limite superior del ROPE (default = 1.00, es decir "al menos 0.70")
#' @param cred_level Nivel de credibilidad para HDI (default = 0.95)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Lista con decision ROPE
#'
#' @details
#' Decisiones posibles:
#' - ACEPTAR: Todo el HDI esta dentro del ROPE
#' - RECHAZAR: Todo el HDI esta fuera del ROPE
#' - INDECISO: El HDI se superpone parcialmente con el ROPE
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 3, 3, 2, 3, 3)
#' result <- coef_V(ratings, l = 0, s = 3)
#' rope_analysis(result, rope_lower = 0.70, rope_upper = 1.00)
#'
#' @export
rope_analysis <- function(result,
                          rope_lower = 0.70,
                          rope_upper = 1.00,
                          cred_level = 0.95,
                          verbose = TRUE) {

  # Extraer parametros
  params <- .extract_posterior_params(result)

  # Calcular HDI
  hdi <- compute_hdi(params$alpha, params$beta, cred_level)

  # Media posterior
  mean_post <- params$alpha / (params$alpha + params$beta)

  # Porcentaje de la posterior dentro del ROPE
  prob_in_rope <- pbeta(rope_upper, params$alpha, params$beta) -
                  pbeta(rope_lower, params$alpha, params$beta)

  # Decision basada en HDI + ROPE
  if (hdi[1] >= rope_lower && hdi[2] <= rope_upper) {
    decision <- "ACEPTAR"
    explanation <- "Todo el HDI esta dentro del ROPE"
  } else if (hdi[2] < rope_lower || hdi[1] > rope_upper) {
    decision <- "RECHAZAR"
    explanation <- "Todo el HDI esta fuera del ROPE"
  } else {
    decision <- "INDECISO"
    explanation <- "El HDI se superpone parcialmente con el ROPE"
  }

  output <- list(
    mean_posterior = mean_post,
    hdi = hdi,
    rope = c(rope_lower, rope_upper),
    prob_in_rope = prob_in_rope,
    decision = decision,
    explanation = explanation,
    cred_level = cred_level,
    post_alpha = params$alpha,
    post_beta = params$beta
  )

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("        ANALISIS ROPE (Region de Equivalencia Practica)\n")
    cat("================================================================\n\n")
    cat(sprintf("  ROPE:                        [%.2f, %.2f]\n", rope_lower, rope_upper))
    cat(sprintf("  Nivel de credibilidad:       %.0f%%\n", cred_level * 100))
    cat("\n")
    cat(sprintf("  Media posterior:             %.3f\n", mean_post))
    cat(sprintf("  %.0f%% HDI:                    [%.3f, %.3f]\n",
                cred_level * 100, hdi[1], hdi[2]))
    cat(sprintf("  P(V en ROPE):                %.3f\n", prob_in_rope))
    cat("\n")
    cat(sprintf("  DECISION: %s\n", decision))
    cat(sprintf("  Explicacion: %s\n", explanation))
    cat("\n")
    cat("================================================================\n")
  }

  class(output) <- "bayes_rope"
  return(output)
}


#' Plot del analisis ROPE
#'
#' @param x Objeto bayes_rope
#' @param ... Argumentos adicionales
#'
#' @export
plot.bayes_rope <- function(x, ...) {

  # Secuencia de valores
  vals <- seq(0.001, 0.999, length.out = 500)
  dens <- dbeta(vals, x$post_alpha, x$post_beta)

  # Configurar grafico
  old_par <- par(mar = c(5, 4, 4, 2))
  on.exit(par(old_par))

  y_max <- max(dens) * 1.1

  plot(vals, dens, type = "n",
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "Coeficiente V",
       ylab = "Densidad",
       main = "Analisis ROPE - HDI + Region de Equivalencia",
       las = 1)

  # ROPE region (sombreado amarillo)
  rect(x$rope[1], 0, x$rope[2], y_max, col = rgb(1, 0.8, 0, 0.3), border = NA)

  # HDI shading (azul)
  x_hdi <- vals[vals >= x$hdi[1] & vals <= x$hdi[2]]
  y_hdi <- dbeta(x_hdi, x$post_alpha, x$post_beta)
  if (length(x_hdi) > 0) {
    polygon(c(x_hdi[1], x_hdi, x_hdi[length(x_hdi)]),
            c(0, y_hdi, 0),
            col = rgb(0.27, 0.51, 0.71, 0.4),
            border = NA)
  }

  # Posterior
  lines(vals, dens, col = "steelblue", lwd = 3)

  # Lineas ROPE
  abline(v = x$rope[1], col = "darkorange", lwd = 2, lty = 2)
  abline(v = x$rope[2], col = "darkorange", lwd = 2, lty = 2)

  # Media
  abline(v = x$mean_posterior, col = "steelblue", lwd = 2, lty = 3)

  # Leyenda
  legend("topleft",
         legend = c(
           "Posterior",
           sprintf("%.0f%% HDI", x$cred_level * 100),
           "ROPE",
           sprintf("Decision: %s", x$decision)
         ),
         col = c("steelblue", rgb(0.27, 0.51, 0.71, 0.5), "darkorange", "black"),
         lty = c(1, NA, 2, NA),
         lwd = c(3, NA, 2, NA),
         pch = c(NA, 15, NA, NA),
         pt.cex = c(NA, 2, NA, NA),
         bty = "n")
}


# =============================================================================
# 5. PLANEAMIENTO DE MUESTRA (Posterior Predictivo)
# =============================================================================

#' Planear numero de jueces usando posterior predictivo
#'
#' Simula que pasaria si se agregan mas jueces y calcula la probabilidad
#' de alcanzar el umbral de validez con diferentes tamanos de panel.
#'
#' @param current_result Resultado actual (objeto bayes_aiken o coef_V)
#' @param target_prob Probabilidad objetivo de superar el umbral (default = 0.95)
#' @param threshold Umbral de validez (default = 0.70)
#' @param additional_judges Vector de jueces adicionales a evaluar (default = 1:10)
#' @param n_simulations Numero de simulaciones (default = 1000)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Data.frame con probabilidades por numero de jueces
#'
#' @examples
#' ratings <- c(3, 2, 3, 3, 2)  # 5 jueces actuales
#' result <- coef_V(ratings, l = 0, s = 3)
#' plan <- sample_size_planning(result, target_prob = 0.95, threshold = 0.70)
#'
#' @export
sample_size_planning <- function(current_result,
                                  target_prob = 0.95,
                                  threshold = 0.70,
                                  additional_judges = 1:10,
                                  n_simulations = 2000,
                                  verbose = TRUE) {

  # Extraer parametros actuales
  params <- .extract_posterior_params(current_result)
  current_mean <- params$alpha / (params$alpha + params$beta)
  current_n <- params$n

  # Escala (asumiendo 0-3 por defecto)
  scale_range <- 3

  # Probabilidad actual de superar el umbral
  current_prob <- 1 - pbeta(threshold, params$alpha, params$beta)

  # Maximo de jueces adicionales
  max_judges <- max(additional_judges)

  # Matrices para resultados
  # prob_matrix[sim, m] = P(V >= threshold | datos actuales + m jueces nuevos)
  prob_matrix <- matrix(NA, nrow = n_simulations, ncol = max_judges)
  true_V_vector <- numeric(n_simulations)

  set.seed(123)

  for (sim in 1:n_simulations) {
    # Muestrear V verdadero de la posterior actual
    true_V <- rbeta(1, params$alpha, params$beta)
    true_V_vector[sim] <- true_V

    # Simular calificaciones de max_judges jueces nuevos
    all_ratings <- rbinom(max_judges, scale_range, true_V)
    cumsum_ratings <- cumsum(all_ratings)

    for (m in 1:max_judges) {
      new_successes <- cumsum_ratings[m]
      new_trials <- m * scale_range

      # Posterior actualizada
      post_alpha_new <- params$alpha + new_successes
      post_beta_new <- params$beta + (new_trials - new_successes)

      # P(V >= threshold | datos actualizados)
      prob_matrix[sim, m] <- 1 - pbeta(threshold, post_alpha_new, post_beta_new)
    }
  }

  # Construir tabla de resultados
  results <- data.frame(
    n_additional = additional_judges,
    n_total = current_n + additional_judges,
    E_prob = numeric(length(additional_judges)),
    P_accept = numeric(length(additional_judges)),
    P_reject = numeric(length(additional_judges)),
    P_unclear = numeric(length(additional_judges)),
    meets_target = logical(length(additional_judges))
  )

  for (i in seq_along(additional_judges)) {
    m <- additional_judges[i]
    probs <- prob_matrix[, m]

    # Esperanza de P(V >= threshold)
    results$E_prob[i] <- mean(probs)
    # Proporcion que lleva a ACEPTAR (P >= target_prob)
    results$P_accept[i] <- mean(probs >= target_prob)
    # Proporcion que lleva a RECHAZAR (P <= 1 - target_prob)
    results$P_reject[i] <- mean(probs <= (1 - target_prob))
    # Proporcion indecisa
    results$P_unclear[i] <- 1 - results$P_accept[i] - results$P_reject[i]
    # Cumple objetivo si P(aceptar) >= 0.80
    results$meets_target[i] <- results$P_accept[i] >= 0.80
  }

  output <- list(
    current = list(
      n = current_n,
      mean = current_mean,
      prob_above_threshold = current_prob
    ),
    target_prob = target_prob,
    threshold = threshold,
    planning_table = results,
    recommended_n = NULL,
    simulation_details = list(
      true_V = true_V_vector,
      prob_matrix = prob_matrix
    )
  )

  # Encontrar n recomendado
  # 1) Primer n donde P_accept >= 80% (ideal)
  # 2) Si no, primer n donde P_accept >= 70% (aceptable)
  # 3) Si no, el n con mayor P_accept
  meets_rows <- results[results$meets_target, ]
  if (nrow(meets_rows) > 0) {
    output$recommended_n <- min(meets_rows$n_additional)
    output$recommendation_type <- "optimal"
  } else {
    # Buscar primer n donde P_accept >= 70%
    acceptable_rows <- results[results$P_accept >= 0.70, ]
    if (nrow(acceptable_rows) > 0) {
      output$recommended_n <- min(acceptable_rows$n_additional)
      output$recommendation_type <- "acceptable"
    } else {
      # Si no, el n con mayor P_accept
      best_idx <- which.max(results$P_accept)
      output$recommended_n <- results$n_additional[best_idx]
      output$recommendation_type <- "best_available"
    }
  }

  if (verbose) {
    cat("\n")
    cat("================================================================\n")
    cat("        PLANIFICACION DE MUESTRA (Posterior Predictivo)\n")
    cat("================================================================\n\n")
    cat(sprintf("  Jueces actuales:             %d\n", current_n))
    cat(sprintf("  V actual (media):            %.3f\n", current_mean))
    cat(sprintf("  P(V >= %.2f) actual:         %.3f\n", threshold, current_prob))
    cat("\n")
    cat(sprintf("  Criterio de decision: P(V >= %.2f) >= %.2f\n", threshold, target_prob))
    cat("\n")
    cat("  Proyeccion:\n")
    cat("  ----------------------------------------------------------------\n")
    cat("  n_add  n_total  P(Aceptar)  P(Rechazar)  P(Indeciso)\n")
    cat("  ----------------------------------------------------------------\n")
    for (i in 1:nrow(results)) {
      cat(sprintf("  %5d  %7d  %10.1f%%  %11.1f%%  %11.1f%%\n",
                  results$n_additional[i],
                  results$n_total[i],
                  results$P_accept[i] * 100,
                  results$P_reject[i] * 100,
                  results$P_unclear[i] * 100))
    }
    cat("  ----------------------------------------------------------------\n")
    cat("\n")

    rec_row <- results[results$n_additional == output$recommended_n, ]

    cat(sprintf("  RECOMENDACION: Agregar %d juez(ces) mas (total = %d)\n",
                output$recommended_n, current_n + output$recommended_n))
    cat(sprintf("                 Probabilidad de aceptacion: %.1f%%\n",
                rec_row$P_accept * 100))

    if (output$recommendation_type == "acceptable") {
      cat("\n")
      cat("  NOTA: Recomendacion basada en P(Aceptar) >= 70%%.\n")
      cat("        El criterio ideal (>= 80%%) no se alcanza en el rango.\n")
    } else if (output$recommendation_type == "best_available") {
      cat("\n")
      cat(sprintf("  NOTA: Con V = %.3f (cercano al umbral %.2f),\n",
                  current_mean, threshold))
      cat("        P(Aceptar) no alcanza 70%% en el rango evaluado.\n")
      cat("        Se recomienda el n con mayor probabilidad disponible.\n")
    }
    cat("\n")
    cat("================================================================\n")
  }

  class(output) <- "bayes_sample_plan"
  return(output)
}


#' Plot del planeamiento de muestra
#'
#' @param x Objeto bayes_sample_plan
#' @param ... Argumentos adicionales
#'
#' @export
plot.bayes_sample_plan <- function(x, ...) {

  df <- x$planning_table

  old_par <- par(mar = c(5, 4, 4, 2))
  on.exit(par(old_par))

  # Plot de probabilidad vs numero de jueces
  plot(df$n_total, df$prob_above_threshold,
       type = "b",
       pch = 19,
       col = "steelblue",
       lwd = 2,
       xlab = "Numero total de jueces",
       ylab = sprintf("P(V >= %.2f)", x$threshold),
       main = "Planificacion de Tamano de Muestra",
       ylim = c(0, 1),
       las = 1)

  # Linea de objetivo
  abline(h = x$target_prob, col = "darkorange", lwd = 2, lty = 2)

  # Punto actual
  points(x$current$n, x$current$prob_above_threshold,
         pch = 17, col = "red", cex = 1.5)

  # Punto recomendado
  if (!is.null(x$recommended_n)) {
    rec_row <- df[df$n_additional == x$recommended_n, ]
    points(rec_row$n_total, rec_row$prob_above_threshold,
           pch = 19, col = "forestgreen", cex = 1.5)
  }

  # Leyenda
  legend("bottomright",
         legend = c(
           "Proyeccion",
           sprintf("Objetivo (%.2f)", x$target_prob),
           "Actual",
           if (!is.null(x$recommended_n)) "Recomendado" else NULL
         ),
         col = c("steelblue", "darkorange", "red",
                 if (!is.null(x$recommended_n)) "forestgreen" else NULL),
         pch = c(19, NA, 17, if (!is.null(x$recommended_n)) 19 else NULL),
         lty = c(1, 2, NA, NA),
         lwd = c(2, 2, NA, NA),
         bty = "n")

  # Grid
  grid(col = "gray90")
}


# =============================================================================
# 6. SIMULACION PARA PLANIFICACION DEL NUMERO DE JUECES (Tabla 3 del articulo)
# =============================================================================

#' Simulacion para planificacion del numero de jueces
#'
#' Realiza una simulacion Monte Carlo para evaluar como cambia la precision
#' de los intervalos y la probabilidad de cumplir una regla decisional
#' al variar el numero de jueces. Replica la metodologia de la Tabla 3
#' del articulo "V de Aiken Bayesiana".
#'
#' @param true_V Vector de valores verdaderos de V a simular (default = c(0.80, 0.90))
#' @param n_judges Vector de numeros de jueces a evaluar (default = c(5, 10, 15))
#' @param k Amplitud de la escala (s - l). Para escala 0-3, k = 3 (default)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad para intervalos (default = 0.95)
#' @param threshold Umbral de validez para la regla decisional (default = 0.70)
#' @param min_prob Probabilidad minima para la regla (default = 0.95)
#' @param n_sim Numero de simulaciones Monte Carlo (default = 1000)
#' @param seed Semilla para reproducibilidad (default = 123)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Objeto de clase 'bayes_simulation' con:
#' \itemize{
#'   \item results: Data.frame con resultados por V verdadero y n de jueces
#'   \item parameters: Lista con parametros de la simulacion
#'   \item raw_data: Lista con datos crudos de cada simulacion
#' }
#'
#' @details
#' Para cada combinacion de V verdadero y numero de jueces, la funcion:
#' \enumerate{
#'   \item Genera n_sim conjuntos de calificaciones simuladas
#'   \item Calcula la posterior Beta para cada simulacion
#'   \item Computa ETI y HDI al nivel cred_level
#'   \item Evalua cobertura (si el intervalo contiene V verdadero)
#'   \item Evalua si se cumple la regla P(V > threshold) >= min_prob
#' }
#'
#' Las metricas reportadas son:
#' \itemize{
#'   \item Cobertura ETI/HDI: Proporcion de simulaciones donde el intervalo contiene V verdadero
#'   \item Ancho ETI/HDI: Ancho promedio del intervalo de credibilidad
#'   \item Pr(regla): Proporcion de simulaciones que cumplen la regla decisional
#' }
#'
#' @examples
#' # Replicar Tabla 3 del articulo
#' sim <- simulate_judges_planning(
#'   true_V = c(0.80, 0.90),
#'   n_judges = c(5, 10, 15),
#'   k = 3,
#'   n_sim = 1000
#' )
#' print(sim)
#' plot(sim)
#'
#' # Simulacion mas extensa
#' sim2 <- simulate_judges_planning(
#'   true_V = seq(0.70, 0.95, by = 0.05),
#'   n_judges = c(5, 7, 10, 12, 15, 20),
#'   n_sim = 2000
#' )
#'
#' @export
simulate_judges_planning <- function(true_V = c(0.80, 0.90),
                                      n_judges = c(5, 10, 15),
                                      k = 3,
                                      prior_alpha = 1,
                                      prior_beta = 1,
                                      cred_level = 0.95,
                                      threshold = 0.70,
                                      min_prob = 0.95,
                                      n_sim = 1000,
                                      seed = 123,
                                      verbose = TRUE) {

  # Validaciones
  if (any(true_V <= 0) || any(true_V >= 1)) {
    stop("true_V debe estar en el intervalo (0, 1)")
  }
  if (any(n_judges < 2)) {
    stop("n_judges debe ser al menos 2")
  }
  if (k < 1) {
    stop("k (amplitud de escala) debe ser al menos 1")
  }

  set.seed(seed)

  # Preparar matriz de resultados
  n_combinations <- length(true_V) * length(n_judges)
  results <- data.frame(
    V_verdadero = numeric(n_combinations),
    n_jueces = integer(n_combinations),
    cobertura_ETI = numeric(n_combinations),
    cobertura_HDI = numeric(n_combinations),
    ancho_ETI = numeric(n_combinations),
    ancho_HDI = numeric(n_combinations),
    prob_regla = numeric(n_combinations),
    V_media_estimada = numeric(n_combinations),
    sesgo_medio = numeric(n_combinations)
  )

  # Lista para guardar datos crudos
  raw_data <- list()

  row_idx <- 1

  for (V_true in true_V) {
    for (n in n_judges) {

      # Total de ensayos binomiales por simulacion
      N <- n * k

      # Vectores para almacenar resultados de cada simulacion
      V_estimates <- numeric(n_sim)
      ETI_lower <- numeric(n_sim)
      ETI_upper <- numeric(n_sim)
      HDI_lower <- numeric(n_sim)
      HDI_upper <- numeric(n_sim)
      prob_above_threshold <- numeric(n_sim)
      covers_ETI <- logical(n_sim)
      covers_HDI <- logical(n_sim)
      meets_rule <- logical(n_sim)

      for (sim in 1:n_sim) {
        # Simular S ~ Binomial(N, V_true)
        S <- rbinom(1, N, V_true)

        # Parametros posterior
        post_alpha <- prior_alpha + S
        post_beta <- prior_beta + (N - S)

        # Estimador puntual (media posterior)
        V_estimates[sim] <- post_alpha / (post_alpha + post_beta)

        # ETI (Equal-Tailed Interval)
        alpha_level <- 1 - cred_level
        ETI_lower[sim] <- qbeta(alpha_level / 2, post_alpha, post_beta)
        ETI_upper[sim] <- qbeta(1 - alpha_level / 2, post_alpha, post_beta)

        # HDI (Highest Density Interval)
        hdi <- compute_hdi(post_alpha, post_beta, cred_level)
        HDI_lower[sim] <- hdi[1]
        HDI_upper[sim] <- hdi[2]

        # Cobertura
        covers_ETI[sim] <- (ETI_lower[sim] <= V_true) && (V_true <= ETI_upper[sim])
        covers_HDI[sim] <- (HDI_lower[sim] <= V_true) && (V_true <= HDI_upper[sim])

        # Probabilidad posterior de superar umbral
        prob_above_threshold[sim] <- 1 - pbeta(threshold, post_alpha, post_beta)

        # Cumple regla decisional
        meets_rule[sim] <- prob_above_threshold[sim] >= min_prob
      }

      # Calcular metricas
      results$V_verdadero[row_idx] <- V_true
      results$n_jueces[row_idx] <- n
      results$cobertura_ETI[row_idx] <- mean(covers_ETI)
      results$cobertura_HDI[row_idx] <- mean(covers_HDI)
      results$ancho_ETI[row_idx] <- mean(ETI_upper - ETI_lower)
      results$ancho_HDI[row_idx] <- mean(HDI_upper - HDI_lower)
      results$prob_regla[row_idx] <- mean(meets_rule)
      results$V_media_estimada[row_idx] <- mean(V_estimates)
      results$sesgo_medio[row_idx] <- mean(V_estimates) - V_true

      # Guardar datos crudos
      raw_data[[paste0("V", V_true, "_n", n)]] <- list(
        V_true = V_true,
        n = n,
        V_estimates = V_estimates,
        ETI = cbind(ETI_lower, ETI_upper),
        HDI = cbind(HDI_lower, HDI_upper),
        prob_above = prob_above_threshold,
        covers_ETI = covers_ETI,
        covers_HDI = covers_HDI,
        meets_rule = meets_rule
      )

      row_idx <- row_idx + 1
    }
  }

  # Crear objeto de salida
  output <- list(
    results = results,
    parameters = list(
      true_V = true_V,
      n_judges = n_judges,
      k = k,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      cred_level = cred_level,
      threshold = threshold,
      min_prob = min_prob,
      n_sim = n_sim,
      seed = seed
    ),
    raw_data = raw_data
  )

  class(output) <- "bayes_simulation"

  if (verbose) {
    print(output)
  }

  return(output)
}


#' Imprimir resultados de simulacion
#'
#' @param x Objeto bayes_simulation
#' @param ... Argumentos adicionales
#'
#' @export
print.bayes_simulation <- function(x, ...) {
  cat("\n")
  cat("================================================================\n")
  cat("   SIMULACION PARA PLANIFICACION DEL NUMERO DE JUECES\n")
  cat("================================================================\n\n")

  cat(sprintf("  Prior:                 Beta(%.1f, %.1f)\n",
              x$parameters$prior_alpha, x$parameters$prior_beta))
  cat(sprintf("  Escala (k):            %d (ej. 0-%d)\n",
              x$parameters$k, x$parameters$k))
  cat(sprintf("  Nivel credibilidad:    %.0f%%\n", x$parameters$cred_level * 100))
  cat(sprintf("  Umbral decision:       V > %.2f\n", x$parameters$threshold))
  cat(sprintf("  Probabilidad minima:   %.2f\n", x$parameters$min_prob))
  cat(sprintf("  N simulaciones:        %d\n", x$parameters$n_sim))
  cat("\n")

  cat("  Resultados:\n")
  cat("  ----------------------------------------------------------------\n")
  cat("  V verdadero   J   Cobertura ETI   Ancho ETI   Pr(regla)\n")
  cat("  ----------------------------------------------------------------\n")

  for (i in 1:nrow(x$results)) {
    cat(sprintf("      %.2f      %2d      %.3f         %.3f       %.3f\n",
                x$results$V_verdadero[i],
                x$results$n_jueces[i],
                x$results$cobertura_ETI[i],
                x$results$ancho_ETI[i],
                x$results$prob_regla[i]))
  }

  cat("  ----------------------------------------------------------------\n")
  cat("\n")

  # Interpretacion
  cat("  Interpretacion:\n")
  cat("  - Cobertura ETI: Proporcion de simulaciones donde el intervalo\n")
  cat("                   contiene el valor verdadero (deberia ser ~0.95)\n")
  cat("  - Ancho ETI: Precision del intervalo (menor = mas preciso)\n")
  cat("  - Pr(regla): Probabilidad de cumplir P(V > 0.70) >= 0.95\n")
  cat("\n")

  # Recomendaciones
  cat("  Observaciones:\n")

  # Para cada V verdadero
  for (V in unique(x$results$V_verdadero)) {
    subset_df <- x$results[x$results$V_verdadero == V, ]

    # Primer n donde Pr(regla) >= 0.80
    good_rows <- subset_df[subset_df$prob_regla >= 0.80, ]
    if (nrow(good_rows) > 0) {
      min_n <- min(good_rows$n_jueces)
      cat(sprintf("  - V = %.2f: Con %d jueces, Pr(regla) = %.1f%% (>= 80%%)\n",
                  V, min_n, good_rows$prob_regla[good_rows$n_jueces == min_n] * 100))
    } else {
      max_prob <- max(subset_df$prob_regla)
      best_n <- subset_df$n_jueces[which.max(subset_df$prob_regla)]
      cat(sprintf("  - V = %.2f: Maximo Pr(regla) = %.1f%% con %d jueces\n",
                  V, max_prob * 100, best_n))
    }
  }

  cat("\n")
  cat("================================================================\n")

  invisible(x)
}


#' Grafico de resultados de simulacion
#'
#' @param x Objeto bayes_simulation
#' @param type Tipo de grafico: "coverage" (cobertura), "width" (ancho),
#'        "power" (Pr regla), o "all" (todos, default)
#' @param ... Argumentos adicionales
#'
#' @export
plot.bayes_simulation <- function(x, type = "all", ...) {

  results <- x$results
  true_V_values <- unique(results$V_verdadero)
  n_V <- length(true_V_values)

  # Colores para diferentes V
  colors <- c("steelblue", "darkorange", "forestgreen", "purple", "red",
              "brown", "pink", "cyan")[1:n_V]

  if (type == "all") {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old_par))
  }

  # 1. Cobertura ETI
  if (type %in% c("all", "coverage")) {
    plot(NULL, xlim = range(results$n_jueces),
         ylim = c(0.90, 1),
         xlab = "Numero de jueces (J)",
         ylab = "Cobertura ETI 95%",
         main = "Cobertura del Intervalo de Credibilidad",
         las = 1)

    abline(h = x$parameters$cred_level, col = "gray50", lty = 2, lwd = 1.5)

    for (i in seq_along(true_V_values)) {
      V <- true_V_values[i]
      subset_df <- results[results$V_verdadero == V, ]
      lines(subset_df$n_jueces, subset_df$cobertura_ETI,
            type = "b", pch = 19, col = colors[i], lwd = 2)
    }

    legend("bottomright",
           legend = paste("V =", true_V_values),
           col = colors,
           pch = 19, lty = 1, lwd = 2,
           bty = "n", cex = 0.8)
  }

  # 2. Ancho ETI
  if (type %in% c("all", "width")) {
    plot(NULL, xlim = range(results$n_jueces),
         ylim = c(0, max(results$ancho_ETI) * 1.1),
         xlab = "Numero de jueces (J)",
         ylab = "Ancho ETI 95%",
         main = "Precision del Intervalo",
         las = 1)

    for (i in seq_along(true_V_values)) {
      V <- true_V_values[i]
      subset_df <- results[results$V_verdadero == V, ]
      lines(subset_df$n_jueces, subset_df$ancho_ETI,
            type = "b", pch = 19, col = colors[i], lwd = 2)
    }

    legend("topright",
           legend = paste("V =", true_V_values),
           col = colors,
           pch = 19, lty = 1, lwd = 2,
           bty = "n", cex = 0.8)
  }

  # 3. Pr(regla)
  if (type %in% c("all", "power")) {
    plot(NULL, xlim = range(results$n_jueces),
         ylim = c(0, 1),
         xlab = "Numero de jueces (J)",
         ylab = sprintf("Pr(V > %.2f) >= %.2f",
                        x$parameters$threshold, x$parameters$min_prob),
         main = "Probabilidad de Cumplir Regla Decisional",
         las = 1)

    abline(h = 0.80, col = "gray50", lty = 2, lwd = 1.5)
    abline(h = 0.95, col = "gray70", lty = 3, lwd = 1)

    for (i in seq_along(true_V_values)) {
      V <- true_V_values[i]
      subset_df <- results[results$V_verdadero == V, ]
      lines(subset_df$n_jueces, subset_df$prob_regla,
            type = "b", pch = 19, col = colors[i], lwd = 2)
    }

    legend("bottomright",
           legend = paste("V =", true_V_values),
           col = colors,
           pch = 19, lty = 1, lwd = 2,
           bty = "n", cex = 0.8)
  }

  # 4. Sesgo
  if (type == "all") {
    plot(NULL, xlim = range(results$n_jueces),
         ylim = c(min(results$sesgo_medio) - 0.01,
                  max(results$sesgo_medio) + 0.01),
         xlab = "Numero de jueces (J)",
         ylab = "Sesgo medio (E[V] - V verdadero)",
         main = "Sesgo del Estimador",
         las = 1)

    abline(h = 0, col = "gray50", lty = 2, lwd = 1.5)

    for (i in seq_along(true_V_values)) {
      V <- true_V_values[i]
      subset_df <- results[results$V_verdadero == V, ]
      lines(subset_df$n_jueces, subset_df$sesgo_medio,
            type = "b", pch = 19, col = colors[i], lwd = 2)
    }

    legend("topright",
           legend = paste("V =", true_V_values),
           col = colors,
           pch = 19, lty = 1, lwd = 2,
           bty = "n", cex = 0.8)
  }
}


#' Convertir resultados de simulacion a data.frame (formato Tabla 3)
#'
#' @param x Objeto bayes_simulation
#' @param format Formato de salida: "wide" o "long" (default = "wide")
#'
#' @return Data.frame con resultados formateados
#'
#' @examples
#' sim <- simulate_judges_planning()
#' tabla <- as.data.frame(sim)
#' print(tabla)
#'
#' @export
as.data.frame.bayes_simulation <- function(x, row.names = NULL,
                                            optional = FALSE,
                                            format = "wide", ...) {
  df <- x$results

  if (format == "wide") {
    # Formato similar a Tabla 3 del articulo
    df_out <- data.frame(
      `V verdadero` = df$V_verdadero,
      J = df$n_jueces,
      `Cobertura ETI 95%` = round(df$cobertura_ETI, 3),
      `Ancho ETI 95%` = round(df$ancho_ETI, 3),
      `Pr(regla)` = round(df$prob_regla, 3),
      check.names = FALSE
    )
  } else {
    df_out <- df
  }

  return(df_out)
}


# =============================================================================
# FUNCIONES AUXILIARES INTERNAS
# =============================================================================

#' Extraer parametros de posterior de un resultado
#' @noRd
.extract_posterior_params <- function(result) {
  if (inherits(result, "bayes_aiken")) {
    coef_name <- names(result$results)[1]
    return(list(
      alpha = result$results[[coef_name]]$post_alpha,
      beta = result$results[[coef_name]]$post_beta,
      n = result$results[[coef_name]]$n_jueces
    ))
  } else if (is.list(result) && !is.null(result$post_alpha)) {
    return(list(
      alpha = result$post_alpha,
      beta = result$post_beta,
      n = if (!is.null(result$n_jueces)) result$n_jueces else 10
    ))
  } else {
    stop("Formato de resultado no reconocido")
  }
}
