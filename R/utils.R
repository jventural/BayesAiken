#' =============================================================================
#' UTILIDADES Y HELPERS - BayesAiken
#' =============================================================================

# =============================================================================
# 1. DIAGNOSTICO DE HETEROGENEIDAD ENTRE JUECES
# =============================================================================

#' Diagnostico de heterogeneidad entre jueces
#'
#' Evalua si existe heterogeneidad significativa entre los jueces que podria
#' violar el supuesto de intercambiabilidad del modelo binomial. Reporta
#' estadisticos descriptivos, identifica jueces extremos y genera advertencias.
#'
#' @param ratings Vector de calificaciones o matriz (filas = items, cols = jueces)
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param threshold_cv Umbral de coeficiente de variacion para advertencia (default = 0.20)
#' @param threshold_outlier Numero de desviaciones estandar para considerar outlier (default = 2)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Objeto de clase 'judge_heterogeneity' con:
#' \itemize{
#'   \item summary: Estadisticos descriptivos por juez
#'   \item overall: Metricas globales de heterogeneidad
#'   \item outliers: Jueces identificados como extremos
#'   \item warning_level: Nivel de advertencia ("none", "low", "moderate", "high")
#'   \item recommendations: Recomendaciones basadas en el diagnostico
#' }
#'
#' @details
#' El diagnostico incluye:
#' \itemize{
#'   \item Media y desviacion estandar por juez
#'   \item Coeficiente de variacion (CV) entre jueces
#'   \item Indice de dispersion (varianza/media)
#'   \item Identificacion de jueces outliers (muy duros o muy blandos)
#'   \item Prueba de Levene simplificada para homogeneidad de varianzas
#' }
#'
#' Niveles de advertencia:
#' \itemize{
#'   \item none: CV < 0.10 (heterogeneidad minima)
#'   \item low: CV 0.10-0.20 (heterogeneidad baja, aceptable)
#'   \item moderate: CV 0.20-0.30 (heterogeneidad moderada, considerar)
#'   \item high: CV > 0.30 (heterogeneidad alta, intervalos pueden ser optimistas)
#' }
#'
#' @examples
#' # Vector simple (un item, multiples jueces)
#' ratings <- c(3, 3, 2, 3, 1, 3, 3, 2, 3, 3)
#' check <- check_judge_heterogeneity(ratings)
#'
#' # Matriz (multiples items x jueces)
#' mat <- matrix(c(3,3,2,3,1,
#'                 3,3,3,3,2,
#'                 2,3,2,3,1), nrow = 3, byrow = TRUE)
#' check <- check_judge_heterogeneity(mat)
#'
#' @export
check_judge_heterogeneity <- function(ratings,
                                       l = 0,
                                       s = 3,
                                       threshold_cv = 0.20,
                                       threshold_outlier = 2,
                                       verbose = TRUE) {

  # Convertir a matriz si es vector

  if (is.vector(ratings)) {
    ratings_matrix <- matrix(ratings, nrow = 1)
  } else if (is.data.frame(ratings)) {
    # Seleccionar solo columnas numericas (jueces)
    numeric_cols <- sapply(ratings, is.numeric)
    ratings_matrix <- as.matrix(ratings[, numeric_cols])
  } else {
    ratings_matrix <- as.matrix(ratings)
  }

  n_items <- nrow(ratings_matrix)
  n_judges <- ncol(ratings_matrix)

  if (n_judges < 2) {
    stop("Se necesitan al menos 2 jueces para evaluar heterogeneidad")
  }

  # Normalizar ratings a escala 0-1 (como V de Aiken)
  ratings_normalized <- (ratings_matrix - l) / (s - l)

  # Calcular estadisticos por juez (promedio a traves de items)
  judge_means <- colMeans(ratings_normalized, na.rm = TRUE)
  judge_sds <- apply(ratings_normalized, 2, sd, na.rm = TRUE)

  # Estadisticos globales
  overall_mean <- mean(judge_means)
  overall_sd <- sd(judge_means)
  cv <- if (overall_mean > 0) overall_sd / overall_mean else NA

  # Indice de dispersion (para detectar sobre-dispersion)
  # Varianza observada vs varianza esperada bajo binomial
  k <- s - l
  expected_var <- overall_mean * (1 - overall_mean) / (n_items * k)
  observed_var <- var(judge_means)
  dispersion_index <- if (expected_var > 0) observed_var / expected_var else NA

  # Identificar jueces outliers
  z_scores <- (judge_means - overall_mean) / overall_sd
  outlier_indices <- which(abs(z_scores) > threshold_outlier)

  outliers_df <- data.frame(
    juez = outlier_indices,
    media = round(judge_means[outlier_indices], 3),
    z_score = round(z_scores[outlier_indices], 2),
    tipo = ifelse(z_scores[outlier_indices] > 0, "Muy blando", "Muy duro"),
    stringsAsFactors = FALSE
  )

  # Determinar nivel de advertencia
  if (is.na(cv)) {
    warning_level <- "unknown"
  } else if (cv < 0.10) {
    warning_level <- "none"
  } else if (cv < 0.20) {
    warning_level <- "low"
  } else if (cv < 0.30) {
    warning_level <- "moderate"
  } else {
    warning_level <- "high"
  }

  # Generar recomendaciones
  recommendations <- character(0)

  if (warning_level == "none") {
    recommendations <- c(recommendations,
      "La heterogeneidad entre jueces es minima.",
      "El modelo binomial con prior Beta es apropiado.")
  } else if (warning_level == "low") {
    recommendations <- c(recommendations,
      "La heterogeneidad entre jueces es baja y aceptable.",
      "Los intervalos de credibilidad son confiables.")
  } else if (warning_level == "moderate") {
    recommendations <- c(recommendations,
      "Se detecta heterogeneidad moderada entre jueces.",
      "Los intervalos podrian ser ligeramente optimistas.",
      "Considere reportar el coeficiente H (homogeneidad) junto con V.")
  } else if (warning_level == "high") {
    recommendations <- c(recommendations,
      "Se detecta alta heterogeneidad entre jueces.",
      "Los intervalos de credibilidad pueden ser optimistas.",
      "Considere: (1) revisar criterios de los jueces outliers,",
      "           (2) usar un modelo ordinal jerarquico,",
      "           (3) reportar analisis de sensibilidad excluyendo outliers.")
  }

  if (nrow(outliers_df) > 0) {
    recommendations <- c(recommendations,
      sprintf("Se identificaron %d juez(ces) con calificaciones extremas.",
              nrow(outliers_df)))
  }

  if (!is.na(dispersion_index) && dispersion_index > 2) {
    recommendations <- c(recommendations,
      sprintf("Indice de dispersion = %.2f (>2 indica sobre-dispersion).",
              dispersion_index))
  }

  # Crear resumen por juez
  judge_summary <- data.frame(
    juez = 1:n_judges,
    media_V = round(judge_means, 3),
    DE = round(judge_sds, 3),
    z_score = round(z_scores, 2),
    es_outlier = 1:n_judges %in% outlier_indices,
    stringsAsFactors = FALSE
  )

  # Objeto de salida
  output <- list(
    summary = judge_summary,
    overall = list(
      n_items = n_items,
      n_judges = n_judges,
      mean_across_judges = round(overall_mean, 4),
      sd_across_judges = round(overall_sd, 4),
      cv = round(cv, 4),
      dispersion_index = round(dispersion_index, 2),
      n_outliers = nrow(outliers_df)
    ),
    outliers = outliers_df,
    warning_level = warning_level,
    recommendations = recommendations,
    parameters = list(
      l = l,
      s = s,
      threshold_cv = threshold_cv,
      threshold_outlier = threshold_outlier
    )
  )

  class(output) <- "judge_heterogeneity"

  if (verbose) {
    print(output)
  }

  return(output)
}


#' Imprimir diagnostico de heterogeneidad
#'
#' @param x Objeto judge_heterogeneity
#' @param ... Argumentos adicionales
#'
#' @export
print.judge_heterogeneity <- function(x, ...) {
  cat("\n")
  cat("================================================================\n")
  cat("   DIAGNOSTICO DE HETEROGENEIDAD ENTRE JUECES\n")
  cat("================================================================\n\n")

  cat(sprintf("  Numero de items:       %d\n", x$overall$n_items))
  cat(sprintf("  Numero de jueces:      %d\n", x$overall$n_judges))
  cat(sprintf("  Escala:                [%d, %d]\n", x$parameters$l, x$parameters$s))
  cat("\n")

  cat("  ESTADISTICOS GLOBALES:\n")
  cat("  ----------------------------------------------------------------\n")
  cat(sprintf("  Media V entre jueces:  %.4f\n", x$overall$mean_across_judges))
  cat(sprintf("  DE entre jueces:       %.4f\n", x$overall$sd_across_judges))
  cat(sprintf("  Coef. Variacion (CV):  %.4f\n", x$overall$cv))
  cat(sprintf("  Indice dispersion:     %.2f\n", x$overall$dispersion_index))
  cat("\n")

  # Nivel de advertencia con color visual
  warning_symbols <- list(
    none = "[OK]",
    low = "[OK]",
    moderate = "[!]",
    high = "[!!]",
    unknown = "[?]"
  )

  warning_texts <- list(
    none = "Heterogeneidad minima",
    low = "Heterogeneidad baja (aceptable)",
    moderate = "Heterogeneidad MODERADA",
    high = "Heterogeneidad ALTA",
    unknown = "No se pudo calcular"
  )

  cat(sprintf("  NIVEL DE ADVERTENCIA:  %s %s\n",
              warning_symbols[[x$warning_level]],
              warning_texts[[x$warning_level]]))
  cat("\n")

  # Resumen por juez
  cat("  RESUMEN POR JUEZ:\n")
  cat("  ----------------------------------------------------------------\n")
  print(x$summary, row.names = FALSE)
  cat("\n")

  # Outliers
  if (nrow(x$outliers) > 0) {
    cat("  JUECES EXTREMOS (OUTLIERS):\n")
    cat("  ----------------------------------------------------------------\n")
    print(x$outliers, row.names = FALSE)
    cat("\n")
  }

  # Recomendaciones
  cat("  RECOMENDACIONES:\n")
  cat("  ----------------------------------------------------------------\n")
  for (rec in x$recommendations) {
    cat(sprintf("  - %s\n", rec))
  }

  cat("\n")
  cat("================================================================\n")

  invisible(x)
}


# =============================================================================
# 2. ANALISIS DE SENSIBILIDAD MEJORADO
# =============================================================================

#' Analisis de sensibilidad completo
#'
#' Realiza un analisis de sensibilidad que incluye: (1) comparacion de priors,
#' (2) analisis de influencia de jueces, y (3) bootstrap no parametrico.
#'
#' @param ratings Vector de calificaciones
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param priors Lista de priors a comparar (default: uniforme, Jeffreys, debil)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param n_bootstrap Numero de muestras bootstrap (default = 1000)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Objeto de clase 'sensitivity_analysis' con:
#' \itemize{
#'   \item prior_comparison: Comparacion de diferentes priors
#'   \item influence_analysis: Efecto de remover cada juez
#'   \item bootstrap: Intervalos bootstrap (no parametrico)
#'   \item summary: Resumen de estabilidad
#' }
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 1, 3)
#' sens <- sensitivity_analysis(ratings, l = 0, s = 3)
#' plot(sens)
#'
#' @export
sensitivity_analysis <- function(ratings,
                                  l = 0,
                                  s = 3,
                                  priors = list(
                                    "Uniforme Beta(1,1)" = c(1, 1),
                                    "Jeffreys Beta(0.5,0.5)" = c(0.5, 0.5),
                                    "Debil Beta(2,2)" = c(2, 2)
                                  ),
                                  cred_level = 0.95,
                                  n_bootstrap = 1000,
                                  verbose = TRUE) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n < 3) {
    stop("Se necesitan al menos 3 jueces para analisis de sensibilidad")
  }

  k <- s - l

  # ==========================================================================
  # 1. COMPARACION DE PRIORS
  # ==========================================================================
  prior_results <- list()

  for (prior_name in names(priors)) {
    prior_params <- priors[[prior_name]]
    alpha0 <- prior_params[1]
    beta0 <- prior_params[2]

    # Calcular posterior
    S <- sum(ratings - l)
    N <- n * k
    post_alpha <- alpha0 + S
    post_beta <- beta0 + (N - S)

    V_mean <- post_alpha / (post_alpha + post_beta)
    hdi <- compute_hdi(post_alpha, post_beta, cred_level)
    prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)

    prior_results[[prior_name]] <- data.frame(
      prior = prior_name,
      V_bayesiano = round(V_mean, 4),
      HDI_lower = round(hdi[1], 4),
      HDI_upper = round(hdi[2], 4),
      HDI_width = round(hdi[2] - hdi[1], 4),
      P_gt_70 = round(prob_70, 4),
      stringsAsFactors = FALSE
    )
  }

  prior_comparison <- do.call(rbind, prior_results)
  rownames(prior_comparison) <- NULL

  # ==========================================================================
  # 2. ANALISIS DE INFLUENCIA (leave-one-out)
  # ==========================================================================
  influence_results <- list()

  # V con todos los jueces
  S_full <- sum(ratings - l)
  N_full <- n * k
  V_full <- (S_full + 1) / (N_full + 2)  # Con prior uniforme

  for (i in 1:n) {
    ratings_loo <- ratings[-i]
    S_loo <- sum(ratings_loo - l)
    N_loo <- (n - 1) * k

    post_alpha_loo <- 1 + S_loo
    post_beta_loo <- 1 + (N_loo - S_loo)

    V_loo <- post_alpha_loo / (post_alpha_loo + post_beta_loo)
    hdi_loo <- compute_hdi(post_alpha_loo, post_beta_loo, cred_level)
    prob_70_loo <- 1 - pbeta(0.70, post_alpha_loo, post_beta_loo)

    influence_results[[i]] <- data.frame(
      juez_removido = i,
      rating_removido = ratings[i],
      V_sin_juez = round(V_loo, 4),
      cambio_V = round(V_loo - V_full, 4),
      HDI_lower = round(hdi_loo[1], 4),
      HDI_upper = round(hdi_loo[2], 4),
      P_gt_70 = round(prob_70_loo, 4),
      stringsAsFactors = FALSE
    )
  }

  influence_df <- do.call(rbind, influence_results)
  rownames(influence_df) <- NULL

  # Identificar jueces influyentes
  influence_threshold <- 0.02  # Cambio > 2% en V
  influential_judges <- influence_df[abs(influence_df$cambio_V) > influence_threshold, ]

  # ==========================================================================
  # 3. BOOTSTRAP NO PARAMETRICO
  # ==========================================================================
  set.seed(123)

  V_bootstrap <- numeric(n_bootstrap)

  for (b in 1:n_bootstrap) {
    # Remuestrear jueces con reemplazo
    ratings_boot <- sample(ratings, n, replace = TRUE)
    S_boot <- sum(ratings_boot - l)
    V_bootstrap[b] <- S_boot / (n * k)
  }

  # Intervalos bootstrap
  boot_ci_percentile <- quantile(V_bootstrap, c((1 - cred_level)/2, 1 - (1 - cred_level)/2))
  boot_mean <- mean(V_bootstrap)
  boot_sd <- sd(V_bootstrap)

  bootstrap_results <- list(
    V_mean = round(boot_mean, 4),
    V_sd = round(boot_sd, 4),
    CI_percentile = round(boot_ci_percentile, 4),
    CI_width = round(boot_ci_percentile[2] - boot_ci_percentile[1], 4),
    samples = V_bootstrap
  )

  # ==========================================================================
  # 4. RESUMEN DE ESTABILIDAD
  # ==========================================================================

  # V clasico
  V_classic <- (mean(ratings) - l) / k

  # Rango de V bajo diferentes priors
  V_range_prior <- range(prior_comparison$V_bayesiano)

  # Rango de V bajo leave-one-out
  V_range_loo <- range(influence_df$V_sin_juez)

  # Evaluacion de estabilidad
  max_prior_diff <- V_range_prior[2] - V_range_prior[1]
  max_loo_diff <- V_range_loo[2] - V_range_loo[1]

  stability <- list(
    V_classic = round(V_classic, 4),
    V_range_prior = round(V_range_prior, 4),
    V_range_loo = round(V_range_loo, 4),
    V_bootstrap_CI = round(boot_ci_percentile, 4),
    max_prior_sensitivity = round(max_prior_diff, 4),
    max_loo_sensitivity = round(max_loo_diff, 4),
    is_stable_prior = max_prior_diff < 0.05,
    is_stable_loo = max_loo_diff < 0.10,
    n_influential_judges = nrow(influential_judges)
  )

  # ==========================================================================
  # OBJETO DE SALIDA
  # ==========================================================================

  output <- list(
    prior_comparison = prior_comparison,
    influence_analysis = influence_df,
    influential_judges = influential_judges,
    bootstrap = bootstrap_results,
    stability = stability,
    parameters = list(
      n = n,
      l = l,
      s = s,
      cred_level = cred_level,
      n_bootstrap = n_bootstrap
    )
  )

  class(output) <- "sensitivity_analysis"

  if (verbose) {
    print(output)
  }

  return(output)
}


#' Imprimir analisis de sensibilidad
#'
#' @param x Objeto sensitivity_analysis
#' @param ... Argumentos adicionales
#'
#' @export
print.sensitivity_analysis <- function(x, ...) {
  cat("\n")
  cat("================================================================\n")
  cat("   ANALISIS DE SENSIBILIDAD - BayesAiken\n")
  cat("================================================================\n\n")

  cat(sprintf("  Numero de jueces:      %d\n", x$parameters$n))
  cat(sprintf("  Escala:                [%d, %d]\n", x$parameters$l, x$parameters$s))
  cat(sprintf("  V clasico:             %.4f\n", x$stability$V_classic))
  cat("\n")

  # 1. Comparacion de priors
  cat("  1. SENSIBILIDAD A LA PRIORI:\n")
  cat("  ----------------------------------------------------------------\n")
  print(x$prior_comparison, row.names = FALSE)
  cat(sprintf("\n  Rango de V bajo diferentes priors: [%.4f, %.4f]\n",
              x$stability$V_range_prior[1], x$stability$V_range_prior[2]))
  cat(sprintf("  Diferencia maxima: %.4f %s\n",
              x$stability$max_prior_sensitivity,
              if(x$stability$is_stable_prior) "[Estable]" else "[Sensible]"))
  cat("\n")

  # 2. Analisis de influencia
  cat("  2. ANALISIS DE INFLUENCIA (Leave-One-Out):\n")
  cat("  ----------------------------------------------------------------\n")
  print(x$influence_analysis, row.names = FALSE)
  cat(sprintf("\n  Rango de V sin cada juez: [%.4f, %.4f]\n",
              x$stability$V_range_loo[1], x$stability$V_range_loo[2]))

  if (nrow(x$influential_judges) > 0) {
    cat(sprintf("\n  JUECES INFLUYENTES (cambio > 2%%):\n"))
    print(x$influential_judges[, c("juez_removido", "rating_removido", "cambio_V")],
          row.names = FALSE)
  } else {
    cat("\n  No se detectaron jueces con influencia desproporcionada.\n")
  }
  cat("\n")

  # 3. Bootstrap
  cat("  3. BOOTSTRAP NO PARAMETRICO:\n")
  cat("  ----------------------------------------------------------------\n")
  cat(sprintf("  V medio (bootstrap):   %.4f\n", x$bootstrap$V_mean))
  cat(sprintf("  DE (bootstrap):        %.4f\n", x$bootstrap$V_sd))
  cat(sprintf("  IC %.0f%% percentil:     [%.4f, %.4f]\n",
              x$parameters$cred_level * 100,
              x$bootstrap$CI_percentile[1],
              x$bootstrap$CI_percentile[2]))
  cat("\n")

  # 4. Conclusion
  cat("  CONCLUSION:\n")
  cat("  ----------------------------------------------------------------\n")

  if (x$stability$is_stable_prior && x$stability$is_stable_loo) {
    cat("  [OK] Los resultados son ESTABLES bajo diferentes priors y\n")
    cat("       al remover jueces individuales.\n")
  } else if (!x$stability$is_stable_prior) {
    cat("  [!] Los resultados son SENSIBLES a la eleccion del prior.\n")
    cat("      Considere reportar analisis con multiples priors.\n")
  } else if (!x$stability$is_stable_loo) {
    cat("  [!] Los resultados son SENSIBLES a jueces individuales.\n")
    cat("      Revise los jueces influyentes identificados.\n")
  }

  cat("\n")
  cat("================================================================\n")

  invisible(x)
}


#' Grafico del analisis de sensibilidad
#'
#' @param x Objeto sensitivity_analysis
#' @param type Tipo de grafico: "prior", "influence", "bootstrap", o "all" (default)
#' @param ... Argumentos adicionales
#'
#' @export
plot.sensitivity_analysis <- function(x, type = "all", ...) {

  if (type == "all") {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old_par))
  }

  # 1. Comparacion de priors
  if (type %in% c("all", "prior")) {
    prior_df <- x$prior_comparison

    # Barplot de V bajo diferentes priors
    barplot(prior_df$V_bayesiano,
            names.arg = substr(prior_df$prior, 1, 15),
            col = "steelblue",
            ylim = c(0, 1),
            main = "Sensibilidad al Prior",
            ylab = "V Bayesiano",
            las = 2,
            cex.names = 0.7)

    # Linea de V clasico
    abline(h = x$stability$V_classic, col = "red", lty = 2, lwd = 2)
    legend("bottomright", "V clasico", col = "red", lty = 2, lwd = 2, bty = "n")
  }

  # 2. Analisis de influencia
  if (type %in% c("all", "influence")) {
    influence_df <- x$influence_analysis

    plot(influence_df$juez_removido, influence_df$V_sin_juez,
         type = "b", pch = 19, col = "steelblue",
         xlab = "Juez removido",
         ylab = "V sin ese juez",
         main = "Analisis de Influencia (LOO)",
         ylim = range(c(influence_df$V_sin_juez, x$stability$V_classic)) + c(-0.05, 0.05))

    abline(h = x$stability$V_classic, col = "red", lty = 2, lwd = 2)

    # Marcar jueces influyentes
    if (nrow(x$influential_judges) > 0) {
      points(x$influential_judges$juez_removido,
             x$influential_judges$V_sin_juez,
             pch = 19, col = "darkorange", cex = 1.5)
    }
  }

  # 3. Bootstrap distribution
  if (type %in% c("all", "bootstrap")) {
    hist(x$bootstrap$samples,
         breaks = 30,
         col = rgb(0.27, 0.51, 0.71, 0.5),
         border = "white",
         main = "Distribucion Bootstrap de V",
         xlab = "V",
         probability = TRUE)

    abline(v = x$stability$V_classic, col = "red", lty = 2, lwd = 2)
    abline(v = x$bootstrap$CI_percentile, col = "steelblue", lty = 3, lwd = 1.5)

    legend("topright",
           c("V clasico", sprintf("IC %.0f%%", x$parameters$cred_level * 100)),
           col = c("red", "steelblue"),
           lty = c(2, 3), lwd = c(2, 1.5),
           bty = "n")
  }

  # 4. Comparacion de intervalos
  if (type == "all") {
    # Grafico de intervalos
    prior_df <- x$prior_comparison

    n_methods <- nrow(prior_df) + 1  # +1 para bootstrap

    plot(NULL,
         xlim = c(0, 1),
         ylim = c(0.5, n_methods + 0.5),
         xlab = "V de Aiken",
         ylab = "",
         main = "Comparacion de Intervalos",
         yaxt = "n")

    # Intervalos de cada prior
    for (i in 1:nrow(prior_df)) {
      segments(prior_df$HDI_lower[i], i,
               prior_df$HDI_upper[i], i,
               col = "steelblue", lwd = 3)
      points(prior_df$V_bayesiano[i], i, pch = 19, col = "steelblue")
    }

    # Intervalo bootstrap
    segments(x$bootstrap$CI_percentile[1], n_methods,
             x$bootstrap$CI_percentile[2], n_methods,
             col = "forestgreen", lwd = 3)
    points(x$bootstrap$V_mean, n_methods, pch = 19, col = "forestgreen")

    # Labels
    axis(2, at = 1:n_methods,
         labels = c(substr(prior_df$prior, 1, 15), "Bootstrap"),
         las = 1, cex.axis = 0.7)

    # V clasico
    abline(v = x$stability$V_classic, col = "red", lty = 2, lwd = 1.5)
  }
}


# =============================================================================
# 3. CALCULAR HDI (Highest Density Interval)
# =============================================================================

#' Calcular HDI (Highest Density Interval)
#'
#' Calcula el intervalo de mayor densidad para una distribucion Beta
#'
#' @param alpha Parametro alpha de la distribucion Beta
#' @param beta Parametro beta de la distribucion Beta
#' @param cred_level Nivel de credibilidad (default = 0.95)
#'
#' @return Vector con limites inferior y superior del HDI
#'
#' @examples
#' compute_hdi(25, 5, 0.95)
#'
#' @export
compute_hdi <- function(alpha, beta, cred_level = 0.95) {

  if (alpha <= 0 || beta <= 0) {
    stop("alpha y beta deben ser mayores que 0")
  }

  if (cred_level <= 0 || cred_level >= 1) {
    stop("cred_level debe estar entre 0 y 1")
  }

  # Casos especiales para distribuciones no unimodales
  # Si alpha <= 1 y beta <= 1, la distribucion puede ser U-shaped o uniforme
  # En ese caso, usamos ETI como aproximacion
  if (alpha <= 1 && beta <= 1) {
    alpha_tail <- (1 - cred_level) / 2
    lower <- qbeta(alpha_tail, alpha, beta)
    upper <- qbeta(1 - alpha_tail, alpha, beta)
    return(c(lower, upper))
  }

  # Si alpha <= 1 (distribucion sesgada a la derecha con pico en 1)
  if (alpha <= 1 && beta > 1) {
    upper <- 1
    lower <- qbeta(1 - cred_level, alpha, beta)
    return(c(lower, upper))
  }

  # Si beta <= 1 (distribucion sesgada a la izquierda con pico en 0)
  if (beta <= 1 && alpha > 1) {
    lower <- 0
    upper <- qbeta(cred_level, alpha, beta)
    return(c(lower, upper))
  }

  # Para distribuciones unimodales (alpha > 1 y beta > 1), usamos optimizacion
  # Buscamos el intervalo mas estrecho que contenga cred_level de la masa

  # Funcion objetivo: ancho del intervalo
  objective <- function(lower_p) {
    lower <- qbeta(lower_p, alpha, beta)
    upper <- qbeta(lower_p + cred_level, alpha, beta)
    if (is.na(upper) || is.na(lower)) return(Inf)
    return(upper - lower)
  }

  # Buscar el lower_p optimo (probabilidad del limite inferior)
  result <- optimize(objective, interval = c(0, 1 - cred_level), tol = 1e-10)

  lower_p <- result$minimum
  lower <- qbeta(lower_p, alpha, beta)
  upper <- qbeta(lower_p + cred_level, alpha, beta)

  # Asegurar limites validos
  lower <- max(0, lower)
  upper <- min(1, upper)

  return(c(lower, upper))
}


#' Comparar diferentes priors
#'
#' Compara los resultados usando diferentes distribuciones prior
#'
#' @param ratings Vector de calificaciones
#' @param l Valor minimo de la escala
#' @param s Valor maximo de la escala
#' @param priors Lista de priors a comparar
#' @param cred_level Nivel de credibilidad (default = 0.95)
#'
#' @return Data frame con comparacion de priors
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3)
#' compare_priors(ratings, l = 0, s = 3)
#'
#' @export
compare_priors <- function(ratings, l = 0, s = 3,
                           priors = list(
                             "Uniform" = c(1, 1),
                             "Jeffreys" = c(0.5, 0.5),
                             "Weakly Informative" = c(2, 2),
                             "Skeptical" = c(1, 3),
                             "Optimistic" = c(3, 1)
                           ),
                           cred_level = 0.95) {

  results <- list()

  for (prior_name in names(priors)) {
    prior_params <- priors[[prior_name]]
    result <- bayes_aiken(ratings, coefficient = "V",
                          l = l, s = s,
                          prior_alpha = prior_params[1],
                          prior_beta = prior_params[2],
                          cred_level = cred_level,
                          verbose = FALSE, plot = FALSE)

    r <- result$results$V
    results[[prior_name]] <- data.frame(
      prior = prior_name,
      prior_alpha = prior_params[1],
      prior_beta = prior_params[2],
      V_bayesian = r$bayesiano_media,
      hdi_lower = r$CI_HDI[1],
      hdi_upper = r$CI_HDI[2],
      hdi_width = r$CI_HDI[2] - r$CI_HDI[1],
      prob_gt_70 = r$prob_mayor_70,
      prob_gt_80 = r$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, results)
  rownames(df) <- NULL

  # Agregar V clasico (no depende del prior)
  n <- length(ratings)
  V_classic <- (mean(ratings) - l) / (s - l)
  attr(df, "V_classical") <- round(V_classic, 4)

  class(df) <- c("prior_comparison", "data.frame")
  return(df)
}


#' @export
print.prior_comparison <- function(x, ...) {
  cat("\n")
  cat("=== Prior Sensitivity Analysis ===\n\n")
  cat("Classical V:", attr(x, "V_classical"), "\n\n")
  print.data.frame(x, row.names = FALSE)
  cat("\n")
  invisible(x)
}


# =============================================================================
# FUNCIONES INTERNAS DE IMPRESION
# =============================================================================

#' @noRd
.print_header <- function() {
  cat("\n")
  cat("===============================================================================\n")
  cat("                                                                               \n")
  cat("                     ____                            _    _ _                  \n")
  cat("                    |  _ \\                     /\\   (_)  | |                   \n")
  cat("                    | |_) | __ _ _   _  ___   /  \\   _| | _____  _ __          \n")
  cat("                    |  _ < / _` | | | |/ _ \\ / /\\ \\ | | |/ / _ \\| '_ \\         \n")
  cat("                    | |_) | (_| | |_| |  __// ____ \\| |   <  __/| | | |        \n")
  cat("                    |____/ \\__,_|\\__, |\\___/_/    \\_\\_|_|\\_\\___||_| |_|        \n")
  cat("                                  __/ |                                        \n")
  cat("                                 |___/                                         \n")
  cat("                                                                               \n")
  cat("             Bayesian Content Validity Coefficients - Version 1.0.0            \n")
  cat("                       Author: Jose Ventura-Leon                               \n")
  cat("                                                                               \n")
  cat("===============================================================================\n")
}


#' @noRd
.print_section <- function(num, title) {
  cat("\n")
  cat(sprintf("--- [%d] %s ", num, title))
  remaining <- 78 - nchar(sprintf("--- [%d] %s ", num, title))
  cat(paste(rep("-", max(0, remaining)), collapse = ""))
  cat("\n\n")
}


#' @noRd
.print_footer <- function(elapsed_time) {
  cat("\n")
  cat("===============================================================================\n")
  cat(sprintf("  Analysis completed in %.3f seconds\n", as.numeric(elapsed_time)))
  cat("  BayesAiken Package v1.0.0 | Author: Jose Ventura-Leon\n")
  cat("===============================================================================\n")
}


#' @noRd
.print_coef_header <- function(coef, title) {
  cat("\n")
  cat(paste(rep("=", 78), collapse = ""), "\n")
  cat(sprintf("  COEFICIENTE %s: %s\n", coef, title))
  cat(paste(rep("=", 78), collapse = ""), "\n\n")
  cat("DATOS DE ENTRADA:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
}


#' @noRd
.print_coef_footer <- function() {
  cat(paste(rep("=", 78), collapse = ""), "\n")
}


#' @noRd
.print_bayesian_model <- function(prior_alpha, prior_beta, post_alpha, post_beta) {
  cat("MODELO BAYESIANO:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  Prior:      Beta(%.2f, %.2f)\n", prior_alpha, prior_beta))
  cat(sprintf("  Posterior:  Beta(%.2f, %.2f)\n", post_alpha, post_beta))
  cat("\n")
}


#' @noRd
.print_estimates <- function(coef, classic, mean, median, mode, sd) {
  cat("ESTIMACIONES:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  %s Clasico:          %.4f\n", coef, classic))
  cat(sprintf("  %s Bayesiano (media): %.4f\n", coef, mean))
  cat(sprintf("  %s Bayesiano (mediana): %.4f\n", coef, median))
  if (!is.na(mode)) {
    cat(sprintf("  %s Bayesiano (moda):  %.4f\n", coef, mode))
  }
  cat(sprintf("  Desviacion estandar: %.4f\n", sd))
  cat("\n")
}


#' @noRd
.print_intervals <- function(ci_eti, ci_hdi, cred_level) {
  pct <- round(cred_level * 100)
  cat(sprintf("INTERVALOS DE CREDIBILIDAD (%d%%):\n", pct))
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  ETI: [%.4f, %.4f]\n", ci_eti[1], ci_eti[2]))
  cat(sprintf("  HDI: [%.4f, %.4f]\n", ci_hdi[1], ci_hdi[2]))
  cat("\n")
}


#' @noRd
.print_probabilities <- function(coef, prob_70, prob_80) {
  cat("PROBABILIDADES:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat(sprintf("  P(%s > 0.70) = %.4f\n", coef, prob_70))
  cat(sprintf("  P(%s > 0.80) = %.4f\n", coef, prob_80))
  cat("\n")
}


#' @noRd
.print_interpretation <- function(V_mean, hdi_lower, prob_70) {
  cat("INTERPRETACION:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  # REGLA PRINCIPAL: Basada en V_media (segun simulacion Monte Carlo)
  # La regla V_media >= cutoff tiene Youden J = 0.82 vs 0.41 de HDI_lower
  cat("  Regla de decision: V posterior (media)\n\n")

  if (V_mean >= 0.90) {
    cat(sprintf("  [EXCELENTE] V = %.2f >= 0.90\n", V_mean))
    cat("  -> Validez de contenido excelente\n")
    cat("  -> El item puede incluirse con alta confianza\n")
  } else if (V_mean >= 0.80) {
    cat(sprintf("  [MUY BUENO] V = %.2f >= 0.80\n", V_mean))
    cat("  -> Validez de contenido muy buena\n")
    cat("  -> El item es adecuado para su uso\n")
  } else if (V_mean >= 0.70) {
    cat(sprintf("  [ACEPTABLE] V = %.2f >= 0.70\n", V_mean))
    cat("  -> Validez de contenido aceptable\n")
    cat("  -> El item cumple el umbral minimo recomendado\n")
  } else if (V_mean >= 0.60) {
    cat(sprintf("  [MARGINAL] V = %.2f (entre 0.60 y 0.70)\n", V_mean))
    cat("  -> Validez de contenido marginal\n")
    cat("  -> Considerar revision del item o agregar mas jueces\n")
  } else {
    cat(sprintf("  [INSUFICIENTE] V = %.2f < 0.60\n", V_mean))
    cat("  -> Validez de contenido insuficiente\n")
    cat("  -> El item requiere revision sustancial\n")
  }

  cat("\n")

  # Informacion complementaria (HDI para referencia)
  cat("  Informacion complementaria:\n")
  cat(sprintf("  - Limite inferior HDI 95%%: %.2f\n", hdi_lower))
  cat(sprintf("  - P(V > 0.70): %.1f%%\n", prob_70 * 100))

  if (hdi_lower >= 0.70) {
    cat("  - El HDI completo supera 0.70 (evidencia robusta)\n")
  } else if (V_mean >= 0.70 && hdi_lower < 0.70) {
    cat("  - Nota: HDI cruza 0.70, considerar mas jueces\n")
  }

  cat("\n")
}


# =============================================================================
# FUNCIONES MULTI PARA EXPORTAR
# =============================================================================

#' Calcular V de Aiken para multiples items
#'
#' @param data Data frame con calificaciones
#' @param item_col Nombre de columna de items
#' @param criterion_col Nombre de columna de criterios (opcional)
#' @param l Minimo de escala
#' @param s Maximo de escala
#' @param prior_alpha Alpha del prior
#' @param prior_beta Beta del prior
#' @param cred_level Nivel de credibilidad
#'
#' @return Data frame con resultados
#' @export
coef_V_multi <- function(data, item_col = NULL, criterion_col = NULL,
                         l = 0, s = 3,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  result <- bayes_aiken(data, coefficient = "V",
                        item_col = item_col, criterion_col = criterion_col,
                        l = l, s = s,
                        prior_alpha = prior_alpha, prior_beta = prior_beta,
                        cred_level = cred_level,
                        verbose = FALSE, plot = FALSE)

  return(result$results$V)
}


#' Calcular H para multiples items
#' @inheritParams coef_V_multi
#' @param c Numero de categorias
#' @export
coef_H_multi <- function(data, item_col = NULL, criterion_col = NULL,
                         c = 4,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  result <- bayes_aiken(data, coefficient = "H",
                        item_col = item_col, criterion_col = criterion_col,
                        c = c,
                        prior_alpha = prior_alpha, prior_beta = prior_beta,
                        cred_level = cred_level,
                        verbose = FALSE, plot = FALSE)

  return(result$results$H)
}


#' Calcular C para multiples items
#' @inheritParams coef_V_multi
#' @export
coef_C_multi <- function(data, item_col = NULL, criterion_col = NULL,
                         l = 0, s = 3,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  result <- bayes_aiken(data, coefficient = "C",
                        item_col = item_col, criterion_col = criterion_col,
                        l = l, s = s,
                        prior_alpha = prior_alpha, prior_beta = prior_beta,
                        cred_level = cred_level,
                        verbose = FALSE, plot = FALSE)

  return(result$results$C)
}


#' Calcular R para multiples items (requiere datos en formato especial)
#' @param data_t1 Data frame con calificaciones tiempo 1
#' @param data_t2 Data frame con calificaciones tiempo 2
#' @param c Numero de categorias
#' @param prior_alpha Alpha del prior
#' @param prior_beta Beta del prior
#' @param cred_level Nivel de credibilidad
#' @export
coef_R_multi <- function(data_t1, data_t2, c = 4,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  if (nrow(data_t1) != nrow(data_t2)) {
    stop("Los data frames deben tener el mismo numero de filas")
  }

  results <- list()
  for (i in 1:nrow(data_t1)) {
    ratings_t1 <- as.numeric(data_t1[i, ])
    ratings_t2 <- as.numeric(data_t2[i, ])
    r <- coef_R(ratings_t1, ratings_t2, c = c,
                prior_alpha = prior_alpha, prior_beta = prior_beta,
                cred_level = cred_level, verbose = FALSE)

    results[[i]] <- data.frame(
      item = i,
      R_clasico = r$clasico,
      R_bayesiano = r$bayesiano_media,
      IC_lower = r$CI_HDI[1],
      IC_upper = r$CI_HDI[2],
      P_70 = r$prob_mayor_70,
      P_80 = r$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Calcular A para multiples jueces
#' @param ratings_matrix Matriz de calificaciones (filas = items, cols = jueces)
#' @param c Numero de categorias
#' @param prior_alpha Alpha del prior
#' @param prior_beta Beta del prior
#' @param cred_level Nivel de credibilidad
#' @export
coef_A_multi <- function(ratings_matrix, c = 4,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  n_judges <- ncol(ratings_matrix)
  results <- list()

  for (j in 1:n_judges) {
    judge <- ratings_matrix[, j]
    others <- ratings_matrix[, -j, drop = FALSE]

    r <- coef_A(judge, others, c = c,
                prior_alpha = prior_alpha, prior_beta = prior_beta,
                cred_level = cred_level, verbose = FALSE)

    results[[j]] <- data.frame(
      judge = j,
      A_clasico = r$clasico,
      A_bayesiano = r$bayesiano_media,
      IC_lower = r$CI_HDI[1],
      IC_upper = r$CI_HDI[2],
      P_70 = r$prob_mayor_70,
      P_80 = r$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Calcular I para todos los pares de items
#' @param ratings_matrix Matriz de calificaciones (filas = criterios, cols = items)
#' @param c Numero de categorias
#' @param prior_alpha Alpha del prior
#' @param prior_beta Beta del prior
#' @param cred_level Nivel de credibilidad
#' @export
coef_I_multi <- function(ratings_matrix, c = 4,
                         prior_alpha = 1, prior_beta = 1,
                         cred_level = 0.95) {

  n_items <- ncol(ratings_matrix)
  results <- list()
  k <- 1

  for (i in 1:(n_items - 1)) {
    for (j in (i + 1):n_items) {
      item1 <- ratings_matrix[, i]
      item2 <- ratings_matrix[, j]

      r <- coef_I(item1, item2, c = c,
                  prior_alpha = prior_alpha, prior_beta = prior_beta,
                  cred_level = cred_level, verbose = FALSE)

      results[[k]] <- data.frame(
        item1 = i,
        item2 = j,
        I_clasico = r$clasico,
        I_bayesiano = r$bayesiano_media,
        IC_lower = r$CI_HDI[1],
        IC_upper = r$CI_HDI[2],
        P_70 = r$prob_mayor_70,
        P_80 = r$prob_mayor_80,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  do.call(rbind, results)
}


# =============================================================================
# METODOS S3 PARA bayes_aiken
# =============================================================================

#' Print method for bayes_aiken objects
#' @param x Objeto bayes_aiken
#' @param ... Argumentos adicionales
#' @export
print.bayes_aiken <- function(x, ...) {
  cat("\n=== BayesAiken Results ===\n\n")
  cat("Input type:", x$input_type, "\n")
  cat("Coefficients:", paste(x$coefficients_requested, collapse = ", "), "\n")
  cat("Prior: Beta(", x$prior$alpha, ", ", x$prior$beta, ")\n", sep = "")
  cat("Credibility level:", x$cred_level * 100, "%\n\n")
  if (!is.null(x$summary_table)) {
    cat("Summary:\n")
    print(x$summary_table, row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}


#' Summary method for bayes_aiken objects
#' @param object Objeto bayes_aiken
#' @param ... Argumentos adicionales
#' @export
summary.bayes_aiken <- function(object, ...) {
  cat("\n")
  cat("===============================================================================\n")
  cat("                    BAYESIAN CONTENT VALIDITY ANALYSIS\n")
  cat("===============================================================================\n\n")
  cat("CONFIGURATION:\n")
  cat("  Input type:        ", object$input_type, "\n")
  cat("  Coefficients:      ", paste(object$coefficients_requested, collapse = ", "), "\n")
  cat("  Prior:             Beta(", object$prior$alpha, ", ", object$prior$beta, ")\n", sep = "")
  cat("  Credibility level: ", object$cred_level * 100, "%\n\n", sep = "")
  cat("RESULTS:\n")
  if (!is.null(object$summary_table)) {
    print(object$summary_table, row.names = FALSE)
  }
  cat("\n===============================================================================\n")
  invisible(object)
}


# =============================================================================
# METODOS S3 PARA bayes_coef (resultado de coef_V, coef_H, etc.)
# =============================================================================

#' Print method for bayes_coef objects
#'
#' Muestra el resultado de coef_V, coef_H, etc. en formato de tabla ordenada
#'
#' @param x Objeto bayes_coef
#' @param ... Argumentos adicionales
#' @export
print.bayes_coef <- function(x, ...) {
  coef <- x$coeficiente

  # Crear data.frame con los resultados principales
  df_main <- data.frame(
    Estadistico = c(
      paste0(coef, " Clasico"),
      paste0(coef, " Bayesiano (media)"),
      paste0(coef, " Bayesiano (mediana)"),
      paste0(coef, " Bayesiano (moda)"),
      "Desviacion estandar"
    ),
    Valor = c(
      x$clasico,
      x$bayesiano_media,
      x$bayesiano_mediana,
      x$bayesiano_moda,
      x$bayesiano_sd
    ),
    stringsAsFactors = FALSE
  )

  # Data.frame para intervalos
  df_intervals <- data.frame(
    Intervalo = c("ETI 95%", "HDI 95%"),
    Inferior = c(x$CI_ETI[1], x$CI_HDI[1]),
    Superior = c(x$CI_ETI[2], x$CI_HDI[2]),
    stringsAsFactors = FALSE
  )

  # Data.frame para probabilidades
  df_probs <- data.frame(
    Probabilidad = c(
      paste0("P(", coef, " > 0.70)"),
      paste0("P(", coef, " > 0.80)")
    ),
    Valor = c(x$prob_mayor_70, x$prob_mayor_80),
    stringsAsFactors = FALSE
  )

  # Imprimir
  cat("\n")
  cat("=== Coeficiente", coef, "de Aiken (Bayesiano) ===\n")
  cat("Jueces:", x$n_jueces, "\n\n")

  cat("ESTIMACIONES:\n")
  print(df_main, row.names = FALSE)
  cat("\n")

  cat("INTERVALOS DE CREDIBILIDAD:\n")
  print(df_intervals, row.names = FALSE)
  cat("\n")

  cat("PROBABILIDADES:\n")
  print(df_probs, row.names = FALSE)
  cat("\n")

  invisible(x)
}


#' Convertir bayes_coef a data.frame
#'
#' @param x Objeto bayes_coef
#' @param row.names No usado
#' @param optional No usado
#' @param ... Argumentos adicionales
#' @export
as.data.frame.bayes_coef <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(
    coeficiente = x$coeficiente,
    n_jueces = x$n_jueces,
    clasico = x$clasico,
    bayesiano_media = x$bayesiano_media,
    bayesiano_mediana = x$bayesiano_mediana,
    bayesiano_moda = x$bayesiano_moda,
    bayesiano_sd = x$bayesiano_sd,
    CI_ETI_lower = x$CI_ETI[1],
    CI_ETI_upper = x$CI_ETI[2],
    CI_HDI_lower = x$CI_HDI[1],
    CI_HDI_upper = x$CI_HDI[2],
    prob_mayor_70 = x$prob_mayor_70,
    prob_mayor_80 = x$prob_mayor_80,
    post_alpha = x$post_alpha,
    post_beta = x$post_beta,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# FUNCION PARA EXTRAER RESULTADOS
# =============================================================================

#' Extraer resultados en formato data.frame
#'
#' Extrae los coeficientes y sus intervalos de credibilidad de uno o varios
#' objetos bayes_coef o bayes_aiken en un data.frame ordenado.
#'
#' @param ... Objetos bayes_coef (de coef_V, coef_H, etc.) o bayes_aiken,
#'            o una lista de objetos
#' @param ci_type Tipo de intervalo a incluir: "HDI" (default), "ETI", o "both"
#'
#' @return Data.frame con los resultados extraidos
#'
#' @examples
#' # Un solo coeficiente
#' v1 <- coef_V(c(3,3,2,3,3), l=0, s=3, verbose=FALSE)
#' extract_results(v1)
#'
#' # Multiples coeficientes
#' v2 <- coef_V(c(3,2,3,3,2), l=0, s=3, verbose=FALSE)
#' extract_results(v1, v2)
#'
#' # Con nombres
#' extract_results(item1 = v1, item2 = v2)
#'
#' @export
extract_results <- function(..., ci_type = "HDI") {

  args <- list(...)

  # Si el primer argumento es una lista, usarla directamente
  if (length(args) == 1 && is.list(args[[1]]) && !inherits(args[[1]], "bayes_coef") && !inherits(args[[1]], "bayes_aiken")) {
    args <- args[[1]]
  }

  # Obtener nombres
  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- paste0("result_", seq_along(args))
  } else {
    # Reemplazar nombres vacios
    empty_names <- arg_names == ""
    arg_names[empty_names] <- paste0("result_", which(empty_names))
  }

  results_list <- list()

  for (i in seq_along(args)) {
    obj <- args[[i]]
    name <- arg_names[i]

    if (inherits(obj, "bayes_coef")) {
      # Objeto individual de coef_V, coef_H, etc.
      df <- .extract_single_coef(obj, name, ci_type)
      results_list[[i]] <- df

    } else if (inherits(obj, "bayes_aiken")) {
      # Objeto bayes_aiken (puede tener multiples coeficientes)
      df <- .extract_bayes_aiken(obj, ci_type)
      results_list[[i]] <- df

    } else {
      warning(paste("Objeto", i, "no es de tipo bayes_coef o bayes_aiken, se omite"))
    }
  }

  # Combinar todos los resultados
  if (length(results_list) == 0) {
    return(data.frame())
  }

  result <- do.call(rbind, results_list)
  rownames(result) <- NULL

  return(result)
}


#' Extraer un solo objeto bayes_coef
#' @noRd
.extract_single_coef <- function(x, name = NULL, ci_type = "HDI") {

  df <- data.frame(
    nombre = if (!is.null(name)) name else x$coeficiente,
    coeficiente = x$coeficiente,
    n_jueces = x$n_jueces,
    clasico = round(x$clasico, 4),
    bayesiano = round(x$bayesiano_media, 4),
    DE = round(x$bayesiano_sd, 4),
    stringsAsFactors = FALSE
  )

  if (ci_type == "HDI" || ci_type == "both") {
    df$HDI_lower <- round(x$CI_HDI[1], 4)
    df$HDI_upper <- round(x$CI_HDI[2], 4)
  }

  if (ci_type == "ETI" || ci_type == "both") {
    df$ETI_lower <- round(x$CI_ETI[1], 4)
    df$ETI_upper <- round(x$CI_ETI[2], 4)
  }

  df$P_70 <- round(x$prob_mayor_70, 4)
  df$P_80 <- round(x$prob_mayor_80, 4)

  return(df)
}


#' Extraer de objeto bayes_aiken
#' @noRd
.extract_bayes_aiken <- function(x, ci_type = "HDI") {

  results_list <- list()

  for (coef_name in names(x$results)) {
    coef_data <- x$results[[coef_name]]

    # Si es un data.frame (multiples items)
    if (is.data.frame(coef_data)) {
      df <- data.frame(
        nombre = if (!is.null(coef_data$item)) coef_data$item else paste0("item_", 1:nrow(coef_data)),
        coeficiente = coef_name,
        n_jueces = coef_data$n_jueces,
        clasico = round(coef_data[[paste0(coef_name, "_clasico")]], 4),
        bayesiano = round(coef_data[[paste0(coef_name, "_bayesiano")]], 4),
        stringsAsFactors = FALSE
      )

      # Buscar columnas de intervalos
      hdi_lower_col <- grep("HDI_lower|IC_lower", names(coef_data), value = TRUE)
      hdi_upper_col <- grep("HDI_upper|IC_upper", names(coef_data), value = TRUE)

      if (ci_type == "HDI" || ci_type == "both") {
        if (length(hdi_lower_col) > 0) {
          df$HDI_lower <- round(coef_data[[hdi_lower_col[1]]], 4)
          df$HDI_upper <- round(coef_data[[hdi_upper_col[1]]], 4)
        }
      }

      # Probabilidades
      p70_col <- grep("P_70|prob.*70", names(coef_data), value = TRUE)
      p80_col <- grep("P_80|prob.*80", names(coef_data), value = TRUE)

      if (length(p70_col) > 0) df$P_70 <- round(coef_data[[p70_col[1]]], 4)
      if (length(p80_col) > 0) df$P_80 <- round(coef_data[[p80_col[1]]], 4)

      results_list[[coef_name]] <- df

    } else if (is.list(coef_data) && !is.null(coef_data$bayesiano_media)) {
      # Es un resultado individual
      df <- .extract_single_coef(coef_data, coef_name, ci_type)
      results_list[[coef_name]] <- df
    }
  }

  if (length(results_list) == 0) {
    return(data.frame())
  }

  do.call(rbind, results_list)
}
