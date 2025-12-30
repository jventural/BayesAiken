#' =============================================================================
#' UTILIDADES Y HELPERS - BayesAiken
#' =============================================================================

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
.print_interpretation <- function(hdi_lower, prob_70) {
  cat("INTERPRETACION:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")

  # Interpretacion basada en HDI
  if (hdi_lower >= 0.80) {
    cat("  [EXCELENTE] El limite inferior del HDI supera 0.80\n")
    cat("  -> Evidencia fuerte de validez de contenido excelente\n")
  } else if (hdi_lower >= 0.70) {
    cat("  [BUENO] El limite inferior del HDI supera 0.70\n")
    cat("  -> Evidencia de validez de contenido adecuada\n")
  } else if (hdi_lower >= 0.50) {
    cat("  [MODERADO] El limite inferior del HDI esta entre 0.50 y 0.70\n")
    cat("  -> Evidencia moderada, considerar revision\n")
  } else {
    cat("  [BAJO] El limite inferior del HDI es menor a 0.50\n")
    cat("  -> Evidencia debil, se requiere revision\n")
  }

  cat("\n")

  # Interpretacion basada en probabilidad
  if (prob_70 >= 0.95) {
    cat("  La probabilidad de superar 0.70 es mayor al 95%\n")
  } else if (prob_70 >= 0.80) {
    cat("  La probabilidad de superar 0.70 esta entre 80% y 95%\n")
  } else {
    cat("  La probabilidad de superar 0.70 es menor al 80%\n")
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
