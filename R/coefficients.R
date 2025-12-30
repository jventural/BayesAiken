#' =============================================================================
#' COEFICIENTES DE VALIDEZ DE CONTENIDO - FUNCIONES MODULARES
#' Cada coeficiente puede usarse de forma independiente
#' =============================================================================

# =============================================================================
# COEFICIENTE V DE AIKEN
# =============================================================================

#' Calcular Coeficiente V de Aiken con intervalos bayesianos
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente V
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' result <- coef_V(ratings, l = 0, s = 3)
#'
#' @export
coef_V <- function(ratings, l = 0, s = 3,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_V_single(ratings, l, s, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_V")
  return(result$result)
}


#' @noRd
.calc_V_single <- function(ratings, l, s, prior_alpha, prior_beta, cred_level, verbose) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n == 0) stop("No hay calificaciones validas")
  if (any(ratings < l) || any(ratings > s)) stop("Las calificaciones deben estar entre l y s")

  k <- s - l
  x_bar <- mean(ratings)

  # V clasica
  V_classic <- (x_bar - l) / k

  # Modelo bayesiano
  ratings_transformed <- ratings - l
  sum_success <- sum(ratings_transformed)
  total_trials <- n * k

  post_alpha <- prior_alpha + sum_success
  post_beta <- prior_beta + (total_trials - sum_success)

  # Estadisticos
  V_mean <- post_alpha / (post_alpha + post_beta)
  V_median <- qbeta(0.5, post_alpha, post_beta)
  V_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  V_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("V", "V DE AIKEN (Validez de Contenido)")
      cat(sprintf("  Calificaciones:        %s\n", paste(ratings, collapse = ", ")))
      cat(sprintf("  Numero de jueces (n):  %d\n", n))
      cat(sprintf("  Escala:                [%d, %d]\n", l, s))
      cat(sprintf("  Suma (S):              %.0f\n", sum_success))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("V", V_classic, V_mean, V_median, V_mode, V_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("V", prob_70, prob_80)
      .print_interpretation(CI_HDI[1], prob_70)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "V",
    n_jueces = n,
    clasico = round(V_classic, 2),
    bayesiano_media = round(V_mean, 2),
    bayesiano_mediana = round(V_median, 2),
    bayesiano_moda = round(V_mode, 2),
    bayesiano_sd = round(V_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


# =============================================================================
# COEFICIENTE H (Homogeneidad)
# =============================================================================

#' Calcular Coeficiente H (Homogeneidad) con intervalos bayesianos
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param c Numero de categorias en la escala (default = 4)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente H
#'
#' @examples
#' ratings <- c(2, 2, 3, 2, 2)
#' result <- coef_H(ratings, c = 4)
#'
#' @export
coef_H <- function(ratings, c = 4,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_H_single(ratings, c, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_H")
  return(result$result)
}


#' @noRd
.calc_H_single <- function(ratings, c, prior_alpha, prior_beta, cred_level, verbose) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n < 2) stop("Se necesitan al menos 2 jueces")

  # Calcular S (suma de diferencias absolutas entre pares)
  S <- 0
  n_pairs <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S <- S + abs(ratings[i] - ratings[j])
      n_pairs <- n_pairs + 1
    }
  }

  # H clasico
  denominator <- (c - 1) * (n^2 - 1)
  H_classic <- 1 - (4 * S) / denominator

  # Modelo bayesiano
  max_disagreement <- n_pairs * (c - 1)
  agreement <- max_disagreement - S

  post_alpha <- prior_alpha + agreement
  post_beta <- prior_beta + S

  # Estadisticos
  H_mean <- post_alpha / (post_alpha + post_beta)
  H_median <- qbeta(0.5, post_alpha, post_beta)
  H_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  H_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("H", "H (Homogeneidad entre Jueces)")
      cat(sprintf("  Calificaciones:        %s\n", paste(ratings, collapse = ", ")))
      cat(sprintf("  Numero de jueces (n):  %d\n", n))
      cat(sprintf("  Numero de pares:       %d\n", n_pairs))
      cat(sprintf("  Categorias (c):        %d\n", c))
      cat(sprintf("  Suma diferencias (S):  %d\n", S))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("H", H_classic, H_mean, H_median, H_mode, H_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("H", prob_70, prob_80)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "H",
    n_jueces = n,
    n_pares = n_pairs,
    S = S,
    clasico = round(H_classic, 2),
    bayesiano_media = round(H_mean, 2),
    bayesiano_mediana = round(H_median, 2),
    bayesiano_moda = round(H_mode, 2),
    bayesiano_sd = round(H_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


# =============================================================================
# COEFICIENTE R (Reproducibilidad)
# =============================================================================

#' Calcular Coeficiente R (Reproducibilidad) con intervalos bayesianos
#'
#' @param ratings_t1 Vector de calificaciones en tiempo 1
#' @param ratings_t2 Vector de calificaciones en tiempo 2
#' @param c Numero de categorias en la escala (default = 4)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente R
#'
#' @examples
#' t1 <- c(3, 3, 3, 3, 3)
#' t2 <- c(3, 3, 3, 2, 3)
#' result <- coef_R(t1, t2, c = 4)
#'
#' @export
coef_R <- function(ratings_t1, ratings_t2, c = 4,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_R_single(ratings_t1, ratings_t2, c, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_R")
  return(result$result)
}


#' @noRd
.calc_R_single <- function(ratings_t1, ratings_t2, c, prior_alpha, prior_beta, cred_level, verbose) {

  if (is.null(ratings_t2)) {
    stop("Se requiere ratings_t2 para el coeficiente R (Reproducibilidad)")
  }
  if (length(ratings_t1) != length(ratings_t2)) {
    stop("Los vectores t1 y t2 deben tener igual longitud")
  }

  valid <- !is.na(ratings_t1) & !is.na(ratings_t2)
  ratings_t1 <- ratings_t1[valid]
  ratings_t2 <- ratings_t2[valid]
  n <- length(ratings_t1)

  # Calcular S
  S <- sum(abs(ratings_t1 - ratings_t2))

  # R clasico
  denominator <- n * (c - 1)
  R_classic <- 1 - S / denominator

  # Modelo bayesiano
  max_diff <- n * (c - 1)
  agreement <- max_diff - S

  post_alpha <- prior_alpha + agreement
  post_beta <- prior_beta + S

  # Estadisticos
  R_mean <- post_alpha / (post_alpha + post_beta)
  R_median <- qbeta(0.5, post_alpha, post_beta)
  R_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  R_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("R", "R (Reproducibilidad Temporal)")
      cat(sprintf("  Tiempo 1:              %s\n", paste(ratings_t1, collapse = ", ")))
      cat(sprintf("  Tiempo 2:              %s\n", paste(ratings_t2, collapse = ", ")))
      cat(sprintf("  Numero de jueces (n):  %d\n", n))
      cat(sprintf("  Categorias (c):        %d\n", c))
      cat(sprintf("  Suma diferencias (S):  %d\n", S))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("R", R_classic, R_mean, R_median, R_mode, R_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("R", prob_70, prob_80)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "R",
    n_jueces = n,
    S = S,
    clasico = round(R_classic, 2),
    bayesiano_media = round(R_mean, 2),
    bayesiano_mediana = round(R_median, 2),
    bayesiano_moda = round(R_mode, 2),
    bayesiano_sd = round(R_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


# =============================================================================
# COEFICIENTE C (Consenso)
# =============================================================================

#' Calcular Coeficiente C (Consenso) con intervalos bayesianos
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente C
#'
#' @examples
#' ratings <- c(2, 2, 3, 2, 2)
#' result <- coef_C(ratings, l = 0, s = 3)
#'
#' @export
coef_C <- function(ratings, l = 0, s = 3,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_C_single(ratings, l, s, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_C")
  return(result$result)
}


#' @noRd
.calc_C_single <- function(ratings, x_low, x_high, prior_alpha, prior_beta, cred_level, verbose) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n < 2) stop("Se necesitan al menos 2 jueces")

  # Varianza muestral
  S2 <- var(ratings)
  range_scale <- x_high - x_low

  # C clasico
  d <- 1
  numerator <- 4 * n * (n - 1) * S2
  denominator <- (n^2 - d) * (range_scale^2)
  C_classic <- 1 - numerator / denominator
  C_classic <- max(0, min(1, C_classic))

  # Modelo bayesiano
  n_comparisons <- n * (n - 1) / 2
  agreement <- round(C_classic * n_comparisons * range_scale)
  disagreement <- round((1 - C_classic) * n_comparisons * range_scale)
  agreement <- max(0, agreement)
  disagreement <- max(0, disagreement)

  post_alpha <- prior_alpha + agreement
  post_beta <- prior_beta + disagreement

  # Estadisticos
  C_mean <- post_alpha / (post_alpha + post_beta)
  C_median <- qbeta(0.5, post_alpha, post_beta)
  C_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  C_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("C", "C (Consenso)")
      cat(sprintf("  Calificaciones:        %s\n", paste(ratings, collapse = ", ")))
      cat(sprintf("  Numero de jueces (n):  %d\n", n))
      cat(sprintf("  Varianza muestral:     %.3f\n", S2))
      cat(sprintf("  Rango escala:          %d\n", range_scale))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("C", C_classic, C_mean, C_median, C_mode, C_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("C", prob_70, prob_80)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "C",
    n_jueces = n,
    varianza = round(S2, 3),
    clasico = round(C_classic, 2),
    bayesiano_media = round(C_mean, 2),
    bayesiano_mediana = round(C_median, 2),
    bayesiano_moda = round(C_mode, 2),
    bayesiano_sd = round(C_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


# =============================================================================
# COEFICIENTE A (Acuerdo del juez)
# =============================================================================

#' Calcular Coeficiente A (Acuerdo del Juez) con intervalos bayesianos
#'
#' @param ratings_judge Vector de calificaciones del juez focal
#' @param ratings_others Matriz de calificaciones de otros jueces
#' @param c Numero de categorias en la escala (default = 4)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente A
#'
#' @examples
#' judge <- c(2, 3, 3)
#' others <- matrix(c(2,3,2,2, 3,3,3,3, 3,3,3,3), nrow=3, byrow=TRUE)
#' result <- coef_A(judge, others, c = 4)
#'
#' @export
coef_A <- function(ratings_judge, ratings_others, c = 4,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_A_single(ratings_judge, ratings_others, c, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_A")
  return(result$result)
}


#' @noRd
.calc_A_single <- function(ratings_judge, ratings_others, c, prior_alpha, prior_beta, cred_level, verbose) {

  if (is.null(ratings_others)) {
    stop("Se requiere ratings_others para el coeficiente A")
  }

  if (is.vector(ratings_others)) {
    ratings_others <- matrix(ratings_others, nrow = 1)
  }

  n_others <- ncol(ratings_others)
  n_criteria <- length(ratings_judge)

  # Calcular S
  S <- 0
  for (crit in 1:n_criteria) {
    for (j in 1:n_others) {
      S <- S + abs(ratings_judge[crit] - ratings_others[crit, j])
    }
  }

  # A clasico
  n <- n_others + 1
  denominator <- (n - 1) * (c - 1)
  A_classic <- 1 - S / denominator

  # Modelo bayesiano
  max_diff <- n_others * n_criteria * (c - 1)
  agreement <- max_diff - S

  post_alpha <- prior_alpha + agreement
  post_beta <- prior_beta + S

  # Estadisticos
  A_mean <- post_alpha / (post_alpha + post_beta)
  A_median <- qbeta(0.5, post_alpha, post_beta)
  A_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  A_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("A", "A (Acuerdo del Juez)")
      cat(sprintf("  Ratings del juez:      %s\n", paste(ratings_judge, collapse = ", ")))
      cat(sprintf("  Numero otros jueces:   %d\n", n_others))
      cat(sprintf("  Numero de criterios:   %d\n", n_criteria))
      cat(sprintf("  Categorias (c):        %d\n", c))
      cat(sprintf("  Suma diferencias (S):  %d\n", S))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("A", A_classic, A_mean, A_median, A_mode, A_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("A", prob_70, prob_80)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "A",
    n_otros_jueces = n_others,
    n_criterios = n_criteria,
    S = S,
    clasico = round(A_classic, 2),
    bayesiano_media = round(A_mean, 2),
    bayesiano_mediana = round(A_median, 2),
    bayesiano_moda = round(A_mode, 2),
    bayesiano_sd = round(A_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


# =============================================================================
# COEFICIENTE I (Consistencia inter-item)
# =============================================================================

#' Calcular Coeficiente I (Consistencia Inter-Item) con intervalos bayesianos
#'
#' @param ratings_item1 Vector de calificaciones del item 1
#' @param ratings_item2 Vector de calificaciones del item 2
#' @param c Numero de categorias en la escala (default = 4)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1)
#' @param prior_beta Parametro beta del prior Beta (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente I
#'
#' @examples
#' item1 <- c(2, 3, 3)
#' item2 <- c(3, 3, 3)
#' result <- coef_I(item1, item2, c = 4)
#'
#' @export
coef_I <- function(ratings_item1, ratings_item2, c = 4,
                   prior_alpha = 1, prior_beta = 1,
                   cred_level = 0.95, verbose = TRUE) {

  result <- .calc_I_single(ratings_item1, ratings_item2, c, prior_alpha, prior_beta, cred_level, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_I")
  return(result$result)
}


#' @noRd
.calc_I_single <- function(ratings_item1, ratings_item2, c, prior_alpha, prior_beta, cred_level, verbose) {

  if (is.null(ratings_item2)) {
    stop("Se requiere ratings_item2 para el coeficiente I")
  }
  if (length(ratings_item1) != length(ratings_item2)) {
    stop("Los vectores deben tener igual longitud")
  }

  valid <- !is.na(ratings_item1) & !is.na(ratings_item2)
  ratings_item1 <- ratings_item1[valid]
  ratings_item2 <- ratings_item2[valid]

  n_criteria <- length(ratings_item1)
  m <- 2  # Numero de items

  # Calcular S
  S <- sum(abs(ratings_item1 - ratings_item2))

  # I clasico
  denominator <- (m - 1) * (c - 1) * n_criteria
  I_classic <- 1 - S / denominator

  # Modelo bayesiano
  max_diff <- n_criteria * (c - 1)
  agreement <- max_diff - S

  post_alpha <- prior_alpha + agreement
  post_beta <- prior_beta + S

  # Estadisticos
  I_mean <- post_alpha / (post_alpha + post_beta)
  I_median <- qbeta(0.5, post_alpha, post_beta)
  I_mode <- ifelse(post_alpha > 1 & post_beta > 1,
                   (post_alpha - 1) / (post_alpha + post_beta - 2), NA)
  I_sd <- sqrt((post_alpha * post_beta) /
               ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1)))

  # Intervalos
  alpha_cred <- 1 - cred_level
  CI_ETI <- c(qbeta(alpha_cred/2, post_alpha, post_beta),
              qbeta(1 - alpha_cred/2, post_alpha, post_beta))
  CI_HDI <- compute_hdi(post_alpha, post_beta, cred_level)

  # Probabilidades
  prob_70 <- 1 - pbeta(0.70, post_alpha, post_beta)
  prob_80 <- 1 - pbeta(0.80, post_alpha, post_beta)

  # Verbose
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      .print_coef_header("I", "I (Consistencia Inter-Item)")
      cat(sprintf("  Item 1:                %s\n", paste(ratings_item1, collapse = ", ")))
      cat(sprintf("  Item 2:                %s\n", paste(ratings_item2, collapse = ", ")))
      cat(sprintf("  Numero de criterios:   %d\n", n_criteria))
      cat(sprintf("  Categorias (c):        %d\n", c))
      cat(sprintf("  Suma diferencias (S):  %d\n", S))
      cat("\n")
      .print_bayesian_model(prior_alpha, prior_beta, post_alpha, post_beta)
      .print_estimates("I", I_classic, I_mean, I_median, I_mode, I_sd)
      .print_intervals(CI_ETI, CI_HDI, cred_level)
      .print_probabilities("I", prob_70, prob_80)
      .print_coef_footer()
    })
  }

  result <- list(
    coeficiente = "I",
    n_criterios = n_criteria,
    S = S,
    clasico = round(I_classic, 2),
    bayesiano_media = round(I_mean, 2),
    bayesiano_mediana = round(I_median, 2),
    bayesiano_moda = round(I_mode, 2),
    bayesiano_sd = round(I_sd, 2),
    CI_ETI = round(CI_ETI, 2),
    CI_HDI = round(CI_HDI, 2),
    prob_mayor_70 = round(prob_70, 2),
    prob_mayor_80 = round(prob_80, 2),
    post_alpha = post_alpha,
    post_beta = post_beta
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}
