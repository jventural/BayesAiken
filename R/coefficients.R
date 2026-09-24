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
      .print_interpretation(V_mean, CI_HDI[1], prob_70)
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


# =============================================================================
# COEFICIENTE V DE AIKEN - MODELO DIRICHLET-MULTINOMIAL
# =============================================================================

#' Calcular Coeficiente V de Aiken con modelo Dirichlet-Multinomial
#'
#' Esta funcion implementa una alternativa al modelo Beta-Binomial estandar,
#' modelando directamente las frecuencias de cada categoria usando una
#' distribucion Dirichlet. Esto respeta mejor la naturaleza categorica de los
#' datos y proporciona intervalos de credibilidad mas honestos.
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Vector de parametros alpha del prior Dirichlet.
#'        Si es un escalar, se usa el mismo valor para todas las categorias.
#'        Default = 1 (prior uniforme)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param n_samples Numero de muestras Monte Carlo (default = 10000)
#' @param seed Semilla para reproducibilidad (default = NULL)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#'
#' @return Lista con resultados del coeficiente V usando modelo Dirichlet
#'
#' @details
#' El modelo Dirichlet-Multinomial:
#' \itemize{
#'   \item Cuenta cuantos jueces eligieron cada categoria (n_0, n_1, ..., n_k)
#'   \item Usa prior Dirichlet(alpha) sobre las probabilidades de categoria
#'   \item La posterior es Dirichlet(alpha + n)
#'   \item V se calcula como suma ponderada: V = sum(pi_c * c/k)
#'   \item Los intervalos se obtienen via simulacion Monte Carlo
#' }
#'
#' Ventajas sobre Beta-Binomial:
#' \itemize{
#'   \item Respeta que los datos son categoricos, no puntos binarios
#'   \item Cada juez contribuye exactamente 1 observacion (sin pseudorreplicacion)
#'   \item Intervalos mas amplios y honestos cuando hay variabilidad
#' }
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' result <- coef_V_dirichlet(ratings, l = 0, s = 3)
#'
#' # Comparar con modelo estandar
#' result_std <- coef_V(ratings, l = 0, s = 3)
#'
#' @export
coef_V_dirichlet <- function(ratings, l = 0, s = 3,
                              prior_alpha = 1,
                              cred_level = 0.95,
                              n_samples = 10000,
                              seed = NULL,
                              verbose = TRUE) {

  result <- .calc_V_dirichlet_single(ratings, l, s, prior_alpha, cred_level,
                                      n_samples, seed, verbose)

  if (verbose) {
    cat(result$verbose)
  }

  class(result$result) <- c("bayes_coef", "bayes_V_dirichlet")
  return(result$result)
}


#' @noRd
.calc_V_dirichlet_single <- function(ratings, l, s, prior_alpha, cred_level,
                                      n_samples, seed, verbose) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n == 0) stop("No hay calificaciones validas")
  if (any(ratings < l) || any(ratings > s)) stop("Las calificaciones deben estar entre l y s")

  k <- s - l  # Numero de categorias - 1
  n_categories <- k + 1

  # Transformar ratings a indices 0:k
  ratings_idx <- ratings - l

  # Contar frecuencias de cada categoria
  counts <- tabulate(ratings_idx + 1, nbins = n_categories)
  names(counts) <- as.character(l:s)

  # Prior Dirichlet
  if (length(prior_alpha) == 1) {
    prior_alpha <- rep(prior_alpha, n_categories)
  }
  if (length(prior_alpha) != n_categories) {
    stop("prior_alpha debe ser escalar o vector de longitud ", n_categories)
  }

  # Posterior Dirichlet: alpha + counts
  post_alpha <- prior_alpha + counts

  # V clasica
  x_bar <- mean(ratings)
  V_classic <- (x_bar - l) / k

  # Simulacion Monte Carlo para obtener distribucion de V
  if (!is.null(seed)) set.seed(seed)

  # Muestrear de Dirichlet usando transformacion Gamma
  # Si X_i ~ Gamma(alpha_i, 1), entonces X_i/sum(X) ~ Dirichlet(alpha)
  V_samples <- numeric(n_samples)

  for (i in 1:n_samples) {
    gamma_samples <- rgamma(n_categories, shape = post_alpha, rate = 1)
    pi_samples <- gamma_samples / sum(gamma_samples)

    # V = sum(pi_c * c/k) donde c va de 0 a k
    category_values <- 0:k
    V_samples[i] <- sum(pi_samples * category_values) / k
  }

  # Estadisticos de la distribucion posterior de V
  V_mean <- mean(V_samples)
  V_median <- median(V_samples)
  V_sd <- sd(V_samples)

  # Calcular moda usando densidad kernel
  dens <- density(V_samples, from = 0, to = 1)
  V_mode <- dens$x[which.max(dens$y)]

  # Intervalos
  alpha_cred <- 1 - cred_level

  # ETI (Equal-Tailed Interval)
  CI_ETI <- quantile(V_samples, c(alpha_cred/2, 1 - alpha_cred/2))

  # HDI (Highest Density Interval)
  CI_HDI <- .compute_hdi_samples(V_samples, cred_level)

  # Probabilidades
  prob_70 <- mean(V_samples > 0.70)
  prob_80 <- mean(V_samples > 0.80)

  # Verbose output
  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      cat("\n")
      cat("===========================================================================\n")
      cat("   COEFICIENTE V DE AIKEN - MODELO DIRICHLET-MULTINOMIAL\n")
      cat("===========================================================================\n")
      cat("\n")
      cat("--- DATOS DE ENTRADA ---\n")
      cat(sprintf("  Calificaciones:        %s\n", paste(ratings, collapse = ", ")))
      cat(sprintf("  Numero de jueces (n):  %d\n", n))
      cat(sprintf("  Escala:                [%d, %d]\n", l, s))
      cat(sprintf("  Categorias:            %d\n", n_categories))
      cat("\n")
      cat("--- FRECUENCIAS POR CATEGORIA ---\n")
      for (cat_i in 1:n_categories) {
        cat(sprintf("  Categoria %d:           %d jueces (%.1f%%)\n",
                    l + cat_i - 1, counts[cat_i], 100 * counts[cat_i] / n))
      }
      cat("\n")
      cat("--- MODELO BAYESIANO ---\n")
      cat(sprintf("  Prior:                 Dirichlet(%s)\n",
                  paste(round(prior_alpha, 2), collapse = ", ")))
      cat(sprintf("  Posterior:             Dirichlet(%s)\n",
                  paste(round(post_alpha, 2), collapse = ", ")))
      cat(sprintf("  Muestras Monte Carlo:  %d\n", n_samples))
      cat("\n")
      cat("--- ESTIMACIONES ---\n")
      cat(sprintf("  V clasica:             %.4f\n", V_classic))
      cat(sprintf("  V posterior (media):   %.4f\n", V_mean))
      cat(sprintf("  V posterior (mediana): %.4f\n", V_median))
      cat(sprintf("  V posterior (moda):    %.4f\n", V_mode))
      cat(sprintf("  Desviacion estandar:   %.4f\n", V_sd))
      cat("\n")
      cat(sprintf("--- INTERVALOS DE CREDIBILIDAD (%.0f%%) ---\n", cred_level * 100))
      cat(sprintf("  ETI:                   [%.4f, %.4f]\n", CI_ETI[1], CI_ETI[2]))
      cat(sprintf("  HDI:                   [%.4f, %.4f]\n", CI_HDI[1], CI_HDI[2]))
      cat(sprintf("  Ancho ETI:             %.4f\n", CI_ETI[2] - CI_ETI[1]))
      cat(sprintf("  Ancho HDI:             %.4f\n", CI_HDI[2] - CI_HDI[1]))
      cat("\n")
      cat("--- PROBABILIDADES POSTERIORES ---\n")
      cat(sprintf("  P(V > 0.70):           %.4f\n", prob_70))
      cat(sprintf("  P(V > 0.80):           %.4f\n", prob_80))
      cat("\n")
      cat("--- INTERPRETACION ---\n")
      cat("  Regla de decision: V posterior (media)\n\n")
      if (V_mean >= 0.90) {
        cat(sprintf("  [EXCELENTE] V = %.2f >= 0.90\n", V_mean))
        cat("  -> Validez de contenido excelente\n")
      } else if (V_mean >= 0.80) {
        cat(sprintf("  [MUY BUENO] V = %.2f >= 0.80\n", V_mean))
        cat("  -> Validez de contenido muy buena\n")
      } else if (V_mean >= 0.70) {
        cat(sprintf("  [ACEPTABLE] V = %.2f >= 0.70\n", V_mean))
        cat("  -> Validez de contenido aceptable\n")
      } else if (V_mean >= 0.60) {
        cat(sprintf("  [MARGINAL] V = %.2f (entre 0.60 y 0.70)\n", V_mean))
        cat("  -> Considerar revision del item\n")
      } else {
        cat(sprintf("  [INSUFICIENTE] V = %.2f < 0.60\n", V_mean))
        cat("  -> Requiere revision sustancial\n")
      }
      cat(sprintf("\n  HDI inferior: %.2f | P(V>0.70): %.1f%%\n", CI_HDI[1], prob_70*100))
      cat("\n")
      cat("--- NOTA SOBRE EL MODELO ---\n")
      cat("  Este modelo trata cada juez como UNA observacion categorica,\n")
      cat("  evitando la pseudorreplicacion del modelo Beta-Binomial.\n")
      cat("  Los intervalos son tipicamente mas amplios y honestos.\n")
      cat("===========================================================================\n")
    })
  }

  result <- list(
    coeficiente = "V_dirichlet",
    modelo = "Dirichlet-Multinomial",
    n_jueces = n,
    n_categories = n_categories,
    counts = counts,
    clasico = round(V_classic, 4),
    bayesiano_media = round(V_mean, 4),
    bayesiano_mediana = round(V_median, 4),
    bayesiano_moda = round(V_mode, 4),
    bayesiano_sd = round(V_sd, 4),
    CI_ETI = round(as.numeric(CI_ETI), 4),
    CI_HDI = round(CI_HDI, 4),
    prob_mayor_70 = round(prob_70, 4),
    prob_mayor_80 = round(prob_80, 4),
    prior_alpha = prior_alpha,
    post_alpha = post_alpha,
    V_samples = V_samples,
    n_samples = n_samples
  )

  return(list(result = result, verbose = paste(verbose_output, collapse = "\n")))
}


#' Calcular HDI desde muestras
#' @noRd
.compute_hdi_samples <- function(samples, cred_level = 0.95) {
  sorted_samples <- sort(samples)
  n <- length(sorted_samples)
  ci_size <- ceiling(cred_level * n)

  # Encontrar el intervalo mas estrecho
  min_width <- Inf
  hdi_lower <- NA
  hdi_upper <- NA

  for (i in 1:(n - ci_size + 1)) {
    width <- sorted_samples[i + ci_size - 1] - sorted_samples[i]
    if (width < min_width) {
      min_width <- width
      hdi_lower <- sorted_samples[i]
      hdi_upper <- sorted_samples[i + ci_size - 1]
    }
  }

  return(c(hdi_lower, hdi_upper))
}


#' Calcular V de Aiken Dirichlet para multiples items
#'
#' @param ratings_matrix Matriz donde cada fila es un item y cada columna un juez
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param prior_alpha Vector de parametros alpha del prior Dirichlet (default = 1)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param n_samples Numero de muestras Monte Carlo (default = 10000)
#' @param seed Semilla para reproducibilidad (default = NULL)
#' @param item_names Nombres de los items (opcional)
#' @param verbose Mostrar resultados (default = TRUE)
#'
#' @return Data frame con resultados para cada item
#'
#' @examples
#' # Matriz de 3 items x 10 jueces
#' ratings <- matrix(c(
#'   3, 3, 2, 3, 3, 2, 3, 3, 2, 3,
#'   2, 2, 2, 3, 2, 2, 2, 3, 2, 2,
#'   3, 3, 3, 3, 3, 3, 3, 3, 3, 2
#' ), nrow = 3, byrow = TRUE)
#' result <- coef_V_dirichlet_multi(ratings, l = 0, s = 3)
#'
#' @export
coef_V_dirichlet_multi <- function(ratings_matrix, l = 0, s = 3,
                                    prior_alpha = 1,
                                    cred_level = 0.95,
                                    n_samples = 10000,
                                    seed = NULL,
                                    item_names = NULL,
                                    verbose = TRUE) {

  if (is.data.frame(ratings_matrix)) {
    ratings_matrix <- as.matrix(ratings_matrix)
  }

  n_items <- nrow(ratings_matrix)

  if (is.null(item_names)) {
    item_names <- paste0("Item_", 1:n_items)
  }

  # Almacenar resultados
  results_list <- vector("list", n_items)
  results_df <- data.frame(
    Item = item_names,
    n_jueces = integer(n_items),
    V_clasica = numeric(n_items),
    V_media = numeric(n_items),
    V_mediana = numeric(n_items),
    V_sd = numeric(n_items),
    ETI_lower = numeric(n_items),
    ETI_upper = numeric(n_items),
    HDI_lower = numeric(n_items),
    HDI_upper = numeric(n_items),
    P_mayor_70 = numeric(n_items),
    P_mayor_80 = numeric(n_items),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n")
    cat("===========================================================================\n")
    cat("   V DE AIKEN DIRICHLET-MULTINOMIAL - ANALISIS MULTIPLE\n")
    cat("===========================================================================\n")
    cat(sprintf("  Numero de items:       %d\n", n_items))
    cat(sprintf("  Escala:                [%d, %d]\n", l, s))
    cat(sprintf("  Nivel de credibilidad: %.0f%%\n", cred_level * 100))
    cat(sprintf("  Muestras Monte Carlo:  %d\n", n_samples))
    cat("===========================================================================\n\n")
  }

  for (i in 1:n_items) {
    ratings_i <- ratings_matrix[i, ]

    result_i <- .calc_V_dirichlet_single(
      ratings = ratings_i,
      l = l, s = s,
      prior_alpha = prior_alpha,
      cred_level = cred_level,
      n_samples = n_samples,
      seed = if (!is.null(seed)) seed + i else NULL,
      verbose = FALSE
    )

    results_list[[i]] <- result_i$result

    results_df$n_jueces[i] <- result_i$result$n_jueces
    results_df$V_clasica[i] <- result_i$result$clasico
    results_df$V_media[i] <- result_i$result$bayesiano_media
    results_df$V_mediana[i] <- result_i$result$bayesiano_mediana
    results_df$V_sd[i] <- result_i$result$bayesiano_sd
    results_df$ETI_lower[i] <- result_i$result$CI_ETI[1]
    results_df$ETI_upper[i] <- result_i$result$CI_ETI[2]
    results_df$HDI_lower[i] <- result_i$result$CI_HDI[1]
    results_df$HDI_upper[i] <- result_i$result$CI_HDI[2]
    results_df$P_mayor_70[i] <- result_i$result$prob_mayor_70
    results_df$P_mayor_80[i] <- result_i$result$prob_mayor_80
  }

  if (verbose) {
    cat("RESULTADOS POR ITEM:\n")
    cat("---------------------------------------------------------------------------\n")
    cat(sprintf("%-12s %5s %7s %7s %7s %12s %12s %7s\n",
                "Item", "n", "V_clas", "V_media", "V_sd", "HDI_95%", "P(V>.70)", "Valid"))
    cat("---------------------------------------------------------------------------\n")

    for (i in 1:n_items) {
      valid_flag <- ifelse(results_df$HDI_lower[i] >= 0.70, "Si",
                           ifelse(results_df$P_mayor_70[i] >= 0.95, "Prob", "No"))
      cat(sprintf("%-12s %5d %7.3f %7.3f %7.3f [%.3f,%.3f] %7.3f   %s\n",
                  results_df$Item[i],
                  results_df$n_jueces[i],
                  results_df$V_clasica[i],
                  results_df$V_media[i],
                  results_df$V_sd[i],
                  results_df$HDI_lower[i],
                  results_df$HDI_upper[i],
                  results_df$P_mayor_70[i],
                  valid_flag))
    }
    cat("---------------------------------------------------------------------------\n")
    cat("\nLeyenda: Valid = 'Si' si HDI_lower >= 0.70, 'Prob' si P(V>.70) >= 0.95\n")
    cat("===========================================================================\n")
  }

  output <- list(
    summary = results_df,
    details = results_list,
    parameters = list(
      l = l,
      s = s,
      prior_alpha = prior_alpha,
      cred_level = cred_level,
      n_samples = n_samples
    )
  )

  class(output) <- c("bayes_V_dirichlet_multi", "list")
  return(output)
}


#' Comparar modelos Beta-Binomial vs Dirichlet-Multinomial
#'
#' Esta funcion calcula V de Aiken usando ambos modelos y compara
#' los intervalos de credibilidad resultantes.
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param n_samples Numero de muestras para Dirichlet (default = 10000)
#' @param seed Semilla para reproducibilidad (default = NULL)
#' @param plot Generar grafico comparativo (default = TRUE)
#'
#' @return Lista con resultados de ambos modelos y comparacion
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' comp <- compare_V_models(ratings, l = 0, s = 3)
#'
#' @export
compare_V_models <- function(ratings, l = 0, s = 3,
                              cred_level = 0.95,
                              n_samples = 10000,
                              seed = NULL,
                              plot = TRUE) {

  # Modelo Beta-Binomial
  result_bb <- .calc_V_single(ratings, l, s,
                               prior_alpha = 1, prior_beta = 1,
                               cred_level = cred_level, verbose = FALSE)

  # Modelo Dirichlet-Multinomial
  result_dir <- .calc_V_dirichlet_single(ratings, l, s,
                                          prior_alpha = 1,
                                          cred_level = cred_level,
                                          n_samples = n_samples,
                                          seed = seed, verbose = FALSE)

  # Comparacion
  comparison <- data.frame(
    Modelo = c("Beta-Binomial", "Dirichlet-Multinomial"),
    V_media = c(result_bb$result$bayesiano_media, result_dir$result$bayesiano_media),
    V_sd = c(result_bb$result$bayesiano_sd, result_dir$result$bayesiano_sd),
    HDI_lower = c(result_bb$result$CI_HDI[1], result_dir$result$CI_HDI[1]),
    HDI_upper = c(result_bb$result$CI_HDI[2], result_dir$result$CI_HDI[2]),
    Ancho_HDI = c(
      result_bb$result$CI_HDI[2] - result_bb$result$CI_HDI[1],
      result_dir$result$CI_HDI[2] - result_dir$result$CI_HDI[1]
    ),
    P_mayor_70 = c(result_bb$result$prob_mayor_70, result_dir$result$prob_mayor_70)
  )

  cat("\n")
  cat("===========================================================================\n")
  cat("   COMPARACION: BETA-BINOMIAL vs DIRICHLET-MULTINOMIAL\n")
  cat("===========================================================================\n")
  cat(sprintf("  Calificaciones: %s\n", paste(ratings, collapse = ", ")))
  cat(sprintf("  n = %d jueces, escala [%d, %d]\n", length(ratings), l, s))
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("%-25s %8s %8s %15s %10s %10s\n",
              "Modelo", "V_media", "V_sd", "HDI_95%", "Ancho", "P(V>.70)"))
  cat("---------------------------------------------------------------------------\n")
  cat(sprintf("%-25s %8.4f %8.4f [%.4f,%.4f] %10.4f %10.4f\n",
              "Beta-Binomial",
              comparison$V_media[1], comparison$V_sd[1],
              comparison$HDI_lower[1], comparison$HDI_upper[1],
              comparison$Ancho_HDI[1], comparison$P_mayor_70[1]))
  cat(sprintf("%-25s %8.4f %8.4f [%.4f,%.4f] %10.4f %10.4f\n",
              "Dirichlet-Multinomial",
              comparison$V_media[2], comparison$V_sd[2],
              comparison$HDI_lower[2], comparison$HDI_upper[2],
              comparison$Ancho_HDI[2], comparison$P_mayor_70[2]))
  cat("---------------------------------------------------------------------------\n")

  ratio_ancho <- comparison$Ancho_HDI[2] / comparison$Ancho_HDI[1]
  cat(sprintf("\n  Ratio de anchos (Dir/BB): %.2f\n", ratio_ancho))

  if (ratio_ancho > 1.1) {
    cat("  El modelo Dirichlet produce intervalos mas amplios (mas conservador).\n")
  } else if (ratio_ancho < 0.9) {
    cat("  El modelo Dirichlet produce intervalos mas estrechos.\n")
  } else {
    cat("  Ambos modelos producen intervalos similares.\n")
  }
  cat("===========================================================================\n")

  # Plot
  if (plot) {
    .plot_model_comparison(result_bb$result, result_dir$result, ratings, l, s, cred_level)
  }

  output <- list(
    beta_binomial = result_bb$result,
    dirichlet = result_dir$result,
    comparison = comparison
  )

  class(output) <- c("bayes_model_comparison", "list")
  invisible(output)
}


#' @noRd
.plot_model_comparison <- function(result_bb, result_dir, ratings, l, s, cred_level) {

  # Generar densidad Beta-Binomial
  x <- seq(0, 1, length.out = 500)
  y_bb <- dbeta(x, result_bb$post_alpha, result_bb$post_beta)

  # Densidad Dirichlet via kernel
  dens_dir <- density(result_dir$V_samples, from = 0, to = 1, n = 500)

  # Normalizar para comparacion visual
  y_dir <- approx(dens_dir$x, dens_dir$y, xout = x)$y
  y_dir[is.na(y_dir)] <- 0

  # Plot
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mar = c(5, 4, 4, 2) + 0.1)

  y_max <- max(c(y_bb, y_dir), na.rm = TRUE) * 1.1

  plot(x, y_bb, type = "l", col = "blue", lwd = 2,
       xlim = c(0, 1), ylim = c(0, y_max),
       xlab = "V de Aiken",
       ylab = "Densidad",
       main = "Comparacion de Modelos: Beta-Binomial vs Dirichlet")

  lines(x, y_dir, col = "red", lwd = 2)

  # HDI regions
  # Beta-Binomial HDI
  idx_bb <- x >= result_bb$CI_HDI[1] & x <= result_bb$CI_HDI[2]
  polygon(c(x[idx_bb][1], x[idx_bb], x[idx_bb][sum(idx_bb)]),
          c(0, y_bb[idx_bb], 0),
          col = rgb(0, 0, 1, 0.2), border = NA)

  # Dirichlet HDI
  idx_dir <- x >= result_dir$CI_HDI[1] & x <= result_dir$CI_HDI[2]
  polygon(c(x[idx_dir][1], x[idx_dir], x[idx_dir][sum(idx_dir)]),
          c(0, y_dir[idx_dir], 0),
          col = rgb(1, 0, 0, 0.2), border = NA)

  # Linea de umbral
  abline(v = 0.70, lty = 2, col = "darkgray", lwd = 1.5)
  text(0.70, y_max * 0.95, "V = 0.70", pos = 4, cex = 0.8)

  # V clasica
  V_classic <- (mean(ratings) - l) / (s - l)
  abline(v = V_classic, lty = 3, col = "darkgreen", lwd = 1.5)
  text(V_classic, y_max * 0.85, sprintf("V clasica = %.3f", V_classic),
       pos = 4, cex = 0.8, col = "darkgreen")

  # Leyenda
  legend("topleft",
         legend = c(
           sprintf("Beta-Binomial (HDI: [%.3f, %.3f])",
                   result_bb$CI_HDI[1], result_bb$CI_HDI[2]),
           sprintf("Dirichlet (HDI: [%.3f, %.3f])",
                   result_dir$CI_HDI[1], result_dir$CI_HDI[2])
         ),
         col = c("blue", "red"),
         lwd = 2,
         fill = c(rgb(0, 0, 1, 0.2), rgb(1, 0, 0, 0.2)),
         border = c("blue", "red"),
         cex = 0.8,
         bty = "n")

  # Info adicional
  mtext(sprintf("n = %d jueces | Escala [%d, %d] | %.0f%% credibilidad",
                length(ratings), l, s, cred_level * 100),
        side = 1, line = 4, cex = 0.8)
}


#' Recomendar modelo para V de Aiken
#'
#' Esta funcion analiza las caracteristicas de los datos y recomienda
#' cual modelo usar: Beta-Binomial (estandar) o Dirichlet-Multinomial
#' (conservador).
#'
#' @param ratings Vector de calificaciones de los jueces
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param threshold_n Umbral de jueces para considerar muestra pequena (default = 7)
#' @param threshold_cv Umbral de CV para considerar alta heterogeneidad (default = 0.15)
#' @param threshold_consensus Umbral de proporcion en categoria modal (default = 0.80)
#' @param verbose Mostrar analisis detallado (default = TRUE)
#'
#' @return Lista con recomendacion y justificacion
#'
#' @details
#' La funcion evalua tres criterios:
#' \itemize{
#'   \item Tamano muestral: n < threshold_n sugiere Dirichlet
#'   \item Heterogeneidad: CV > threshold_cv sugiere Dirichlet
#'   \item Consenso: proporcion modal < threshold_consensus sugiere Dirichlet
#' }
#'
#' Regla de decision:
#' \itemize{
#'   \item 0 banderas rojas: Beta-Binomial
#'   \item 1 bandera roja: Beta-Binomial con precaucion
#'   \item 2+ banderas rojas: Dirichlet-Multinomial
#' }
#'
#' @examples
#' # Caso con alto consenso - recomienda Beta-Binomial
#' ratings_good <- c(3, 3, 3, 3, 3, 3, 3, 3, 2, 3)
#' recommend_model(ratings_good)
#'
#' # Caso con baja n y heterogeneidad - recomienda Dirichlet
#' ratings_risky <- c(1, 2, 3, 2, 3)
#' recommend_model(ratings_risky)
#'
#' @export
recommend_model <- function(ratings, l = 0, s = 3,
                            threshold_n = 7,
                            threshold_cv = 0.15,
                            threshold_consensus = 0.80,
                            verbose = TRUE) {

  ratings <- ratings[!is.na(ratings)]
  n <- length(ratings)

  if (n == 0) stop("No hay calificaciones validas")
  if (any(ratings < l) || any(ratings > s)) stop("Las calificaciones deben estar entre l y s")

  k <- s - l
  n_categories <- k + 1

  # =========================================================================

  # CRITERIO 1: Tamano muestral
  # =========================================================================
  flag_n <- n < threshold_n
  msg_n <- if (flag_n) {
    sprintf("ADVERTENCIA: Muestra pequena (n = %d < %d)", n, threshold_n)
  } else {
    sprintf("OK: Tamano muestral adecuado (n = %d >= %d)", n, threshold_n)
  }

  # =========================================================================
  # CRITERIO 2: Heterogeneidad (CV)
  # =========================================================================
  mean_r <- mean(ratings)
  sd_r <- sd(ratings)
  cv <- if (mean_r > 0) sd_r / mean_r else 0

  flag_cv <- cv > threshold_cv
  msg_cv <- if (flag_cv) {
    sprintf("ADVERTENCIA: Alta heterogeneidad (CV = %.3f > %.2f)", cv, threshold_cv)
  } else {
    sprintf("OK: Heterogeneidad aceptable (CV = %.3f <= %.2f)", cv, threshold_cv)
  }

  # =========================================================================
  # CRITERIO 3: Consenso (proporcion en categoria modal)
  # =========================================================================
  ratings_idx <- ratings - l
  counts <- tabulate(ratings_idx + 1, nbins = n_categories)
  prop_modal <- max(counts) / n
  modal_category <- l + which.max(counts) - 1

  flag_consensus <- prop_modal < threshold_consensus
  msg_consensus <- if (flag_consensus) {
    sprintf("ADVERTENCIA: Bajo consenso (%.1f%% en categoria %d < %.0f%%)",
            prop_modal * 100, modal_category, threshold_consensus * 100)
  } else {
    sprintf("OK: Alto consenso (%.1f%% en categoria %d >= %.0f%%)",
            prop_modal * 100, modal_category, threshold_consensus * 100)
  }

  # =========================================================================
  # DECISION
  # =========================================================================
  n_flags <- sum(c(flag_n, flag_cv, flag_consensus))

  if (n_flags == 0) {
    recommendation <- "Beta-Binomial"
    confidence <- "Alta"
    reason <- "Todos los indicadores son favorables. El modelo estandar es apropiado."
  } else if (n_flags == 1) {
    recommendation <- "Beta-Binomial"
    confidence <- "Moderada"
    reason <- "Un indicador sugiere precaucion. Considere comparar ambos modelos."
  } else {
    recommendation <- "Dirichlet-Multinomial"
    confidence <- "Alta"
    reason <- "Multiples indicadores sugieren que el modelo conservador es mas apropiado."
  }

  # =========================================================================
  # CALCULO COMPARATIVO (para mostrar diferencia)
  # =========================================================================
  # Beta-Binomial
  sum_success <- sum(ratings - l)
  total_trials <- n * k
  post_alpha_bb <- 1 + sum_success
  post_beta_bb <- 1 + (total_trials - sum_success)
  V_mean_bb <- post_alpha_bb / (post_alpha_bb + post_beta_bb)
  V_sd_bb <- sqrt((post_alpha_bb * post_beta_bb) /
                  ((post_alpha_bb + post_beta_bb)^2 * (post_alpha_bb + post_beta_bb + 1)))

  # Dirichlet (aproximacion rapida sin MC completo)
  post_alpha_dir <- 1 + counts
  # Simulacion rapida
  set.seed(123)
  n_sim <- 5000
  V_samples <- numeric(n_sim)
  for (i in 1:n_sim) {
    gamma_samples <- rgamma(n_categories, shape = post_alpha_dir, rate = 1)
    pi_samples <- gamma_samples / sum(gamma_samples)
    V_samples[i] <- sum(pi_samples * (0:k)) / k
  }
  V_mean_dir <- mean(V_samples)
  V_sd_dir <- sd(V_samples)

  ratio_sd <- V_sd_dir / V_sd_bb

  # =========================================================================
  # OUTPUT
  # =========================================================================
  if (verbose) {
    cat("\n")
    cat("===========================================================================\n")
    cat("   RECOMENDACION DE MODELO PARA V DE AIKEN\n")
    cat("===========================================================================\n")
    cat("\n")
    cat("--- DATOS ---\n")
    cat(sprintf("  Calificaciones:        %s\n", paste(ratings, collapse = ", ")))
    cat(sprintf("  Numero de jueces:      %d\n", n))
    cat(sprintf("  Escala:                [%d, %d]\n", l, s))
    cat(sprintf("  V clasica:             %.4f\n", (mean(ratings) - l) / k))
    cat("\n")
    cat("--- EVALUACION DE CRITERIOS ---\n")
    cat(sprintf("  [%s] %s\n", if(flag_n) "X" else "v", msg_n))
    cat(sprintf("  [%s] %s\n", if(flag_cv) "X" else "v", msg_cv))
    cat(sprintf("  [%s] %s\n", if(flag_consensus) "X" else "v", msg_consensus))
    cat("\n")
    cat(sprintf("  Banderas rojas: %d de 3\n", n_flags))
    cat("\n")
    cat("--- COMPARACION RAPIDA ---\n")
    cat(sprintf("  %-25s V = %.4f (SD = %.4f)\n", "Beta-Binomial:", V_mean_bb, V_sd_bb))
    cat(sprintf("  %-25s V = %.4f (SD = %.4f)\n", "Dirichlet-Multinomial:", V_mean_dir, V_sd_dir))
    cat(sprintf("  Ratio SD (Dir/BB):     %.2fx\n", ratio_sd))
    cat("\n")
    cat("---------------------------------------------------------------------------\n")
    cat(sprintf("  RECOMENDACION:         %s\n", recommendation))
    cat(sprintf("  CONFIANZA:             %s\n", confidence))
    cat(sprintf("  RAZON:                 %s\n", reason))
    cat("---------------------------------------------------------------------------\n")
    cat("\n")
    if (n_flags >= 1) {
      cat("  SUGERENCIA: Use compare_V_models() para ver la diferencia en detalle.\n")
    }
    cat("===========================================================================\n")
  }

  result <- list(
    recommendation = recommendation,
    confidence = confidence,
    reason = reason,
    n_flags = n_flags,
    flags = list(
      small_n = flag_n,
      high_heterogeneity = flag_cv,
      low_consensus = flag_consensus
    ),
    diagnostics = list(
      n = n,
      cv = round(cv, 4),
      prop_modal = round(prop_modal, 4),
      modal_category = modal_category
    ),
    comparison = list(
      V_mean_bb = round(V_mean_bb, 4),
      V_sd_bb = round(V_sd_bb, 4),
      V_mean_dir = round(V_mean_dir, 4),
      V_sd_dir = round(V_sd_dir, 4),
      ratio_sd = round(ratio_sd, 2)
    )
  )

  class(result) <- c("model_recommendation", "list")
  invisible(result)
}


#' @export
print.model_recommendation <- function(x, ...) {
  cat("\n")
  cat("Recomendacion de Modelo para V de Aiken\n")
  cat("---------------------------------------\n")
  cat(sprintf("Modelo recomendado: %s\n", x$recommendation))
  cat(sprintf("Confianza:          %s\n", x$confidence))
  cat(sprintf("Banderas rojas:     %d de 3\n", x$n_flags))
  cat(sprintf("Razon:              %s\n", x$reason))
  cat("\n")
  invisible(x)
}
