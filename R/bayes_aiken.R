#' =============================================================================
#' BAYES AIKEN - FUNCION PRINCIPAL
#' Calcula todos los coeficientes de validez de contenido con enfoque bayesiano
#'
#' @author Jose Ventura-Leon
#' @version 1.0.0
#' =============================================================================

#' Calcular coeficientes de validez de contenido con intervalos bayesianos
#'
#' Funcion principal que calcula uno o mas coeficientes de validez de contenido
#' (V, H, R, C, A, I) utilizando un enfoque bayesiano con intervalos de credibilidad.
#' Acepta tanto vectores como data.frames.
#'
#' @param data Vector de calificaciones o data.frame con multiples items/criterios
#' @param coefficient Coeficiente(s) a calcular: "V", "H", "R", "C", "A", "I" o "all"
#' @param l Valor minimo de la escala (default = 0)
#' @param s Valor maximo de la escala (default = 3)
#' @param c Numero de categorias en la escala (default = 4, para H, R, A, I)
#' @param ratings_t2 Segundo vector de ratings para coeficiente R (reproducibilidad)
#' @param ratings_others Matriz de ratings de otros jueces para coeficiente A
#' @param ratings_item2 Ratings del segundo item para coeficiente I
#' @param item_col Nombre de columna de items (si data es data.frame)
#' @param criterion_col Nombre de columna de criterios (si data es data.frame)
#' @param prior_alpha Parametro alpha del prior Beta (default = 1, uniforme)
#' @param prior_beta Parametro beta del prior Beta (default = 1, uniforme)
#' @param cred_level Nivel de credibilidad (default = 0.95)
#' @param verbose Mostrar resultados detallados (default = TRUE)
#' @param plot Generar graficos (default = FALSE)
#'
#' @return Objeto de clase 'bayes_aiken' con todos los resultados
#'
#' @details
#' ## Coeficientes disponibles:
#' - **V**: V de Aiken - Validez de contenido

#' - **H**: Homogeneidad entre jueces
#' - **R**: Reproducibilidad temporal (test-retest)
#' - **C**: Consenso (basado en varianza)
#' - **A**: Acuerdo de un juez vs los demas
#' - **I**: Consistencia inter-item
#'
#' ## Modelo Bayesiano:
#' - Prior: Beta(alpha, beta), default uniforme Beta(1,1)
#' - Posterior: Beta(alpha + exitos, beta + fracasos)
#' - Intervalos: HDI (Highest Density Interval) y ETI (Equal-Tailed Interval)
#'
#' ## Tipos de Prior:
#' - Beta(1, 1): Uniforme (no informativo)
#' - Beta(0.5, 0.5): Jeffreys
#' - Beta(2, 2): Debilmente informativo, centrado en 0.5
#'
#' @examples
#' # Ejemplo 1: Vector simple - Coeficiente V
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' result <- bayes_aiken(ratings, coefficient = "V", l = 0, s = 3)
#' print(result)
#' plot(result)
#'
#' # Ejemplo 2: Multiples coeficientes
#' ratings <- c(2, 2, 3, 2, 2)
#' result <- bayes_aiken(ratings, coefficient = c("V", "H", "C"), l = 0, s = 3, c = 4)
#'
#' # Ejemplo 3: Data.frame con multiples items
#' data <- data.frame(
#'   item = rep(1:3, each = 3),
#'   criterio = rep(c("relevancia", "representatividad", "claridad"), 3),
#'   J1 = c(3,3,3, 2,3,3, 3,2,3),
#'   J2 = c(3,3,2, 3,3,3, 3,3,2),
#'   J3 = c(2,3,3, 3,2,3, 2,3,3),
#'   J4 = c(3,3,3, 3,3,2, 3,3,3),
#'   J5 = c(3,2,3, 2,3,3, 3,2,3)
#' )
#' result <- bayes_aiken(data, coefficient = "V",
#'                       item_col = "item", criterion_col = "criterio",
#'                       l = 0, s = 3)
#'
#' @export
bayes_aiken <- function(data,
                         coefficient = "V",
                         l = 0,
                         s = 3,
                         c = 4,
                         ratings_t2 = NULL,
                         ratings_others = NULL,
                         ratings_item2 = NULL,
                         item_col = NULL,
                         criterion_col = NULL,
                         prior_alpha = 1,
                         prior_beta = 1,
                         cred_level = 0.95,
                         verbose = TRUE,
                         plot = FALSE) {

  # Inicio del analisis
  start_time <- Sys.time()

  # Validar coeficientes
  valid_coefficients <- c("V", "H", "R", "C", "A", "I", "all")
  coefficient <- toupper(coefficient)

  if ("ALL" %in% coefficient) {
    coefficient <- c("V", "H", "C")  # Solo los que no requieren datos adicionales
  }

  invalid <- setdiff(coefficient, valid_coefficients)
  if (length(invalid) > 0) {
    stop("Coeficiente(s) no valido(s): ", paste(invalid, collapse = ", "),
         "\nUse: V, H, R, C, A, I o 'all'")
  }

  # Determinar si es data.frame o vector
  is_dataframe <- is.data.frame(data)

  # Inicializar objeto resultado
  results <- list(
    call = match.call(),
    coefficients_requested = coefficient,
    input_type = ifelse(is_dataframe, "data.frame", "vector"),
    prior = list(alpha = prior_alpha, beta = prior_beta),
    cred_level = cred_level,
    results = list(),
    summary_table = NULL,
    verbose_text = NULL,
    plots = list(),
    metadata = list(
      start_time = start_time,
      package_version = "1.0.0"
    )
  )

  # Capturar verbose
  if (verbose) {
    verbose_output <- capture.output({
      .print_header()
      cat("\n")
    })
  } else {
    verbose_output <- character(0)
  }

  # ==========================================================================
  # PROCESAR SEGUN TIPO DE INPUT
  # ==========================================================================

 if (is_dataframe) {
    # --- DATA.FRAME: Multiples items ---
    if (verbose) {
      verbose_output <- c(verbose_output, capture.output({
        .print_section(1, "DATOS DE ENTRADA")
        cat(sprintf("  Tipo de entrada:       data.frame\n"))
        cat(sprintf("  Filas (items):         %d\n", nrow(data)))
        cat(sprintf("  Columnas:              %d\n", ncol(data)))
        if (!is.null(item_col)) cat(sprintf("  Columna de items:      %s\n", item_col))
        if (!is.null(criterion_col)) cat(sprintf("  Columna de criterios:  %s\n", criterion_col))
        cat("\n")
      }))
    }

    # Procesar cada coeficiente solicitado
    for (coef in coefficient) {
      if (verbose) {
        verbose_output <- c(verbose_output, capture.output({
          .print_section(2, paste("COEFICIENTE", coef))
        }))
      }

      result_coef <- switch(coef,
        "V" = .process_V_multi(data, item_col, criterion_col, l, s,
                               prior_alpha, prior_beta, cred_level, verbose),
        "H" = .process_H_multi(data, item_col, criterion_col, c,
                               prior_alpha, prior_beta, cred_level, verbose),
        "C" = .process_C_multi(data, item_col, criterion_col, l, s,
                               prior_alpha, prior_beta, cred_level, verbose),
        NULL
      )

      if (!is.null(result_coef)) {
        results$results[[coef]] <- result_coef$results
        if (verbose) {
          verbose_output <- c(verbose_output, result_coef$verbose)
        }
      }
    }

  } else {
    # --- VECTOR: Un solo item ---
    if (verbose) {
      verbose_output <- c(verbose_output, capture.output({
        .print_section(1, "DATOS DE ENTRADA")
        cat(sprintf("  Tipo de entrada:       vector\n"))
        cat(sprintf("  Calificaciones:        %s\n", paste(data, collapse = ", ")))
        cat(sprintf("  Numero de jueces:      %d\n", length(data)))
        cat(sprintf("  Escala:                [%d, %d]\n", l, s))
        cat("\n")
      }))
    }

    # Procesar cada coeficiente solicitado
    for (coef in coefficient) {
      if (verbose) {
        verbose_output <- c(verbose_output, capture.output({
          .print_section(2, paste("COEFICIENTE", coef))
        }))
      }

      result_coef <- switch(coef,
        "V" = .calc_V_single(data, l, s, prior_alpha, prior_beta, cred_level, verbose),
        "H" = .calc_H_single(data, c, prior_alpha, prior_beta, cred_level, verbose),
        "R" = .calc_R_single(data, ratings_t2, c, prior_alpha, prior_beta, cred_level, verbose),
        "C" = .calc_C_single(data, l, s, prior_alpha, prior_beta, cred_level, verbose),
        "A" = .calc_A_single(data, ratings_others, c, prior_alpha, prior_beta, cred_level, verbose),
        "I" = .calc_I_single(data, ratings_item2, c, prior_alpha, prior_beta, cred_level, verbose),
        NULL
      )

      if (!is.null(result_coef)) {
        results$results[[coef]] <- result_coef$result
        if (verbose) {
          verbose_output <- c(verbose_output, result_coef$verbose)
        }
      }
    }
  }

  # ==========================================================================
  # CREAR TABLA RESUMEN
  # ==========================================================================
  if (verbose) {
    verbose_output <- c(verbose_output, capture.output({
      .print_section(3, "RESUMEN")
    }))
  }

  results$summary_table <- .create_summary_table(results$results, is_dataframe)

  if (verbose) {
    verbose_output <- c(verbose_output, capture.output({
      print(results$summary_table)
      cat("\n")
    }))
  }

  # ==========================================================================
  # GENERAR PLOTS SI SE SOLICITA
  # ==========================================================================
  if (plot && !is_dataframe) {
    for (coef in names(results$results)) {
      results$plots[[coef]] <- .generate_plot(results$results[[coef]], coef, cred_level)
    }
  }

  # ==========================================================================
  # FINALIZAR
  # ==========================================================================
  end_time <- Sys.time()
  results$metadata$end_time <- end_time
  results$metadata$elapsed_time <- difftime(end_time, start_time, units = "secs")

  if (verbose) {
    verbose_output <- c(verbose_output, capture.output({
      .print_footer(results$metadata$elapsed_time)
    }))
    results$verbose_text <- paste(verbose_output, collapse = "\n")
    cat(results$verbose_text)
  }

  # Asignar clase
 class(results) <- "bayes_aiken"

  # Mostrar plots si se solicitan
  if (plot && !is_dataframe && length(results$plots) > 0) {
    if (length(results$plots) > 1) {
      n_plots <- length(results$plots)
      par(mfrow = c(ceiling(n_plots/2), min(n_plots, 2)))
    }
    for (coef in names(results$plots)) {
      .plot_single(results$results[[coef]], coef, cred_level, prior_alpha, prior_beta)
    }
    if (length(results$plots) > 1) par(mfrow = c(1, 1))
  }

  return(results)
}


# =============================================================================
# FUNCIONES INTERNAS PARA PROCESAMIENTO
# =============================================================================

#' @noRd
.process_V_multi <- function(data, item_col, criterion_col, l, s,
                              prior_alpha, prior_beta, cred_level, verbose) {

  # Identificar columnas de jueces
  non_judge_cols <- c(item_col, criterion_col)
  judge_cols <- setdiff(names(data), non_judge_cols)

  results_list <- list()
  verbose_output <- character(0)

  for (i in 1:nrow(data)) {
    ratings <- as.numeric(data[i, judge_cols])
    result <- .calc_V_single(ratings, l, s, prior_alpha, prior_beta, cred_level, FALSE)

    item_id <- if (!is.null(item_col)) data[i, item_col] else i
    crit_id <- if (!is.null(criterion_col)) data[i, criterion_col] else NA

    results_list[[i]] <- data.frame(
      item = item_id,
      criterio = crit_id,
      n_jueces = result$result$n_jueces,
      V_clasico = result$result$clasico,
      V_bayesiano = result$result$bayesiano_media,
      IC_lower = result$result$CI_HDI[1],
      IC_upper = result$result$CI_HDI[2],
      P_70 = result$result$prob_mayor_70,
      P_80 = result$result$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)

  if (verbose) {
    verbose_output <- capture.output({
      print(results_df, row.names = FALSE)
      cat("\n")
    })
  }

  return(list(results = results_df, verbose = verbose_output))
}


#' @noRd
.process_H_multi <- function(data, item_col, criterion_col, c,
                              prior_alpha, prior_beta, cred_level, verbose) {

  non_judge_cols <- c(item_col, criterion_col)
  judge_cols <- setdiff(names(data), non_judge_cols)

  results_list <- list()

  for (i in 1:nrow(data)) {
    ratings <- as.numeric(data[i, judge_cols])
    result <- .calc_H_single(ratings, c, prior_alpha, prior_beta, cred_level, FALSE)

    item_id <- if (!is.null(item_col)) data[i, item_col] else i
    crit_id <- if (!is.null(criterion_col)) data[i, criterion_col] else NA

    results_list[[i]] <- data.frame(
      item = item_id,
      criterio = crit_id,
      n_jueces = result$result$n_jueces,
      S = result$result$S,
      H_clasico = result$result$clasico,
      H_bayesiano = result$result$bayesiano_media,
      IC_lower = result$result$CI_HDI[1],
      IC_upper = result$result$CI_HDI[2],
      P_70 = result$result$prob_mayor_70,
      P_80 = result$result$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)

  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      print(results_df, row.names = FALSE)
      cat("\n")
    })
  }

  return(list(results = results_df, verbose = verbose_output))
}


#' @noRd
.process_C_multi <- function(data, item_col, criterion_col, l, s,
                              prior_alpha, prior_beta, cred_level, verbose) {

  non_judge_cols <- c(item_col, criterion_col)
  judge_cols <- setdiff(names(data), non_judge_cols)

  results_list <- list()

  for (i in 1:nrow(data)) {
    ratings <- as.numeric(data[i, judge_cols])
    result <- .calc_C_single(ratings, l, s, prior_alpha, prior_beta, cred_level, FALSE)

    item_id <- if (!is.null(item_col)) data[i, item_col] else i
    crit_id <- if (!is.null(criterion_col)) data[i, criterion_col] else NA

    results_list[[i]] <- data.frame(
      item = item_id,
      criterio = crit_id,
      n_jueces = result$result$n_jueces,
      varianza = result$result$varianza,
      C_clasico = result$result$clasico,
      C_bayesiano = result$result$bayesiano_media,
      IC_lower = result$result$CI_HDI[1],
      IC_upper = result$result$CI_HDI[2],
      P_70 = result$result$prob_mayor_70,
      P_80 = result$result$prob_mayor_80,
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)

  verbose_output <- character(0)
  if (verbose) {
    verbose_output <- capture.output({
      print(results_df, row.names = FALSE)
      cat("\n")
    })
  }

  return(list(results = results_df, verbose = verbose_output))
}


#' @noRd
.create_summary_table <- function(results, is_dataframe) {
  if (is_dataframe) {
    # Para data.frame, combinar todas las tablas
    tables <- list()
    for (coef in names(results)) {
      if (is.data.frame(results[[coef]])) {
        df <- results[[coef]]
        df$coeficiente <- coef
        tables[[coef]] <- df
      }
    }
    if (length(tables) > 0) {
      return(do.call(rbind, tables))
    }
  } else {
    # Para vector, crear tabla resumen simple
    summary_rows <- list()
    for (coef in names(results)) {
      r <- results[[coef]]
      summary_rows[[coef]] <- data.frame(
        coeficiente = coef,
        clasico = r$clasico,
        bayesiano = r$bayesiano_media,
        IC_lower = r$CI_HDI[1],
        IC_upper = r$CI_HDI[2],
        P_70 = r$prob_mayor_70,
        P_80 = r$prob_mayor_80,
        stringsAsFactors = FALSE
      )
    }
    return(do.call(rbind, summary_rows))
  }
  return(NULL)
}
