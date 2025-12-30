#' =============================================================================
#' FUNCIONES DE EXPORTACION Y CAPTURA - BayesAiken
#' =============================================================================

#' Capturar salida verbose de bayes_aiken
#'
#' Ejecuta bayes_aiken y captura toda la salida verbose como texto
#'
#' @param ... Argumentos para bayes_aiken
#'
#' @return Objeto de clase 'bayes_aiken_capture' con texto capturado
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' output <- bayes_aiken_capture(ratings, coefficient = "V", l = 0, s = 3)
#' print(output)  # Muestra el texto formateado
#'
#' @export
bayes_aiken_capture <- function(...) {

  # Capturar la salida
  captured <- capture.output({
    result <- bayes_aiken(..., verbose = TRUE, plot = FALSE)
  })

  # Crear objeto
  output <- list(
    text = captured,
    result = result
  )

  class(output) <- "bayes_aiken_capture"

  return(output)
}


#' @export
print.bayes_aiken_capture <- function(x, ...) {
  cat(paste(x$text, collapse = "\n"))
  cat("\n")
  invisible(x)
}


#' Exportar resultados de bayes_aiken
#'
#' Exporta los resultados a diferentes formatos
#'
#' @param x Objeto bayes_aiken
#' @param file Ruta del archivo de salida
#' @param format Formato: "txt", "csv", "both" (default = "txt")
#' @param include_verbose Incluir salida verbose en txt (default = TRUE)
#'
#' @return Invisible NULL
#'
#' @examples
#' ratings <- c(3, 3, 2, 3, 3, 2, 3, 3, 2, 3)
#' result <- bayes_aiken(ratings, coefficient = "V", l = 0, s = 3, verbose = FALSE)
#' # bayes_aiken_export(result, "output.txt")
#'
#' @export
bayes_aiken_export <- function(x, file, format = "txt", include_verbose = TRUE) {

  if (!inherits(x, "bayes_aiken")) {
    stop("El objeto debe ser de clase 'bayes_aiken'")
  }

  # Obtener nombre base
  base_name <- tools::file_path_sans_ext(file)

  if (format %in% c("txt", "both")) {
    export_to_txt(x, paste0(base_name, ".txt"), include_verbose)
  }

  if (format %in% c("csv", "both")) {
    export_to_csv(x, paste0(base_name, ".csv"))
  }

  invisible(NULL)
}


#' Exportar a archivo de texto
#'
#' @param x Objeto bayes_aiken
#' @param file Ruta del archivo
#' @param include_verbose Incluir verbose (default = TRUE)
#'
#' @export
export_to_txt <- function(x, file, include_verbose = TRUE) {

  sink(file)
  on.exit(sink())

  cat("===============================================================================\n")
  cat("                    BAYESIAN CONTENT VALIDITY ANALYSIS\n")
  cat("                           BayesAiken Package v1.0.0\n")
  cat("===============================================================================\n\n")

  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Informacion general
  cat("ANALYSIS CONFIGURATION:\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("  Coefficients analyzed:", paste(x$coefficients_requested, collapse = ", "), "\n")
  cat("  Input type:", x$input_type, "\n")
  cat("  Prior: Beta(", x$prior$alpha, ", ", x$prior$beta, ")\n", sep = "")
  cat("  Credibility level:", x$cred_level * 100, "%\n")
  cat("\n")

  # Resultados
  if (!is.null(x$verbose_text) && include_verbose) {
    cat("DETAILED OUTPUT:\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(x$verbose_text)
    cat("\n")
  }

  # Tabla resumen
  if (!is.null(x$summary_table)) {
    cat("SUMMARY TABLE:\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    print(x$summary_table)
    cat("\n")
  }

  cat("===============================================================================\n")
  cat("                          End of Report\n")
  cat("===============================================================================\n")

  message("Exported to: ", file)
  invisible(NULL)
}


#' Exportar a archivo CSV
#'
#' @param x Objeto bayes_aiken
#' @param file Ruta del archivo CSV
#'
#' @export
export_to_csv <- function(x, file) {

  if (!is.null(x$summary_table)) {
    write.csv(x$summary_table, file, row.names = FALSE)
    message("Exported to: ", file)
  } else if (x$input_type == "vector") {
    # Crear data.frame desde resultados
    df_list <- list()
    for (coef in names(x$results)) {
      r <- x$results[[coef]]
      df_list[[coef]] <- data.frame(
        coefficient = coef,
        classical = r$clasico,
        bayesian_mean = r$bayesiano_media,
        bayesian_median = r$bayesiano_mediana,
        bayesian_mode = r$bayesiano_moda,
        bayesian_sd = r$bayesiano_sd,
        ci_eti_lower = r$CI_ETI[1],
        ci_eti_upper = r$CI_ETI[2],
        ci_hdi_lower = r$CI_HDI[1],
        ci_hdi_upper = r$CI_HDI[2],
        prob_gt_70 = r$prob_mayor_70,
        prob_gt_80 = r$prob_mayor_80,
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, df_list)
    write.csv(df, file, row.names = FALSE)
    message("Exported to: ", file)
  } else {
    message("No summary table available to export")
  }

  invisible(NULL)
}


#' Crear tabla resumen formateada
#'
#' @param x Objeto bayes_aiken
#' @param digits Numero de decimales (default = 3)
#'
#' @return Data frame con resumen formateado
#'
#' @export
summary_table <- function(x, digits = 3) {

  if (!inherits(x, "bayes_aiken")) {
    stop("El objeto debe ser de clase 'bayes_aiken'")
  }

  if (!is.null(x$summary_table)) {
    return(round(x$summary_table[sapply(x$summary_table, is.numeric)], digits))
  }

  # Crear desde resultados
  rows <- list()
  for (coef in names(x$results)) {
    r <- x$results[[coef]]
    rows[[coef]] <- data.frame(
      Coefficient = coef,
      Classical = round(r$clasico, digits),
      Bayesian = round(r$bayesiano_media, digits),
      HDI_Lower = round(r$CI_HDI[1], digits),
      HDI_Upper = round(r$CI_HDI[2], digits),
      `P(>0.70)` = round(r$prob_mayor_70, digits),
      `P(>0.80)` = round(r$prob_mayor_80, digits),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}
