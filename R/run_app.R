# =============================================================================
# LANZADOR DE LA APLICACION SHINY
# =============================================================================

#' Abrir la aplicacion AikenBayes
#'
#' Lanza en local la interfaz grafica para calcular los seis coeficientes de
#' Aiken bajo el modelo Dirichlet-Multinomial. No requiere conexion ni
#' despliegue: la aplicacion corre en la propia sesion de R.
#'
#' @param launch.browser Si \code{TRUE} (default) abre el navegador.
#' @param port Puerto. Si es \code{NULL} lo elige Shiny.
#' @param ... Argumentos adicionales para \code{\link[shiny]{runApp}}.
#'
#' @return Se invoca por su efecto: abre la aplicacion.
#'
#' @examples
#' if (interactive()) {
#'   run_aiken_app()
#' }
#'
#' @export
run_aiken_app <- function(launch.browser = TRUE, port = NULL, ...) {
  faltan <- character(0)
  for (p in c("shiny", "readxl", "writexl")) {
    if (!requireNamespace(p, quietly = TRUE)) faltan <- c(faltan, p)
  }
  if (length(faltan)) {
    stop("Faltan estos paquetes: ", paste(faltan, collapse = ", "),
         "\nInstalelos con:  install.packages(c(\"",
         paste(faltan, collapse = "\", \""), "\"))", call. = FALSE)
  }

  ruta <- system.file("shiny", package = "BayesAiken")
  if (ruta == "" || !file.exists(file.path(ruta, "app.R")))
    stop("No se encontro la aplicacion. Reinstale el paquete BayesAiken.", call. = FALSE)

  message("Abriendo AikenBayes...  (pulse Esc en la consola para cerrarla)")
  shiny::runApp(ruta, launch.browser = launch.browser, port = port, ...)
}
