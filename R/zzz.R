#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "BayesAiken: Bayesian Aiken's V Coefficient\n",
    "Version 0.1.0 | https://github.com/jventural/BayesAiken\n",
    "Based on: Ventura-Leon (2025) - Dirichlet-Multinomial Framework"
  )
}
