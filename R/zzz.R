#' @keywords internal
.onAttach <- function(libname, pkgname) {

  v <- tryCatch(as.character(utils::packageVersion(pkgname)), error = function(e) "?")

  logo <- paste0("
 ======================================================================

  ____                            _     _ _
 | __ )  __ _ _   _  ___  ___   / \\   (_) | _____ _ __
 |  _ \\ / _` | | | |/ _ \\/ __| / _ \\  | | |/ / _ \\ '_ \\
 | |_) | (_| | |_| |  __/\\__ \\/ ___ \\ | |   <  __/ | | |
 |____/ \\__,_|\\__, |\\___||___/_/   \\_\\|_|_|\\_\\___|_| |_|
              |___/

 ======================================================================
   Los seis coeficientes de Aiken con inferencia bayesiana
   Modelo Dirichlet-Multinomial
 ----------------------------------------------------------------------
   Dr. Jose Ventura-Leon                            Version ", v, "

   Empiece por:   aiken_bayes(c(3,3,2,3,3), coef = \"V\")
   Interfaz:      run_aiken_app()
 ======================================================================
")

  packageStartupMessage(logo)
}
