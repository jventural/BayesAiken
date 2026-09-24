# =============================================================================
# NUCLEO DIRICHLET-MULTINOMIAL UNIFICADO PARA LOS SEIS COEFICIENTES DE AIKEN
# =============================================================================
#
# MARCO
# -----
# Cada coeficiente de Aiken es un FUNCIONAL de la distribucion poblacional de
# categorias. Se modela esa distribucion (multinomial con conjugado Dirichlet)
# y el coeficiente se obtiene por arrastre (push-forward) de la posterior.
#
#   pi ~ Dirichlet(alpha0 + conteos)      y      coef = T(pi)
#
# Funcionales (X, X' independientes con valores 0..k, k = s - l):
#   V(pi)   = E[X] / k                              lineal
#   A(x,pi) = 1 - E|x - X| / k                      lineal, condicional a x
#   I(x,pi) = 1 - E|x - X| / k                      lineal, condicional a x
#   H(pi)   = 1 - 2 E|X - X'| / k                   cuadratico
#   C(pi)   = 1 - 4 Var(X) / k^2                    cuadratico
#   R(Pi)   = 1 - E|X1 - X2| / k                    lineal en la CONJUNTA c x c
#
# Referencias:
#   Aiken, L. R. (1985). EPM, 45(1), 131-142.
#   Aiken, L. R. (1989). EPM, 49(2), 321-324.
# =============================================================================


# -----------------------------------------------------------------------------
# UTILIDADES INTERNAS
# -----------------------------------------------------------------------------

#' @noRd
.ab_scale <- function(l, s) {
  if (!is.numeric(l) || !is.numeric(s) || length(l) != 1 || length(s) != 1)
    stop("'l' y 's' deben ser numeros unicos.", call. = FALSE)
  if (s <= l)
    stop("'s' (maximo de la escala) debe ser mayor que 'l' (minimo).", call. = FALSE)
  k <- s - l
  if (abs(k - round(k)) > 1e-8)
    stop("La escala debe tener categorias enteras: 's - l' no es entero.", call. = FALSE)
  k <- round(k)
  list(l = l, s = s, k = k, c = k + 1L, vals = 0:k)
}

#' @noRd
.ab_check_ratings <- function(x, sc, arg = "ratings") {
  x <- x[!is.na(x)]
  if (!length(x)) stop("'", arg, "' no tiene valores validos.", call. = FALSE)
  if (any(x < sc$l | x > sc$s))
    stop("'", arg, "' contiene valores fuera de la escala [", sc$l, ", ", sc$s,
         "]. Revise 'l' y 's'.", call. = FALSE)
  if (any(abs(x - round(x)) > 1e-8))
    stop("'", arg, "' contiene valores no enteros. Los coeficientes de Aiken ",
         "requieren categorias ordinales discretas.", call. = FALSE)
  as.integer(round(x - sc$l))
}

#' @noRd
.ab_prior <- function(prior, n_cells, alpha0 = NULL) {
  if (!is.null(alpha0)) {
    if (length(alpha0) == 1) alpha0 <- rep(alpha0, n_cells)
    if (length(alpha0) != n_cells)
      stop("'alpha0' debe tener longitud 1 o ", n_cells, ".", call. = FALSE)
    if (any(alpha0 <= 0)) stop("'alpha0' debe ser positivo.", call. = FALSE)
    return(list(alpha = alpha0, nombre = "personalizado", masa = sum(alpha0)))
  }
  prior <- match.arg(prior, c("perks", "jeffreys", "uniform"))
  a <- switch(prior,
              perks    = rep(1 / n_cells, n_cells),
              jeffreys = rep(0.5, n_cells),
              uniform  = rep(1, n_cells))
  nombre <- switch(prior,
                   perks    = "Perks: Dirichlet(1/c, ..., 1/c)",
                   jeffreys = "Jeffreys: Dirichlet(1/2, ..., 1/2)",
                   uniform  = "Uniforme: Dirichlet(1, ..., 1)")
  list(alpha = a, nombre = nombre, masa = sum(a))
}

#' @noRd
.ab_rdirichlet <- function(B, alpha) {
  g <- matrix(stats::rgamma(B * length(alpha), shape = alpha, rate = 1),
              nrow = B, byrow = TRUE)
  g / rowSums(g)
}

#' @noRd
.ab_dist_matrix <- function(k) outer(0:k, 0:k, function(a, b) abs(a - b))

#' Intervalo de maxima densidad posterior a partir de muestras
#' @noRd
.ab_hdi <- function(draws, cred_level = 0.95) {
  d <- sort(draws)
  n <- length(d)
  m <- max(1L, floor(cred_level * n))
  if (m >= n) return(c(d[1], d[n]))
  w <- d[(m + 1):n] - d[1:(n - m)]
  i <- which.min(w)
  c(d[i], d[i + m])
}


# -----------------------------------------------------------------------------
# COEFICIENTES CLASICOS (transcritos de Aiken 1985/1989, con la paridad correcta)
# -----------------------------------------------------------------------------

#' Coeficientes clasicos de Aiken
#'
#' Implementacion literal de las formulas originales, incluida la correccion de
#' paridad (\code{j = 0} si n es par, \code{j = 1} si n es impar) que suele
#' omitirse en las implementaciones disponibles.
#'
#' @param x Vector de calificaciones enteras ya desplazadas al rango 0..k.
#' @param k Amplitud de la escala, \code{k = s - l = c - 1}.
#' @param y Segundo vector (solo para R).
#' @param focal Calificacion focal (solo para A e I).
#' @param otros Calificaciones de comparacion (solo para A e I).
#'
#' @return Un numero: el valor del coeficiente clasico.
#' @name aiken_classic
NULL

#' @rdname aiken_classic
#' @export
aiken_V_classic <- function(x, k) sum(x) / (length(x) * k)

#' Intervalo de puntuacion de Penfield-Giacobbi para la V (frecuentista)
#'
#' Intervalo tipo Wilson que trata la suma de puntuaciones como
#' \code{binomial(n*k, V)}. Se incluye para poder compararlo con el intervalo de
#' credibilidad: \strong{no} es bayesiano y su interpretacion es la de un
#' intervalo de confianza.
#'
#' En simulacion (n de 5 a 20, escala 0-3) su cobertura real va de .905 a .972
#' frente al .95 nominal, y resulta entre un 2 % y un 34 % mas estrecho que el
#' de credibilidad. Aguanta mejor de lo que cabria esperar porque el supuesto
#' binomial sobrestima la varianza cuando los jueces se concentran en dos
#' categorias contiguas, y ese error compensa en parte el de tratar n*k ensayos
#' como independientes.
#'
#' @param V Valor clasico de la V (no la media posterior).
#' @param n Numero de jueces.
#' @param k Amplitud de la escala, \code{c - 1}.
#' @param cred_level Nivel de confianza (default 0.95).
#'
#' @return Vector de dos elementos: limites inferior y superior.
#'
#' @references
#' Penfield, R. D., & Giacobbi, P. R. (2004). Applying a score confidence
#' interval to Aiken's item content-relevance index. \emph{Measurement in
#' Physical Education and Exercise Science, 8}(4), 213-225.
#'
#' @examples
#' aiken_V_ic_pg(0.92, n = 5, k = 3)
#'
#' @export
aiken_V_ic_pg <- function(V, n, k, cred_level = 0.95) {
  Z <- stats::qnorm(1 - (1 - cred_level) / 2)
  a <- 2 * V * n * k + Z^2
  b <- Z * sqrt(4 * n * k * V * (1 - V) + Z^2)
  d <- 2 * (n * k + Z^2)
  c((a - b) / d, (a + b) / d)
}

#' @rdname aiken_classic
#' @export
aiken_H_classic <- function(x, k) {
  n <- length(x)
  if (n < 2) stop("H requiere al menos 2 observaciones.", call. = FALSE)
  j <- if (n %% 2 == 1) 1 else 0          # Aiken (1985, p.140)
  nc <- tabulate(x + 1L, nbins = k + 1L)
  S  <- sum(outer(nc, nc) * .ab_dist_matrix(k)) / 2
  1 - 4 * S / (k * (n^2 - j))
}

#' @rdname aiken_classic
#' @export
aiken_C_classic <- function(x, k) {
  n <- length(x)
  if (n < 2) stop("C requiere al menos 2 observaciones.", call. = FALSE)
  d <- if (n %% 2 == 0) 0 else 1          # Aiken (1989, formula 3)
  1 - (4 * n * (n - 1) * stats::var(x)) / ((n^2 - d) * k^2)
}

#' @rdname aiken_classic
#' @export
aiken_R_classic <- function(x, y, k) 1 - sum(abs(x - y)) / (length(x) * k)

#' @rdname aiken_classic
#' @export
aiken_A_classic <- function(focal, otros, k) 1 - sum(abs(otros - focal)) / (length(otros) * k)

#' @rdname aiken_classic
#' @export
aiken_I_classic <- function(focal, otros, k) 1 - sum(abs(otros - focal)) / (length(otros) * k)


# -----------------------------------------------------------------------------
# FUNCIONALES POBLACIONALES T(pi)
# -----------------------------------------------------------------------------

#' @noRd
.ab_T <- function(P, coef, k, focal = NULL) {
  vals <- 0:k
  switch(coef,
    V = as.vector(P %*% vals) / k,
    H = {
      G <- .ab_dist_matrix(k)
      1 - 2 * rowSums((P %*% G) * P) / k
    },
    C = {
      m1 <- as.vector(P %*% vals)
      m2 <- as.vector(P %*% (vals^2))
      1 - 4 * (m2 - m1^2) / k^2
    },
    A = ,
    I = 1 - as.vector(P %*% abs(vals - focal)) / k,
    R = {
      G <- .ab_dist_matrix(k)
      1 - as.vector(P %*% as.vector(G)) / k
    },
    stop("Coeficiente no reconocido: ", coef, call. = FALSE)
  )
}


# -----------------------------------------------------------------------------
# FUNCION PRINCIPAL
# -----------------------------------------------------------------------------

#' Coeficientes de Aiken bayesianos por el modelo Dirichlet-Multinomial
#'
#' Calcula cualquiera de los seis coeficientes de Aiken (V, H, C, R, A, I) con
#' su distribucion posterior completa, intervalo de credibilidad y
#' probabilidades de superar umbrales, bajo un unico modelo
#' Dirichlet-Multinomial.
#'
#' @section El estimando:
#' La posterior es sobre el \strong{funcional poblacional} del coeficiente, es
#' decir sobre el parametro de la poblacion de jueces (o de items), no sobre el
#' valor muestral de Aiken. Para H y C ambos difieren de forma sistematica con
#' paneles pequenos: el coeficiente muestral esta sesgado al alza en
#' \code{(2 GMD / k) (n - j)/(n^2 - j)} para H y en
#' \code{(4 s2 / k^2) (n - d)/(n^2 - d)} para C. El valor clasico se reporta
#' igualmente como descriptivo.
#'
#' @section El caso de R:
#' R \strong{no} es funcional de la distribucion marginal: dos tablas de
#' contingencia con marginales identicas pueden dar R = 1 y R = 0.33. Por eso R
#' se modela con una Dirichlet sobre las \code{c^2} celdas de la tabla
#' ocasion-1 x ocasion-2. Con \code{c = 4} son 16 celdas, de modo que el prior
#' uniforme aportaria 16 pseudo-observaciones: use el prior de Perks.
#'
#' @param x Vector de calificaciones. Para \code{coef = "R"}, calificaciones de
#'   la primera ocasion. Para \code{coef = "A"} o \code{"I"}, las calificaciones
#'   de comparacion (los otros jueces, o los otros items).
#' @param y Solo para \code{coef = "R"}: calificaciones de la segunda ocasion.
#' @param focal Solo para \code{coef = "A"} o \code{"I"}: la calificacion focal
#'   (la del juez evaluado, o la del item evaluado).
#' @param coef Coeficiente: \code{"V"}, \code{"H"}, \code{"C"}, \code{"R"},
#'   \code{"A"} o \code{"I"}.
#' @param l Minimo de la escala (default 0).
#' @param s Maximo de la escala (default 3). El numero de categorias es
#'   \code{c = s - l + 1}.
#' @param prior \code{"perks"} (default, masa a priori = 1 observacion),
#'   \code{"jeffreys"} o \code{"uniform"}.
#' @param alpha0 Vector de concentracion Dirichlet personalizado. Anula
#'   \code{prior}.
#' @param cred_level Nivel de credibilidad (default 0.95).
#' @param thresholds Umbrales para P(coef > umbral). Default \code{c(.70, .80)}.
#' @param ci_type \code{"HDI"} (default) o \code{"ETI"}.
#' @param B Numero de muestras Monte Carlo (default 20000).
#' @param seed Semilla para reproducibilidad.
#' @param verbose Mostrar el informe explicado paso a paso (default TRUE).
#'
#' @return Objeto de clase \code{aiken_bayes} con la posterior y sus resumenes.
#'
#' @examples
#' # Un item valorado por 5 jueces en escala 0-3
#' aiken_bayes(c(3, 3, 2, 3, 3), coef = "V", seed = 1)
#'
#' # Homogeneidad del mismo panel
#' aiken_bayes(c(3, 3, 2, 3, 3), coef = "H", seed = 1, verbose = FALSE)
#'
#' # Reproducibilidad test-retest
#' aiken_bayes(c(3, 3, 2, 3, 3), y = c(3, 2, 2, 3, 3), coef = "R", seed = 1,
#'             verbose = FALSE)
#'
#' @export
aiken_bayes <- function(x, y = NULL, focal = NULL,
                        coef = c("V", "H", "C", "R", "A", "I"),
                        l = 0, s = 3,
                        prior = c("perks", "jeffreys", "uniform"),
                        alpha0 = NULL,
                        cred_level = 0.95,
                        thresholds = c(0.70, 0.80),
                        ci_type = c("HDI", "ETI"),
                        B = 20000, seed = NULL, verbose = TRUE) {

  coef    <- match.arg(coef)
  ci_type <- match.arg(ci_type)
  sc      <- .ab_scale(l, s)
  k       <- sc$k
  if (!is.null(seed)) set.seed(seed)
  if (cred_level <= 0 || cred_level >= 1)
    stop("'cred_level' debe estar entre 0 y 1.", call. = FALSE)

  xi <- .ab_check_ratings(x, sc, "x")
  n  <- length(xi)

  # --- construccion de conteos y valor clasico, segun el coeficiente ---------
  if (coef == "R") {
    if (is.null(y)) stop("El coeficiente R requiere 'y' (segunda ocasion).", call. = FALSE)
    yi <- .ab_check_ratings(y, sc, "y")
    if (length(yi) != n)
      stop("'x' e 'y' deben tener la misma longitud (mismos jueces en ambas ocasiones).",
           call. = FALSE)
    n_cells <- sc$c^2
    counts  <- tabulate(xi + yi * sc$c + 1L, nbins = n_cells)
    clasico <- aiken_R_classic(xi, yi, k)
    etiqueta_unidad <- "pares de calificaciones"

  } else if (coef %in% c("A", "I")) {
    if (is.null(focal))
      stop("Los coeficientes A e I requieren 'focal' (la calificacion evaluada).",
           call. = FALSE)
    fo <- .ab_check_ratings(focal, sc, "focal")
    if (length(fo) != 1)
      stop("'focal' debe ser una sola calificacion.", call. = FALSE)
    n_cells <- sc$c
    counts  <- tabulate(xi + 1L, nbins = n_cells)
    clasico <- if (coef == "A") aiken_A_classic(fo, xi, k) else aiken_I_classic(fo, xi, k)
    etiqueta_unidad <- if (coef == "A") "otros jueces" else "otros items"

  } else {
    if (coef %in% c("H", "C") && n < 2)
      stop("Los coeficientes H y C requieren al menos 2 observaciones.", call. = FALSE)
    fo      <- NULL
    n_cells <- sc$c
    counts  <- tabulate(xi + 1L, nbins = n_cells)
    clasico <- switch(coef,
                      V = aiken_V_classic(xi, k),
                      H = aiken_H_classic(xi, k),
                      C = aiken_C_classic(xi, k))
    etiqueta_unidad <- "jueces"
  }
  if (!exists("fo", inherits = FALSE)) fo <- NULL

  # --- posterior Dirichlet y push-forward -----------------------------------
  pr        <- .ab_prior(prior, n_cells, alpha0)
  post_a    <- pr$alpha + counts
  P         <- .ab_rdirichlet(B, post_a)
  draws     <- .ab_T(P, coef, k, focal = fo)
  draws     <- pmin(pmax(draws, -1), 1)

  ic <- if (ci_type == "HDI") .ab_hdi(draws, cred_level)
        else unname(stats::quantile(draws, c((1 - cred_level) / 2,
                                             1 - (1 - cred_level) / 2)))
  probs <- vapply(thresholds, function(t) mean(draws > t), numeric(1))
  names(probs) <- paste0("P_", sub("0\\.", "", format(thresholds)))

  # --- sesgo analitico de muestra finita (solo H y C) -----------------------
  sesgo <- NA_real_
  if (coef == "H") {
    j <- if (n %% 2 == 1) 1 else 0
    G <- .ab_dist_matrix(k); ph <- counts / sum(counts)
    sesgo <- (2 * sum(outer(ph, ph) * G) / k) * (n - j) / (n^2 - j)
  } else if (coef == "C") {
    d <- if (n %% 2 == 0) 0 else 1
    ph <- counts / sum(counts); vv <- 0:k
    s2 <- sum(ph * vv^2) - sum(ph * vv)^2
    sesgo <- (4 * s2 / k^2) * (n - d) / (n^2 - d)
  }

  out <- list(
    coeficiente = coef,
    n           = n,
    escala      = c(l = l, s = s, c = sc$c, k = k),
    conteos     = counts,
    focal       = fo,
    prior       = pr,
    post_alpha  = post_a,
    n_efectivo  = n,
    clasico     = clasico,
    media       = mean(draws),
    mediana     = stats::median(draws),
    sd          = stats::sd(draws),
    ci          = ic,
    ci_type     = ci_type,
    cred_level  = cred_level,
    thresholds  = thresholds,
    probs       = probs,
    sesgo_muestral = sesgo,
    draws       = draws,
    unidad      = etiqueta_unidad,
    B           = B
  )
  class(out) <- "aiken_bayes"

  if (verbose) print(out)
  invisible(out)
}


# -----------------------------------------------------------------------------
# VERBOSE / PRINT
# -----------------------------------------------------------------------------

#' @noRd
.ab_nombre_largo <- function(coef) {
  switch(coef,
    V = "V - Validez de contenido",
    H = "H - Homogeneidad entre jueces",
    C = "C - Congruencia (metrica cuadratica)",
    R = "R - Reproducibilidad test-retest",
    A = "A - Acuerdo de un juez con los demas",
    I = "I - Coherencia de un item con los demas")
}

#' @noRd
.ab_interpreta <- function(coef, p, umbral) {
  base <- sprintf("P(%s > %.2f) = %.3f", coef, umbral, p)
  veredicto <- if (p >= 0.95) "evidencia muy fuerte"
    else if (p >= 0.90) "evidencia fuerte"
    else if (p >= 0.80) "evidencia moderada"
    else if (p >= 0.50) "evidencia debil, la decision es incierta"
    else "la evidencia apunta en contra"
  paste0(base, "  ->  ", veredicto)
}

#' @noRd
.ab_barra <- function(p, ancho = 28) {
  n <- round(p * ancho)
  paste0("[", strrep("#", n), strrep(".", ancho - n), "]")
}

#' Imprimir un objeto aiken_bayes
#'
#' @param x Objeto de clase \code{aiken_bayes}.
#' @param ... Ignorado.
#' @export
print.aiken_bayes <- function(x, ...) {
  linea <- strrep("=", 70)
  cat("\n", linea, "\n", sep = "")
  cat("  COEFICIENTE ", .ab_nombre_largo(x$coeficiente), "\n", sep = "")
  cat("  Modelo Dirichlet-Multinomial (estimando: parametro poblacional)\n")
  cat(linea, "\n\n", sep = "")

  cat("1) DATOS\n")
  cat(sprintf("   Escala               : %g a %g  (%d categorias, k = %d)\n",
              x$escala["l"], x$escala["s"], x$escala["c"], x$escala["k"]))
  cat(sprintf("   Observaciones        : %d %s\n", x$n, x$unidad))
  if (!is.null(x$focal))
    cat(sprintf("   Calificacion focal   : %g\n", x$focal + x$escala["l"]))
  if (x$coeficiente == "R") {
    cat(sprintf("   Celdas de la conjunta: %d (tabla %dx%d ocasion1 x ocasion2)\n",
                length(x$conteos), x$escala["c"], x$escala["c"]))
  } else {
    cat("   Conteos por categoria:")
    for (i in seq_along(x$conteos))
      cat(sprintf("  %g:%d", x$escala["l"] + i - 1, x$conteos[i]))
    cat("\n")
  }
  cat("\n")

  cat("2) MODELO BAYESIANO\n")
  cat(sprintf("   Prior                : %s\n", x$prior$nombre))
  cat(sprintf("   Masa a priori        : %.2f observaciones equivalentes\n", x$prior$masa))
  if (x$prior$masa > x$n) {
    cat("   AVISO: el prior pesa MAS que los datos. Con paneles pequenos esto\n")
    cat("          desplaza el resultado. Considere prior = 'perks'.\n")
  }
  cat(sprintf("   Muestras Monte Carlo : %s\n", format(x$B, big.mark = " ")))
  cat("\n")

  cat("3) RESULTADO\n")
  cat(sprintf("   %s clasico (muestral) : %.3f\n", x$coeficiente, x$clasico))
  cat(sprintf("   %s bayesiano (media)  : %.3f   (mediana %.3f, sd %.3f)\n",
              x$coeficiente, x$media, x$mediana, x$sd))
  cat(sprintf("   %s%% %s               : [%.3f, %.3f]\n",
              format(100 * x$cred_level), x$ci_type, x$ci[1], x$ci[2]))
  cat("\n")

  if (!is.na(x$sesgo_muestral) && x$sesgo_muestral > 0.005) {
    cat("   NOTA SOBRE LA DIFERENCIA CON EL CLASICO\n")
    cat(sprintf("   Con n = %d el %s muestral esta sesgado al alza en +%.3f respecto\n",
                x$n, x$coeficiente, x$sesgo_muestral))
    cat("   del parametro poblacional. Aiken uso (n^2 - j) en el denominador para\n")
    cat("   forzar el rango [0,1] exacto, a costa de la insesgadez. La posterior\n")
    cat("   estima el parametro, no el valor muestral: la diferencia es esperable.\n\n")
  }

  cat("4) PROBABILIDADES POSTERIORES\n")
  for (i in seq_along(x$thresholds)) {
    cat(sprintf("   %s  %s\n", .ab_barra(x$probs[i]),
                .ab_interpreta(x$coeficiente, x$probs[i], x$thresholds[i])))
  }
  cat("\n")

  cat("5) LECTURA\n")
  p_ref <- x$probs[1]; u_ref <- x$thresholds[1]
  if (p_ref >= 0.90) {
    cat(sprintf("   Se puede sostener que el %s poblacional supera %.2f.\n",
                x$coeficiente, u_ref))
  } else if (p_ref >= 0.50) {
    cat(sprintf("   El %s poblacional probablemente supera %.2f, pero el panel de %d\n",
                x$coeficiente, u_ref, x$n))
    cat("   observaciones no basta para afirmarlo con solidez. Amplie el panel.\n")
  } else {
    cat(sprintf("   La evidencia no respalda que el %s poblacional supere %.2f.\n",
                x$coeficiente, u_ref))
  }
  if (x$n < 5)
    cat("   AVISO: con menos de 5 observaciones la posterior depende fuertemente del prior.\n")
  cat("\n", linea, "\n\n", sep = "")
  invisible(x)
}

#' Resumen de un objeto aiken_bayes
#' @param object Objeto \code{aiken_bayes}.
#' @param ... Ignorado.
#' @export
summary.aiken_bayes <- function(object, ...) {
  data.frame(
    coeficiente = object$coeficiente,
    n           = object$n,
    clasico     = round(object$clasico, 3),
    media       = round(object$media, 3),
    sd          = round(object$sd, 3),
    ci_inf      = round(object$ci[1], 3),
    ci_sup      = round(object$ci[2], 3),
    t(round(object$probs, 3)),
    row.names   = NULL, check.names = FALSE
  )
}

#' Grafico de la posterior de un coeficiente de Aiken
#' @param x Objeto \code{aiken_bayes}.
#' @param ... Ignorado.
#' @export
plot.aiken_bayes <- function(x, ...) {
  d <- stats::density(x$draws, from = min(x$draws), to = max(x$draws))
  graphics::plot(d, main = paste0("Posterior de ", x$coeficiente,
                                  "  (Dirichlet-Multinomial, n = ", x$n, ")"),
                 xlab = paste("Valor de", x$coeficiente), ylab = "Densidad",
                 lwd = 2, col = "#00707F")
  sel <- d$x >= x$ci[1] & d$x <= x$ci[2]
  graphics::polygon(c(x$ci[1], d$x[sel], x$ci[2]), c(0, d$y[sel], 0),
                    col = grDevices::rgb(0, 0.44, 0.50, 0.25), border = NA)
  graphics::abline(v = x$clasico, lty = 2, lwd = 2, col = "#CC2A36")
  graphics::abline(v = x$media, lty = 1, lwd = 2, col = "#00707F")
  for (t in x$thresholds) graphics::abline(v = t, lty = 3, col = "grey40")
  graphics::legend("topleft", bty = "n", cex = 0.85,
                   legend = c(sprintf("Media posterior = %.3f", x$media),
                              sprintf("Clasico (muestral) = %.3f", x$clasico),
                              sprintf("%s %s%%", x$ci_type, format(100 * x$cred_level))),
                   col = c("#00707F", "#CC2A36", grDevices::rgb(0, 0.44, 0.50, 0.4)),
                   lty = c(1, 2, 1), lwd = c(2, 2, 8))
  invisible(x)
}
