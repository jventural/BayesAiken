# =============================================================================
# PRESENTACION APA DE LA TABLA DE RESULTADOS
# =============================================================================
# La tabla de calculo tiene seis columnas por coeficiente, pensadas para usarse
# desde codigo. Para leerla hacen falta 19 columnas, con el intervalo partido en
# dos y sufijos como "_P_7" que no se entienden sin manual.
#
# Aqui se colapsa cada coeficiente a UNA celda con la notacion de siempre:
#
#     .75 [.58, .93]
#     P = .75
#     clasico .80
#
# De 19 columnas a 5, y en un formato que se copia al manuscrito casi tal cual.
# =============================================================================

#' Formatear un numero acotado en 0-1 al estilo APA
#'
#' Sin cero a la izquierda, como exige APA 7 para valores que no pueden superar
#' la unidad.
#'
#' @param x Vector numerico.
#' @param d Decimales (default 2).
#' @return Vector de caracteres.
#' @examples
#' fmt_apa(c(0.75, -0.04, 1))
#' @export
fmt_apa <- function(x, d = 2) {
  s <- formatC(x, format = "f", digits = d)
  s <- sub("^(-?)0\\.", "\\1.", s)
  s[is.na(x)] <- "\u2014"
  s
}

#' Tabla de resultados en formato APA compacto
#'
#' Convierte la salida de \code{\link{aiken_bayes_table}} en una tabla de pocas
#' columnas, con un coeficiente por celda.
#'
#' @param x Objeto \code{aiken_table}, o su elemento \code{global}.
#' @param umbral Umbral cuya probabilidad se muestra. Si es \code{NULL} se usa
#'   el primero disponible.
#' @param clasico Si \code{TRUE} (default) anade la linea con el valor clasico.
#' @param sep Separador entre lineas dentro de la celda. \code{"\\n"} para
#'   consola y Excel; \code{"<br>"} para HTML.
#' @param digits Decimales (default 2).
#' @param html Si \code{TRUE}, envuelve cada linea en \code{<span>} con clase
#'   (\code{ab-val}, \code{ab-p}, \code{ab-clas}, \code{ab-pg}) para poder darle
#'   estilo. Fuerza \code{sep = "<br>"}.
#' @param pg Si \code{TRUE}, anade a la celda de la V el intervalo frecuentista
#'   de Penfield-Giacobbi, como referencia comparativa. Ver
#'   \code{\link{aiken_V_ic_pg}}.
#'
#' @return Un \code{data.frame} con las columnas identificadoras y una columna
#'   por coeficiente.
#'
#' @examples
#' d <- data.frame(Item = c("I1", "I2"), J1 = c(3, 2), J2 = c(3, 3),
#'                 J3 = c(2, 3), J4 = c(3, 2))
#' r <- aiken_bayes_table(d, seed = 1, B = 500, local = FALSE, verbose = FALSE)
#' aiken_apa(r)
#'
#' @export
aiken_apa <- function(x, umbral = NULL, clasico = TRUE, sep = "\n", digits = 2,
                      html = FALSE, pg = FALSE) {
  if (html) sep <- "<br>"
  esObj <- inherits(x, "aiken_table")
  g     <- if (esObj) x$global else as.data.frame(x)
  id    <- if (esObj && length(x$id_cols)) x$id_cols else
             names(g)[!grepl("_(clasico|li|ls|p[0-9]+)$", names(g)) &
                      !names(g) %in% c("V", "H", "C", "R", "A", "I")]
  coefs <- if (esObj) x$coefs else intersect(c("V", "H", "C"), names(g))

  # umbral: se toma el primero que exista, salvo que se pida otro
  cols_p <- grep("_p[0-9]+$", names(g), value = TRUE)
  if (!length(cols_p)) stop("La tabla no trae columnas de probabilidad.", call. = FALSE)
  suf <- unique(sub("^.*_(p[0-9]+)$", "\\1", cols_p))
  sp  <- if (is.null(umbral)) suf[1] else paste0("p", sub("0\\.", "", format(umbral)))
  if (!sp %in% suf)
    stop("Umbral no disponible. Hay: ", paste(sub("^p", ".", suf), collapse = ", "), call. = FALSE)

  out <- g[, intersect(id, names(g)), drop = FALSE]
  if (!ncol(out)) out <- data.frame(Fila = seq_len(nrow(g)))

  env <- function(txt, clase) if (html) paste0("<span class='", clase, "'>", txt, "</span>") else txt

  for (cf in coefs) {
    if (!cf %in% names(g)) next
    val <- paste0(fmt_apa(g[[cf]], digits),
                  " [", fmt_apa(g[[paste0(cf, "_li")]], digits),
                  ", ", fmt_apa(g[[paste0(cf, "_ls")]], digits), "]")
    p   <- paste0("P = ", fmt_apa(g[[paste0(cf, "_", sp)]], digits))
    celda <- paste0(env(val, "ab-val"), sep, env(p, "ab-p"))
    if (clasico && paste0(cf, "_clasico") %in% names(g))
      celda <- paste0(celda, sep,
                      env(paste0("cl\u00e1sico ", fmt_apa(g[[paste0(cf, "_clasico")]], digits)), "ab-clas"))
    # Comparativa frecuentista, solo para la V y solo si se pide
    if (pg && cf == "V" && all(c("V_pg_li", "V_pg_ls") %in% names(g)))
      celda <- paste0(celda, sep,
                      env(paste0("IC PG [", fmt_apa(g$V_pg_li, digits), ", ",
                                 fmt_apa(g$V_pg_ls, digits), "]"), "ab-pg"))
    out[[cf]] <- celda
  }
  attr(out, "umbral") <- as.numeric(sub("^p", "0.", sp))
  out
}


#' Imprimir la tabla APA en consola
#'
#' Coloca cada celda en tres lineas alineadas, que es lo que hace legible el
#' formato. Un \code{data.frame} normal no lo consigue.
#'
#' @param x Objeto \code{aiken_table} o data.frame de \code{aiken_apa}.
#' @param ... Pasados a \code{\link{aiken_apa}}.
#' @return Se invoca por su efecto.
#' @export
print_apa <- function(x, ...) {
  t <- if (inherits(x, "aiken_table") || "V" %in% names(x)) aiken_apa(x, sep = "\n", ...) else x
  u <- attr(t, "umbral")
  cols  <- names(t)
  filas <- lapply(seq_len(nrow(t)), function(i) lapply(cols, function(cc) strsplit(t[i, cc], "\n", fixed = TRUE)[[1]]))
  anchos <- vapply(seq_along(cols), function(j)
    max(nchar(cols[j]), max(vapply(filas, function(f) max(nchar(f[[j]])), numeric(1)))), numeric(1))

  linea <- paste(mapply(function(a, n) formatC(n, width = a, flag = "-"), anchos, cols), collapse = "   ")
  cat("\n", linea, "\n", strrep("\u2500", nchar(linea)), "\n", sep = "")
  for (f in filas) {
    alto <- max(vapply(f, length, numeric(1)))
    for (k in seq_len(alto)) {
      cat(paste(mapply(function(a, v) formatC(if (k <= length(v)) v[k] else "", width = a, flag = "-"),
                       anchos, f), collapse = "   "), "\n", sep = "")
    }
    cat("\n")
  }
  if (!is.null(u))
    cat("Nota. Media posterior [intervalo de CREDIBILIDAD, no de confianza];\n",
        "P = P(coeficiente > ", fmt_apa(u), " | datos, prior).\n\n", sep = "")
  invisible(t)
}
