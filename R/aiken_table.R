# =============================================================================
# ANALISIS COMPLETO DE UNA MATRIZ DE JUECES (formato Item x Criterio x Jueces)
# =============================================================================

#' Detectar las columnas de jueces en un data.frame
#'
#' @param data Data.frame de entrada.
#' @param judge_cols Nombres o indices de las columnas de jueces. Si es
#'   \code{NULL}, se detectan las columnas numericas que no son identificadores.
#' @param id_cols Nombres de las columnas identificadoras (Item, Criterio, ...).
#'
#' @return Lista con \code{id} y \code{judges} (vectores de nombres).
#' @export
detect_judge_cols <- function(data, judge_cols = NULL, id_cols = NULL) {
  nm <- names(data)
  if (!is.null(judge_cols)) {
    if (is.numeric(judge_cols)) judge_cols <- nm[judge_cols]
    faltan <- setdiff(judge_cols, nm)
    if (length(faltan))
      stop("No existen estas columnas de jueces: ", paste(faltan, collapse = ", "),
           call. = FALSE)
    id <- if (is.null(id_cols)) setdiff(nm, judge_cols) else id_cols
    return(list(id = id, judges = judge_cols))
  }
  if (is.null(id_cols)) {
    es_num <- vapply(data, function(cc) is.numeric(cc) || is.integer(cc), logical(1))
    id <- nm[!es_num]
    # una columna numerica con todos los valores distintos suele ser un ID
    for (v in nm[es_num]) {
      if (grepl("^(item|id|n|num|orden)$", tolower(v))) id <- c(id, v)
    }
  } else id <- id_cols
  judges <- setdiff(nm, id)
  if (!length(judges))
    stop("No se identificaron columnas de jueces. Use 'judge_cols'.", call. = FALSE)
  list(id = id, judges = judges)
}


#' Analisis bayesiano completo de una tabla de calificaciones de jueces
#'
#' Calcula V, H y C por fila (item o combinacion item-criterio) y, opcionalmente,
#' A por juez y I por item, todo bajo el modelo Dirichlet-Multinomial.
#'
#' @param data Data.frame con columnas identificadoras y columnas de jueces.
#' @param judge_cols Columnas de jueces. Si es \code{NULL} se detectan solas.
#' @param id_cols Columnas identificadoras. Si es \code{NULL} se detectan solas.
#' @param coefs Coeficientes a calcular por fila. Default \code{c("V","H","C")}.
#' @param local Si \code{TRUE} (default) calcula ademas A por juez-fila e I por
#'   item-juez.
#' @param item_col Nombre de la columna de item, necesaria para I cuando hay
#'   varias filas por item. Si es \code{NULL} se usa la primera de \code{id_cols}.
#' @param l,s Minimo y maximo de la escala.
#' @param prior,alpha0,cred_level,thresholds,ci_type,B Ver \code{\link{aiken_bayes}}.
#' @param seed Semilla para reproducibilidad.
#' @param verbose Mostrar el progreso y el resumen final.
#'
#' @return Objeto de clase \code{aiken_table} con los elementos \code{global}
#'   (una fila por unidad evaluada), \code{acuerdo} (coeficiente A) y
#'   \code{coherencia} (coeficiente I).
#'
#' @examples
#' d <- data.frame(
#'   Item = paste0("It", rep(1:3, each = 2)),
#'   Criterio = rep(c("Relevancia", "Claridad"), 3),
#'   J1 = c(3,3,2,3,3,3), J2 = c(3,2,3,3,2,3),
#'   J3 = c(2,3,3,2,3,3), J4 = c(3,3,3,3,3,2)
#' )
#' res <- aiken_bayes_table(d, seed = 1, verbose = FALSE)
#' res$global
#'
#' @export
aiken_bayes_table <- function(data,
                              judge_cols = NULL, id_cols = NULL,
                              coefs = c("V", "H", "C"),
                              local = TRUE,
                              item_col = NULL,
                              l = 0, s = 3,
                              prior = c("perks", "jeffreys", "uniform"),
                              alpha0 = NULL,
                              cred_level = 0.95,
                              thresholds = c(0.70, 0.80),
                              ci_type = c("HDI", "ETI"),
                              B = 10000, seed = NULL, verbose = TRUE) {

  prior   <- match.arg(prior)
  ci_type <- match.arg(ci_type)
  coefs   <- match.arg(coefs, c("V", "H", "C"), several.ok = TRUE)
  data    <- as.data.frame(data)
  cols    <- detect_judge_cols(data, judge_cols, id_cols)
  J       <- cols$judges
  ID      <- cols$id
  sc      <- .ab_scale(l, s)
  if (!is.null(seed)) set.seed(seed)

  if (length(J) < 2)
    stop("Se necesitan al menos 2 columnas de jueces. Detectadas: ",
         paste(J, collapse = ", "), call. = FALSE)

  M <- as.matrix(data[, J, drop = FALSE])
  storage.mode(M) <- "double"

  if (verbose) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("  ANALISIS BAYESIANO DE COEFICIENTES DE AIKEN\n")
    cat(strrep("=", 70), "\n\n", sep = "")
    cat(sprintf("  Filas evaluadas   : %d\n", nrow(data)))
    cat(sprintf("  Jueces detectados : %d  (%s)\n", length(J), paste(J, collapse = ", ")))
    cat(sprintf("  Identificadores   : %s\n", if (length(ID)) paste(ID, collapse = ", ") else "(ninguno)"))
    cat(sprintf("  Escala            : %g a %g  (%d categorias)\n", l, s, sc$c))
    cat(sprintf("  Prior             : %s\n", prior))
    cat(sprintf("  Coeficientes      : %s%s\n", paste(coefs, collapse = ", "),
                if (local) " (+ A e I locales)" else ""))
    cat("\n")
  }

  # ---- fuera de rango -------------------------------------------------------
  fuera <- which(M < l | M > s, arr.ind = TRUE)
  if (nrow(fuera)) {
    stop("Hay ", nrow(fuera), " calificaciones fuera de la escala [", l, ", ", s,
         "]. Primera en la fila ", fuera[1, 1], ", juez '", J[fuera[1, 2]],
         "' (valor ", M[fuera[1, 1], fuera[1, 2]], "). Revise 'l' y 's'.",
         call. = FALSE)
  }

  # ---- GLOBAL: un coeficiente por fila --------------------------------------
  glob <- data[, ID, drop = FALSE]
  if (!length(ID)) glob <- data.frame(Fila = seq_len(nrow(data)))

  # Nombres cortos y ESTABLES, pensados para usar desde codigo:
  #   <coef>_clasico  valor muestral de Aiken
  #   <coef>          media posterior
  #   <coef>_li / _ls limites del intervalo
  #   <coef>_p70 ...  probabilidad de superar cada umbral (el sufijo es el umbral)
  # La presentacion bonita se genera aparte con aiken_apa(); asi el codigo no
  # depende del formato de la tabla.
  suf_p <- paste0("p", sub("0\\.", "", format(thresholds)))
  for (cf in coefs) {
    est <- matrix(NA_real_, nrow(data), 4 + length(thresholds))
    for (i in seq_len(nrow(data))) {
      r <- M[i, ]
      r <- r[!is.na(r)]
      if (length(r) < 2) next
      z <- aiken_bayes(r, coef = cf, l = l, s = s, prior = prior, alpha0 = alpha0,
                       cred_level = cred_level, thresholds = thresholds,
                       ci_type = ci_type, B = B, verbose = FALSE)
      est[i, ] <- c(z$clasico, z$media, z$ci[1], z$ci[2], z$probs)
    }
    colnames(est) <- c(paste0(cf, "_clasico"), cf,
                       paste0(cf, "_li"), paste0(cf, "_ls"),
                       paste0(cf, "_", suf_p))
    glob <- cbind(glob, round(as.data.frame(est), 3))
  }

  # Intervalo frecuentista de Penfield-Giacobbi, solo para la V y solo como
  # referencia comparativa. Es barato y permite ensenar la diferencia entre
  # credibilidad y confianza sin recalcular nada.
  if ("V" %in% coefs) {
    nj <- rowSums(!is.na(M))
    pg <- t(mapply(function(v, n) if (is.na(v) || n < 2) c(NA, NA)
                                  else aiken_V_ic_pg(v, n, sc$k, cred_level),
                   glob$V_clasico, nj))
    glob$V_pg_li <- round(pg[, 1], 3)
    glob$V_pg_ls <- round(pg[, 2], 3)
  }

  # ---- LOCAL: A (juez vs. los demas) e I (item vs. los demas) ---------------
  acuerdo <- coherencia <- NULL
  if (local) {
    if (verbose) cat("  Calculando coeficientes locales A e I...\n\n")
    filas_id <- if (length(ID)) apply(data[, ID, drop = FALSE], 1,
                                      function(z) paste(z, collapse = " | "))
                else as.character(seq_len(nrow(data)))

    # --- A: para cada fila y cada juez, acuerdo con los OTROS jueces ---------
    ac <- list()
    for (i in seq_len(nrow(data))) {
      for (jj in seq_along(J)) {
        foc <- M[i, jj]; oth <- M[i, -jj]
        oth <- oth[!is.na(oth)]
        if (is.na(foc) || length(oth) < 1) next
        z <- aiken_bayes(oth, focal = foc, coef = "A", l = l, s = s, prior = prior,
                         alpha0 = alpha0, cred_level = cred_level,
                         thresholds = thresholds, ci_type = ci_type, B = B,
                         verbose = FALSE)
        ac[[length(ac) + 1]] <- data.frame(
          Unidad = filas_id[i], Juez = J[jj],
          A_clasico = round(z$clasico, 3), A = round(z$media, 3),
          A_li = round(z$ci[1], 3), A_ls = round(z$ci[2], 3),
          stringsAsFactors = FALSE)
      }
    }
    acuerdo <- if (length(ac)) do.call(rbind, ac) else NULL

    # --- I: para cada juez, coherencia de cada item con los OTROS items ------
    ic_col <- if (!is.null(item_col)) item_col else if (length(ID)) ID[1] else NULL
    if (!is.null(ic_col) && ic_col %in% names(data)) {
      items <- unique(data[[ic_col]])
      if (length(items) >= 2) {
        co <- list()
        for (jj in seq_along(J)) {
          # una calificacion por item (media redondeada si hay varios criterios)
          por_item <- vapply(items, function(it) {
            v <- M[data[[ic_col]] == it, jj]
            v <- v[!is.na(v)]
            if (!length(v)) NA_real_ else round(mean(v))
          }, numeric(1))
          for (a in seq_along(items)) {
            foc <- por_item[a]; oth <- por_item[-a]
            oth <- oth[!is.na(oth)]
            if (is.na(foc) || length(oth) < 1) next
            z <- aiken_bayes(oth, focal = foc, coef = "I", l = l, s = s,
                             prior = prior, alpha0 = alpha0, cred_level = cred_level,
                             thresholds = thresholds, ci_type = ci_type, B = B,
                             verbose = FALSE)
            co[[length(co) + 1]] <- data.frame(
              Item = as.character(items[a]), Juez = J[jj],
              I_clasico = round(z$clasico, 3), I = round(z$media, 3),
              I_li = round(z$ci[1], 3), I_ls = round(z$ci[2], 3),
              stringsAsFactors = FALSE)
          }
        }
        coherencia <- if (length(co)) do.call(rbind, co) else NULL
      }
    }
  }

  out <- list(global = glob, acuerdo = acuerdo, coherencia = coherencia,
              jueces = J, id_cols = ID, escala = c(l = l, s = s, c = sc$c),
              prior = prior, coefs = coefs, thresholds = thresholds,
              cred_level = cred_level, ci_type = ci_type)
  class(out) <- "aiken_table"

  if (verbose) print(out)
  invisible(out)
}


#' Imprimir un objeto aiken_table
#' @param x Objeto \code{aiken_table}.
#' @param ... Ignorado.
#' @export
print.aiken_table <- function(x, ...) {
  cat(strrep("-", 70), "\n", sep = "")
  cat("  MEDIDAS GLOBALES (una fila por unidad evaluada)\n")
  cat(strrep("-", 70), "\n", sep = "")
  print_apa(x)

  if ("V" %in% x$coefs) {
    u <- x$thresholds[1]
    pv <- x$global[[paste0("V_p", sub("0\\.", "", format(u)))]]
    if (!is.null(pv)) {
      cat("\n  Sintesis de la validez de contenido\n")
      n_ok  <- sum(pv >= 0.90, na.rm = TRUE)
      n_med <- sum(pv >= 0.50 & pv < 0.90, na.rm = TRUE)
      n_no  <- sum(pv < 0.50, na.rm = TRUE)
      cat(sprintf("   %2d unidades con evidencia fuerte de V > %.2f\n", n_ok, u))
      cat(sprintf("   %2d unidades dudosas (revisar redaccion o ampliar el panel)\n", n_med))
      cat(sprintf("   %2d unidades sin respaldo (candidatas a eliminar)\n", n_no))
      if (n_no > 0) {
        idx <- which(pv < 0.50)
        et <- if (length(x$id_cols))
          apply(x$global[idx, x$id_cols, drop = FALSE], 1, paste, collapse = " | ")
          else as.character(idx)
        cat("   Sin respaldo: ", paste(et, collapse = "; "), "\n", sep = "")
      }
    }
  }

  if (!is.null(x$acuerdo)) {
    cat("\n", strrep("-", 70), "\n", sep = "")
    cat("  ACUERDO POR JUEZ (coeficiente A): promedio sobre las unidades\n")
    cat(strrep("-", 70), "\n", sep = "")
    ag <- aggregate(cbind(A_clasico, A) ~ Juez, data = x$acuerdo, FUN = mean)
    ag[, -1] <- round(ag[, -1], 3)
    ag <- ag[order(ag$A), ]
    print(ag, row.names = FALSE)
    peor <- ag$Juez[1]
    cat(sprintf("\n   El juez con menor acuerdo con el resto es '%s' (A = %.3f).\n",
                peor, ag$A[1]))
    cat("   Un A bajo de forma consistente sugiere revisar su criterio o su formacion.\n")
  }

  if (!is.null(x$coherencia)) {
    cat("\n", strrep("-", 70), "\n", sep = "")
    cat("  COHERENCIA POR ITEM (coeficiente I): promedio sobre los jueces\n")
    cat(strrep("-", 70), "\n", sep = "")
    ag <- aggregate(cbind(I_clasico, I) ~ Item, data = x$coherencia, FUN = mean)
    ag[, -1] <- round(ag[, -1], 3)
    print(ag[order(ag$I), ], row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}


#' Convertir un aiken_table a data.frame
#' @param x Objeto \code{aiken_table}.
#' @param row.names,optional,... Ignorados.
#' @export
as.data.frame.aiken_table <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$global
}

