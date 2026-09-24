# =============================================================================
# AikenBayes -- aplicacion Shiny para los seis coeficientes de Aiken bayesianos
# Se lanza con:  BayesAiken::run_aiken_app()
# =============================================================================

library(shiny)
library(BayesAiken)

PAL <- list(azul = "#00707F", rojo = "#CC2A36", gris = "#586E75",
            fondo = "#F7F9FA", verde = "#2E8B57", ambar = "#D68910")

CSS <- "
body { background: #F7F9FA; font-family: 'Segoe UI', Roboto, sans-serif; }
.titulo { background: #00707F; color: #fff; padding: 16px 22px; border-radius: 8px;
          margin-bottom: 18px; }
.titulo h2 { margin: 0; font-weight: 700; letter-spacing: .5px; }
.titulo p  { margin: 4px 0 0; opacity: .9; font-size: 14px; }
.panelbox { background: #fff; border: 1px solid #E1E8EB; border-radius: 8px;
            padding: 16px; margin-bottom: 16px; }
.panelbox h4 { margin-top: 0; color: #00707F; font-weight: 700;
               border-bottom: 2px solid #E1E8EB; padding-bottom: 8px; }
.paso { display: inline-block; background: #00707F; color: #fff; width: 24px;
        height: 24px; border-radius: 50%; text-align: center; line-height: 24px;
        margin-right: 8px; font-weight: 700; font-size: 13px; }
.consola { background: #10202A; color: #C8E6E8; font-family: Consolas, monospace;
           font-size: 12.5px; padding: 14px; border-radius: 6px;
           white-space: pre; overflow-x: auto; max-height: 640px; }
.aviso { background: #FFF6E5; border-left: 4px solid #D68910; padding: 10px 14px;
         border-radius: 4px; margin: 10px 0; font-size: 13px; }
.ok    { background: #EAF7F0; border-left: 4px solid #2E8B57; padding: 10px 14px;
         border-radius: 4px; margin: 10px 0; font-size: 13px; }
.btn-primary { background: #00707F !important; border-color: #00707F !important; }
table.data { width: 100%; border-collapse: collapse; font-size: 13px; }
table.data th { background: #00707F; color: #fff; padding: 7px 9px; text-align: left; }
table.data td { padding: 6px 9px; border-bottom: 1px solid #EDF1F3; }
table.data tr:nth-child(even) td { background: #FAFCFC; }
.badge-ok  { background:#2E8B57; color:#fff; padding:2px 8px; border-radius:10px; font-size:11px;}
.badge-med { background:#D68910; color:#fff; padding:2px 8px; border-radius:10px; font-size:11px;}
.badge-no  { background:#CC2A36; color:#fff; padding:2px 8px; border-radius:10px; font-size:11px;}
table.apa td { padding: 10px 14px; vertical-align: top; border-bottom: 1px solid #E8EEF0; }
table.apa th { background:#00707F; color:#fff; padding:9px 14px; text-align:left; }
.ab-val  { font-size:15px; font-weight:600; color:#12232B; letter-spacing:.2px; }
.ab-p    { font-size:12.5px; color:#00707F; }
.ab-clas { font-size:11.5px; color:#9AA9AE; }
.ab-pg   { font-size:11.5px; color:#B07A2A; font-style:italic; }
.nota    { font-size:12px; color:#586E75; font-style:italic; margin-top:10px; }
"

# ---- datos de demostracion ---------------------------------------------------
demo_data <- function() {
  data.frame(
    Item = paste0("Item", rep(1:6, each = 3)),
    Criterio = rep(c("Relevancia", "Representatividad", "Claridad"), 6),
    Juez1 = c(3,3,3, 3,3,2, 2,3,3, 3,3,3, 2,2,3, 3,3,3),
    Juez2 = c(2,3,3, 3,2,3, 3,2,3, 3,3,2, 3,2,2, 3,3,3),
    Juez3 = c(3,3,3, 2,3,3, 3,3,2, 3,2,3, 2,3,3, 3,2,3),
    Juez4 = c(2,3,3, 3,3,3, 2,3,3, 2,3,3, 3,2,2, 3,3,3),
    Juez5 = c(2,3,3, 3,2,2, 3,2,3, 3,3,3, 1,1,2, 3,3,3),
    stringsAsFactors = FALSE
  )
}

tabla_html <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return(HTML("<p><em>Sin datos.</em></p>"))
  df <- as.data.frame(df)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(z) formatC(z, format = "f", digits = digits))
  th <- paste0("<th>", names(df), "</th>", collapse = "")
  td <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  HTML(paste0("<table class='data'><thead><tr>", th, "</tr></thead><tbody>",
              paste(td, collapse = ""), "</tbody></table>"))
}

# =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML(CSS))),
  div(class = "titulo",
      h2("AikenBayes"),
      p("Los seis coeficientes de Aiken (V, H, C, R, A, I) con el modelo Dirichlet-Multinomial")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "panelbox",
        h4(span(class = "paso", "1"), "Datos"),
        radioButtons("origen", NULL,
                     c("Usar datos de demostracion" = "demo",
                       "Subir mi archivo" = "archivo"), selected = "demo"),
        conditionalPanel("input.origen == 'archivo'",
          fileInput("f1", "Archivo (.xlsx, .xls o .csv)",
                    accept = c(".xlsx", ".xls", ".csv"),
                    buttonLabel = "Examinar", placeholder = "Ningun archivo"),
          helpText("Formato: columnas identificadoras (Item, Criterio...) seguidas de una columna por juez.")
        ),
        uiOutput("ui_jueces")
      ),
      div(class = "panelbox",
        h4(span(class = "paso", "2"), "Escala"),
        fluidRow(
          column(6, numericInput("l", "Minimo", 0, step = 1)),
          column(6, numericInput("s", "Maximo", 3, step = 1))
        ),
        uiOutput("info_escala")
      ),
      div(class = "panelbox",
        h4(span(class = "paso", "3"), "Modelo"),
        selectInput("prior", "Prior Dirichlet",
                    c("Perks 1/c (recomendado)" = "perks",
                      "Jeffreys 1/2" = "jeffreys",
                      "Uniforme 1" = "uniform"), selected = "perks"),
        sliderInput("cred", "Nivel de credibilidad", 0.80, 0.99, 0.95, 0.01),
        textInput("umbrales", "Umbrales de decision", "0.70, 0.80"),
        selectInput("ci", "Tipo de intervalo", c("HDI", "ETI")),
        numericInput("B", "Muestras Monte Carlo", 10000, 1000, 100000, 1000),
        numericInput("semilla", "Semilla", 2026, step = 1),
        checkboxInput("local", "Calcular A e I locales", TRUE),
        checkboxInput("pg", "Mostrar tambien el IC clasico (Penfield-Giacobbi)", FALSE)
      ),
      actionButton("calcular", "Calcular", class = "btn btn-primary btn-block",
                   icon = icon("play")),
      br(),
      downloadButton("descargar", "Descargar informe (.xlsx)", class = "btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "pills",

        tabPanel("Guia",
          div(class = "panelbox",
            h4("Que hace esta aplicacion"),
            p(HTML("Calcula los seis coeficientes que propuso Lewis Aiken para analizar
                   calificaciones de jueces, pero con inferencia bayesiana en lugar de
                   tablas de probabilidad exacta.")),
            tags$table(class = "data",
              tags$thead(tags$tr(tags$th("Coef."), tags$th("Que mide"), tags$th("Se calcula sobre"))),
              tags$tbody(
                tags$tr(tags$td(tags$b("V")), tags$td("Validez de contenido de un item"), tags$td("Las calificaciones de los jueces a ese item")),
                tags$tr(tags$td(tags$b("H")), tags$td("Homogeneidad: cuanto coinciden los jueces entre si"), tags$td("Todas las parejas de jueces")),
                tags$tr(tags$td(tags$b("C")), tags$td("Congruencia: lo mismo que H pero penalizando mas las discrepancias grandes"), tags$td("La varianza de las calificaciones")),
                tags$tr(tags$td(tags$b("R")), tags$td("Reproducibilidad test-retest"), tags$td("Dos ocasiones de medicion")),
                tags$tr(tags$td(tags$b("A")), tags$td("Acuerdo de un juez concreto con los demas"), tags$td("Un juez contra los otros n-1")),
                tags$tr(tags$td(tags$b("I")), tags$td("Coherencia de un item con el resto del test"), tags$td("Un item contra los otros m-1"))
              )
            ),
            br(),
            h4("Como leer los resultados"),
            p(HTML("Cada coeficiente trae tres cosas. El <b>valor clasico</b> es la formula original
                   de Aiken, util como descriptivo. La <b>media bayesiana</b> estima el parametro de la
                   poblacion de jueces, no el valor de este panel concreto. La <b>probabilidad
                   posterior</b> responde la pregunta que de verdad importa: dada la evidencia,
                   cual es la probabilidad de que el coeficiente supere el umbral.")),
            div(class = "aviso", HTML("<b>Por que la media bayesiana no coincide con el valor clasico.</b>
                Con paneles pequenos el coeficiente muestral de H y C esta sesgado al alza,
                porque Aiken puso (n&sup2; - j) en el denominador para forzar el rango [0,1]
                exacto. La posterior estima el parametro poblacional. La diferencia es esperable
                y la aplicacion la cuantifica.")),
            div(class = "aviso", HTML("<b>El prior importa.</b> Con 5 jueces y 4 categorias, un prior
                Dirichlet(1,...,1) aporta 4 pseudo-observaciones, casi tanto peso como el panel
                entero. Por eso el valor por defecto es el prior de Perks, cuya masa total
                equivale a una sola observacion. Revise la pestana de sensibilidad."))
          )
        ),

        tabPanel("Datos",
          div(class = "panelbox", h4("Vista previa"), uiOutput("preview")),
          div(class = "panelbox", h4("Diagnostico"), uiOutput("diagnostico"))
        ),

        tabPanel("Medidas globales",
          div(class = "panelbox", h4("V, H y C por unidad evaluada"), uiOutput("tabla_global")),
          div(class = "panelbox", h4("Decision por item"), uiOutput("decision"))
        ),

        tabPanel("Medidas locales",
          div(class = "panelbox", h4("Acuerdo por juez (A)"), uiOutput("tabla_a")),
          div(class = "panelbox", h4("Coherencia por item (I)"), uiOutput("tabla_i"))
        ),

        tabPanel("Informe explicado",
          div(class = "panelbox",
            h4("Analisis paso a paso"),
            p("Seleccione una unidad y un coeficiente para ver el razonamiento completo."),
            fluidRow(
              column(5, uiOutput("ui_fila")),
              column(3, selectInput("coef_det", "Coeficiente", c("V", "H", "C"))),
              column(4, br(), actionButton("ver", "Ver informe", class = "btn btn-primary"))
            ),
            div(class = "consola", verbatimTextOutput("verbose", placeholder = TRUE))
          ),
          div(class = "panelbox", h4("Distribucion posterior"), plotOutput("plot_post", height = "420px"))
        ),

        tabPanel("Test-retest (R)",
          div(class = "panelbox",
            h4("Coeficiente R"),
            p("R necesita dos ocasiones de medicion. Suba el segundo archivo con la misma
              estructura y los mismos jueces."),
            div(class = "aviso", HTML("R <b>no</b> es funcional de la distribucion marginal:
                dos tablas con marginales identicas pueden dar R = 1.00 y R = 0.33. Por eso se
                modela la tabla conjunta completa de c&sup2; celdas. Con 4 categorias son 16
                celdas, asi que el prior de Perks es aqui especialmente importante.")),
            fileInput("f2", "Archivo de la segunda ocasion", accept = c(".xlsx", ".xls", ".csv")),
            actionButton("calc_r", "Calcular R", class = "btn btn-primary"),
            br(), br(),
            uiOutput("tabla_r")
          )
        ),

        tabPanel("Sensibilidad al prior",
          div(class = "panelbox",
            h4("Cuanto cambia la conclusion segun el prior"),
            p("Si las tres filas coinciden, el resultado esta gobernado por los datos.
              Si difieren, el panel es demasiado pequeno para ser concluyente."),
            uiOutput("ui_fila_sens"),
            actionButton("calc_sens", "Comparar priors", class = "btn btn-primary"),
            br(), br(),
            uiOutput("tabla_sens")
          )
        )
      )
    )
  )
)

# =============================================================================
server <- function(input, output, session) {

  vals <- reactiveValues(res = NULL, res_r = NULL, sens = NULL)

  datos <- reactive({
    if (input$origen == "demo") return(demo_data())
    req(input$f1)
    ext <- tolower(tools::file_ext(input$f1$name))
    if (ext == "csv") {
      utils::read.csv(input$f1$datapath, stringsAsFactors = FALSE)
    } else {
      as.data.frame(readxl::read_excel(input$f1$datapath))
    }
  })

  datos2 <- reactive({
    req(input$f2)
    ext <- tolower(tools::file_ext(input$f2$name))
    if (ext == "csv") utils::read.csv(input$f2$datapath, stringsAsFactors = FALSE)
    else as.data.frame(readxl::read_excel(input$f2$datapath))
  })

  output$ui_jueces <- renderUI({
    d <- try(datos(), silent = TRUE)
    if (inherits(d, "try-error") || is.null(d)) return(NULL)
    det <- try(detect_judge_cols(d), silent = TRUE)
    sel <- if (inherits(det, "try-error")) NULL else det$judges
    tagList(
      hr(),
      checkboxGroupInput("jueces", "Columnas de jueces", choices = names(d), selected = sel),
      helpText("Se detectan solas. Corrija si hace falta.")
    )
  })

  output$info_escala <- renderUI({
    l <- input$l; s <- input$s
    if (is.null(l) || is.null(s) || s <= l)
      return(div(class = "aviso", "El maximo debe ser mayor que el minimo."))
    cc <- s - l + 1
    nj <- length(input$jueces)
    # La granularidad de V es un hecho aritmetico, no un resultado de simulacion:
    # con n jueces y c categorias, V solo puede tomar n(c-1)+1 valores distintos.
    txt <- sprintf("<b>%d categorias</b> (k = %d).", cc, s - l)
    if (nj >= 2)
      txt <- paste0(txt, sprintf(" Con %d jueces, V solo puede tomar <b>%d valores
                     distintos</b> entre 0 y 1 (saltos de %.3f).",
                     nj, nj * (s - l) + 1, 1 / (nj * (s - l))))
    if (cc == 2)
      return(div(class = "aviso", HTML(paste0(txt, " Con solo 2 categorias la rejilla
        es muy gruesa y H y C pierden casi toda su capacidad de discriminar."))))
    div(class = "ok", HTML(txt))
  })

  output$preview <- renderUI({
    d <- try(datos(), silent = TRUE)
    if (inherits(d, "try-error")) return(div(class = "aviso", "No se pudo leer el archivo."))
    tabla_html(utils::head(d, 15), digits = 0)
  })

  output$diagnostico <- renderUI({
    d <- try(datos(), silent = TRUE)
    if (inherits(d, "try-error")) return(NULL)
    J <- input$jueces
    if (is.null(J) || length(J) < 2)
      return(div(class = "aviso", "Seleccione al menos 2 columnas de jueces."))
    M <- suppressWarnings(as.matrix(sapply(d[, J, drop = FALSE], as.numeric)))
    msgs <- list()
    fuera <- sum(M < input$l | M > input$s, na.rm = TRUE)
    if (fuera > 0)
      msgs <- c(msgs, list(div(class = "aviso", sprintf(
        "Hay %d calificaciones fuera de la escala [%g, %g]. Corrija la escala o los datos.",
        fuera, input$l, input$s))))
    nas <- sum(is.na(M))
    if (nas > 0)
      msgs <- c(msgs, list(div(class = "aviso", sprintf("Hay %d valores perdidos. Se omiten fila a fila.", nas))))
    if (length(J) < 5)
      msgs <- c(msgs, list(div(class = "aviso", sprintf(
        "Solo %d jueces. Con paneles tan pequenos la posterior depende bastante del prior.", length(J)))))
    usadas <- sort(unique(as.vector(M)))
    msgs <- c(msgs, list(div(class = "ok", HTML(sprintf(
      "<b>%d filas x %d jueces.</b> Categorias observadas: %s.",
      nrow(d), length(J), paste(usadas[!is.na(usadas)], collapse = ", "))))))
    do.call(tagList, msgs)
  })

  umbrales <- reactive({
    u <- suppressWarnings(as.numeric(trimws(strsplit(input$umbrales, ",")[[1]])))
    u <- u[!is.na(u) & u > 0 & u < 1]
    if (!length(u)) c(0.70, 0.80) else u
  })

  observeEvent(input$calcular, {
    d <- datos(); J <- input$jueces
    if (is.null(J) || length(J) < 2) {
      showNotification("Seleccione al menos 2 columnas de jueces.", type = "error"); return()
    }
    withProgress(message = "Calculando posteriores...", value = 0.3, {
      r <- try(aiken_bayes_table(
        d, judge_cols = J, id_cols = setdiff(names(d), J),
        l = input$l, s = input$s, prior = input$prior,
        cred_level = input$cred, thresholds = umbrales(),
        ci_type = input$ci, B = input$B, seed = input$semilla,
        local = input$local, verbose = FALSE), silent = TRUE)
      incProgress(0.7)
    })
    if (inherits(r, "try-error")) {
      showNotification(paste("Error:", conditionMessage(attr(r, "condition"))),
                       type = "error", duration = 12)
      return()
    }
    vals$res <- r
    updateTabsetPanel(session, "tabs", "Medidas globales")
    showNotification("Listo.", type = "message")
  })

  output$tabla_global <- renderUI({
    req(vals$res)
    t <- aiken_apa(vals$res, html = TRUE, pg = isTRUE(input$pg))
    u <- attr(t, "umbral")
    th <- paste0("<th>", names(t), "</th>", collapse = "")
    td <- apply(t, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
    HTML(paste0(
      "<table class='data apa'><thead><tr>", th, "</tr></thead><tbody>",
      paste(td, collapse = ""), "</tbody></table>",
      "<p class='nota'>Nota. Media posterior [intervalo de <b>credibilidad</b> del ",
      round(100 * vals$res$cred_level), " %, ", vals$res$ci_type,
      "]; P = probabilidad posterior de superar ", fmt_apa(u),
      ". No es un intervalo de confianza: contiene el valor del parámetro con esa ",
      "probabilidad, dados los datos y el prior. El valor clásico es la fórmula ",
      "original de Aiken, descriptiva.",
      if (isTRUE(input$pg))
        paste0(" <b>IC PG</b> es el intervalo de puntuación de Penfield-Giacobbi (2004), ",
               "frecuentista: supone n×k ensayos binomiales independientes. Se muestra ",
               "solo para comparar; suele salir algo más estrecho.")
      else "",
      "</p>"))
  })

  output$decision <- renderUI({
    req(vals$res)
    g <- vals$res$global; u <- umbrales()[1]
    cn <- paste0("V_p", sub("0\\.", "", format(u)))
    if (!cn %in% names(g)) return(NULL)
    p <- g[[cn]]
    et <- if (length(vals$res$id_cols))
      apply(g[, vals$res$id_cols, drop = FALSE], 1, paste, collapse = " | ")
      else as.character(seq_along(p))
    badge <- ifelse(p >= 0.90, "<span class='badge-ok'>RETENER</span>",
             ifelse(p >= 0.50, "<span class='badge-med'>REVISAR</span>",
                    "<span class='badge-no'>ELIMINAR</span>"))
    ic <- paste0("[", fmt_apa(g$V_li), ", ", fmt_apa(g$V_ls), "]")
    filas <- paste0("<tr><td>", et, "</td><td>", fmt_apa(g$V),
                    "</td><td>", ic, "</td><td>", fmt_apa(p),
                    "</td><td>", badge, "</td></tr>", collapse = "")
    HTML(paste0("<table class='data'><thead><tr><th>Unidad</th><th>V</th>",
                "<th>ICr ", round(100 * vals$res$cred_level), " %</th><th>P(V &gt; ",
                fmt_apa(u), ")</th><th>Decisión</th></tr></thead><tbody>",
                filas, "</tbody></table>",
                "<p class='nota'>",
                "RETENER: P &ge; .90 &nbsp;|&nbsp; REVISAR: .50 &le; P &lt; .90 &nbsp;|&nbsp; ",
                "ELIMINAR: P &lt; .50</p>"))
  })

  # Promedia sobre las unidades y devuelve la tabla en notacion APA.
  # OJO: los nombres son <coef>, <coef>_clasico, _li, _ls (antes eran _bayes,
  # _inf, _sup). Este bloque se quedo con los viejos al renombrar y reventaba
  # las dos pestanas locales.
  resumen_local <- function(d, por, coef) {
    f <- stats::as.formula(paste0("cbind(", coef, "_clasico, ", coef, ", ",
                                  coef, "_li, ", coef, "_ls) ~ ", por))
    ag <- aggregate(f, data = d, FUN = mean)
    ag <- ag[order(ag[[coef]]), ]
    data.frame(
      setNames(list(ag[[por]]), por),
      Coeficiente = paste0(fmt_apa(ag[[coef]]), " [", fmt_apa(ag[[paste0(coef, "_li")]]),
                           ", ", fmt_apa(ag[[paste0(coef, "_ls")]]), "]"),
      Clasico = fmt_apa(ag[[paste0(coef, "_clasico")]]),
      check.names = FALSE, stringsAsFactors = FALSE)
  }

  output$tabla_a <- renderUI({
    req(vals$res)
    if (is.null(vals$res$acuerdo))
      return(HTML("<p><em>Active el cálculo de medidas locales en el panel 3.</em></p>"))
    t <- resumen_local(vals$res$acuerdo, "Juez", "A")
    names(t) <- c("Juez", "A [ICr 95 %]", "Clásico")
    tagList(tabla_html(t),
            HTML(paste0("<p class='nota'>Promedio sobre las unidades evaluadas, de menor a ",
                        "mayor acuerdo. El juez de arriba es el que más se desvía del resto.</p>")))
  })

  output$tabla_i <- renderUI({
    req(vals$res)
    if (is.null(vals$res$coherencia))
      return(HTML("<p><em>Se necesitan al menos 2 ítems para calcular la I.</em></p>"))
    t <- resumen_local(vals$res$coherencia, "Item", "I")
    names(t) <- c("Ítem", "I [ICr 95 %]", "Clásico")
    tagList(tabla_html(t),
            HTML(paste0("<p class='nota'>Promedio sobre los jueces, de menor a mayor ",
                        "coherencia. El ítem de arriba es el que peor encaja con el resto.</p>")))
  })

  etiquetas_filas <- reactive({
    d <- datos(); J <- input$jueces
    id <- setdiff(names(d), J)
    if (length(id)) apply(d[, id, drop = FALSE], 1, paste, collapse = " | ")
    else paste("Fila", seq_len(nrow(d)))
  })

  output$ui_fila <- renderUI({
    e <- try(etiquetas_filas(), silent = TRUE)
    if (inherits(e, "try-error")) return(NULL)
    selectInput("fila", "Unidad evaluada", stats::setNames(seq_along(e), e))
  })

  output$ui_fila_sens <- renderUI({
    e <- try(etiquetas_filas(), silent = TRUE)
    if (inherits(e, "try-error")) return(NULL)
    fluidRow(
      column(6, selectInput("fila_s", "Unidad evaluada", stats::setNames(seq_along(e), e))),
      column(6, selectInput("coef_s", "Coeficiente", c("V", "H", "C")))
    )
  })

  detalle <- eventReactive(input$ver, {
    d <- datos(); J <- input$jueces
    req(J, input$fila)
    r <- suppressWarnings(as.numeric(d[as.integer(input$fila), J]))
    aiken_bayes(r, coef = input$coef_det, l = input$l, s = input$s,
                prior = input$prior, cred_level = input$cred,
                thresholds = umbrales(), ci_type = input$ci,
                B = input$B, seed = input$semilla, verbose = FALSE)
  })

  output$verbose <- renderText({
    z <- try(detalle(), silent = TRUE)
    if (inherits(z, "try-error"))
      return("Pulse 'Ver informe' despues de elegir una unidad.")
    paste(utils::capture.output(print(z)), collapse = "\n")
  })

  output$plot_post <- renderPlot({
    z <- try(detalle(), silent = TRUE)
    if (inherits(z, "try-error")) return(NULL)
    plot(z)
  })

  observeEvent(input$calc_r, {
    d1 <- datos(); d2 <- try(datos2(), silent = TRUE)
    if (inherits(d2, "try-error")) { showNotification("Suba el segundo archivo.", type = "error"); return() }
    J <- input$jueces
    if (!all(J %in% names(d2))) { showNotification("El segundo archivo no tiene los mismos jueces.", type = "error"); return() }
    if (nrow(d1) != nrow(d2)) { showNotification("Los dos archivos deben tener el mismo numero de filas.", type = "error"); return() }
    id <- setdiff(names(d1), J)
    et <- if (length(id)) apply(d1[, id, drop = FALSE], 1, paste, collapse = " | ")
          else paste("Fila", seq_len(nrow(d1)))
    withProgress(message = "Calculando R...", value = 0.3, {
      out <- do.call(rbind, lapply(seq_len(nrow(d1)), function(i) {
        x <- suppressWarnings(as.numeric(d1[i, J]))
        y <- suppressWarnings(as.numeric(d2[i, J]))
        z <- try(aiken_bayes(x, y = y, coef = "R", l = input$l, s = input$s,
                             prior = input$prior, cred_level = input$cred,
                             thresholds = umbrales(), ci_type = input$ci,
                             B = input$B, seed = input$semilla, verbose = FALSE),
                 silent = TRUE)
        if (inherits(z, "try-error")) return(NULL)
        data.frame(Unidad = et[i], R_clasico = round(z$clasico, 3),
                   R = round(z$media, 3), R_li = round(z$ci[1], 3),
                   R_ls = round(z$ci[2], 3), stringsAsFactors = FALSE)
      }))
      incProgress(0.7)
    })
    vals$res_r <- out
  })

  output$tabla_r <- renderUI({
    if (is.null(vals$res_r)) return(HTML("<p><em>Sin resultados todavia.</em></p>"))
    tabla_html(vals$res_r)
  })

  observeEvent(input$calc_sens, {
    d <- datos(); J <- input$jueces
    req(J, input$fila_s)
    r <- suppressWarnings(as.numeric(d[as.integer(input$fila_s), J]))
    out <- do.call(rbind, lapply(c("perks", "jeffreys", "uniform"), function(pp) {
      z <- aiken_bayes(r, coef = input$coef_s, l = input$l, s = input$s, prior = pp,
                       cred_level = input$cred, thresholds = umbrales(),
                       ci_type = input$ci, B = input$B, seed = input$semilla,
                       verbose = FALSE)
      data.frame(Prior = z$prior$nombre,
                 Masa_a_priori = round(z$prior$masa, 2),
                 Media = round(z$media, 3),
                 Inferior = round(z$ci[1], 3), Superior = round(z$ci[2], 3),
                 P_umbral = round(z$probs[1], 3), stringsAsFactors = FALSE)
    }))
    vals$sens <- out
  })

  output$tabla_sens <- renderUI({
    if (is.null(vals$sens)) return(HTML("<p><em>Sin resultados todavia.</em></p>"))
    rango <- diff(range(vals$sens$P_umbral))
    nota <- if (rango > 0.15)
      div(class = "aviso", sprintf("La probabilidad varia %.2f entre priors. El panel es
          demasiado pequeno para una conclusion firme: amplie el numero de jueces.", rango))
    else
      div(class = "ok", sprintf("La probabilidad solo varia %.2f entre priors. El resultado
          esta gobernado por los datos.", rango))
    tagList(tabla_html(vals$sens), nota)
  })

  output$descargar <- downloadHandler(
    filename = function() paste0("AikenBayes_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(vals$res)
      # Dos hojas de resultados: la legible para el anexo del artículo y los
      # números crudos para recalcular.
      apa <- aiken_apa(vals$res, sep = "\n", pg = isTRUE(input$pg))
      hojas <- list(`Tabla APA` = apa, `Medidas globales` = vals$res$global)
      if (!is.null(vals$res$acuerdo))    hojas$`Acuerdo (A)`     <- vals$res$acuerdo
      if (!is.null(vals$res$coherencia)) hojas$`Coherencia (I)`  <- vals$res$coherencia
      if (!is.null(vals$res_r))          hojas$`Test-retest (R)` <- vals$res_r
      if (!is.null(vals$sens))           hojas$`Sensibilidad`    <- vals$sens
      hojas$`Configuracion` <- data.frame(
        Parametro = c("Escala minima", "Escala maxima", "Categorias", "Prior",
                      "Nivel de credibilidad", "Tipo de intervalo (credibilidad)",
                      "Umbrales", "Muestras Monte Carlo", "Semilla",
                      "Coeficientes A e I locales", "IC clasico Penfield-Giacobbi",
                      "Paquete", "Generado"),
        Valor = c(input$l, input$s, input$s - input$l + 1, input$prior,
                  input$cred, input$ci, input$umbrales, input$B, input$semilla,
                  if (isTRUE(input$local)) "si" else "no",
                  if (isTRUE(input$pg)) "si" else "no",
                  paste("BayesAiken", utils::packageVersion("BayesAiken")),
                  format(Sys.time(), "%Y-%m-%d %H:%M")),
        stringsAsFactors = FALSE)

      # writexl no aplica formato, de modo que los saltos de linea de la hoja
      # APA se ven pegados en una sola linea. Con openxlsx se activa el ajuste
      # de texto y la hoja queda legible; si no esta instalado, se sustituyen
      # los saltos por un separador para que al menos no salga amontonado.
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        wb <- openxlsx::createWorkbook()
        estilo_cab <- openxlsx::createStyle(fgFill = "#00707F", fontColour = "#FFFFFF",
                                            textDecoration = "bold", halign = "left",
                                            valign = "center", border = "bottom")
        # Sin borde inferior las filas de varias lineas se confunden entre si y
        # la ultima linea parece pertenecer al item siguiente.
        estilo_apa <- openxlsx::createStyle(wrapText = TRUE, valign = "top",
                                            border = "bottom",
                                            borderColour = "#D6E0E3", borderStyle = "thin")
        for (nm in names(hojas)) {
          d <- hojas[[nm]]
          openxlsx::addWorksheet(wb, nm)
          openxlsx::writeData(wb, nm, d, headerStyle = estilo_cab)
          openxlsx::freezePane(wb, nm, firstRow = TRUE)
          if (identical(nm, "Tabla APA")) {
            openxlsx::addStyle(wb, nm, estilo_apa, rows = 2:(nrow(d) + 1),
                               cols = seq_len(ncol(d)), gridExpand = TRUE)
            openxlsx::setColWidths(wb, nm, seq_len(ncol(d)),
                                   widths = c(rep(18, ncol(d) - 3), rep(26, 3)))
            # ~16 puntos por linea + margen; con IC de PG son 4 lineas, no 3.
            n_lineas <- if (isTRUE(input$pg)) 4 else 3
            openxlsx::setRowHeights(wb, nm, 2:(nrow(d) + 1),
                                    heights = n_lineas * 16 + 8)
          } else {
            openxlsx::setColWidths(wb, nm, seq_len(ncol(d)), widths = "auto")
          }
        }
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      } else {
        hojas$`Tabla APA` <- aiken_apa(vals$res, sep = " · ", pg = isTRUE(input$pg))
        writexl::write_xlsx(hojas, path = file)
      }
    }
  )
}

shinyApp(ui, server)
