# =============================================================================
# BayesAiken - Ejemplo Reproducible Completo
# =============================================================================
# Este script demuestra todas las funcionalidades de la libreria BayesAiken
# para calcular intervalos de credibilidad bayesianos para coeficientes de
# validez de contenido.
#
# Autor: Jose Ventura-Leon
# Fecha: 2024-12-30
# =============================================================================

# =============================================================================
# PASO 0: INSTALACION DE LA LIBRERIA
# =============================================================================
# Primero instalar devtools si no lo tienes
# install.packages("devtools")

# Instalar BayesAiken desde la ruta local
devtools::install("D:/14. LIBRERIAS/BayesAiken")

# O si esta en GitHub:
# devtools::install_github("jventural/BayesAiken")

# =============================================================================
# PASO 1: CARGAR LA LIBRERIA
# =============================================================================
library(BayesAiken)

# =============================================================================
# SECCION 1: DATOS DE EJEMPLO
# =============================================================================
# Simulamos un estudio de validez de contenido donde 10 jueces expertos
# evaluan 5 items de un instrumento psicologico en una escala de 0-3
# (0=nada relevante, 1=poco relevante, 2=relevante, 3=muy relevante)

set.seed(2024)  # Para reproducibilidad

# Crear datos simulados realistas
datos_validez <- data.frame(
  item = rep(paste0("Item_", 1:5), each = 3),
  criterio = rep(c("Relevancia", "Representatividad", "Claridad"), 5),
  J1 = c(3,3,3, 3,3,2, 3,2,3, 2,3,3, 3,3,2),
  J2 = c(3,3,2, 3,2,3, 3,3,2, 3,3,3, 2,3,3),
  J3 = c(3,2,3, 2,3,3, 3,3,3, 3,2,3, 3,2,3),
  J4 = c(2,3,3, 3,3,3, 2,3,3, 3,3,2, 3,3,3),
  J5 = c(3,3,3, 3,3,2, 3,2,3, 3,3,3, 3,3,3),
  J6 = c(3,3,2, 3,2,3, 3,3,2, 2,3,3, 3,2,3),
  J7 = c(3,2,3, 2,3,3, 3,3,3, 3,3,3, 2,3,3),
  J8 = c(2,3,3, 3,3,3, 2,3,2, 3,2,3, 3,3,2),
  J9 = c(3,3,3, 3,3,2, 3,2,3, 3,3,3, 3,3,3),
  J10 = c(3,3,2, 3,2,3, 3,3,3, 3,3,2, 3,2,3)
)

cat("\n")
cat("================================================================\n")
cat("           EJEMPLO REPRODUCIBLE - BayesAiken v1.0.0\n")
cat("================================================================\n")
cat("\nDatos de ejemplo:\n")
print(head(datos_validez, 10))

# =============================================================================
# SECCION 2: COEFICIENTE V DE AIKEN (UN SOLO ITEM)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 2: COEFICIENTE V DE AIKEN (ITEM INDIVIDUAL)\n")
cat("================================================================\n")

# Extraer calificaciones del Item 1, criterio Relevancia
ratings_item1 <- as.numeric(datos_validez[datos_validez$item == "Item_1" &
                                           datos_validez$criterio == "Relevancia",
                                           paste0("J", 1:10)])

cat("\nCalificaciones Item 1 (Relevancia):", ratings_item1, "\n")

# Calcular V de Aiken con enfoque bayesiano
resultado_V <- coef_V(ratings_item1, l = 0, s = 3,
                      prior_alpha = 1, prior_beta = 1,
                      cred_level = 0.95, verbose = TRUE)

# =============================================================================
# SECCION 3: TODOS LOS COEFICIENTES PARA UN ITEM
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 3: MULTIPLES COEFICIENTES (V, H, C)\n")
cat("================================================================\n")

# Usar la funcion principal para calcular V, H y C simultaneamente
resultado_multiple <- bayes_aiken(ratings_item1,
                                   coefficient = c("V", "H", "C"),
                                   l = 0, s = 3, c = 4,
                                   verbose = TRUE)

# =============================================================================
# SECCION 4: COEFICIENTE R (REPRODUCIBILIDAD TEST-RETEST)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 4: COEFICIENTE R (REPRODUCIBILIDAD)\n")
cat("================================================================\n")

# Simular calificaciones en dos momentos temporales
tiempo_1 <- c(3, 3, 2, 3, 3, 2, 3, 3, 3, 3)
tiempo_2 <- c(3, 3, 3, 3, 2, 2, 3, 3, 3, 3)

cat("\nCalificaciones Tiempo 1:", tiempo_1)
cat("\nCalificaciones Tiempo 2:", tiempo_2, "\n")

resultado_R <- coef_R(tiempo_1, tiempo_2, c = 4, verbose = TRUE)

# =============================================================================
# SECCION 5: COEFICIENTE A (ACUERDO DEL JUEZ)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 5: COEFICIENTE A (ACUERDO DEL JUEZ)\n")
cat("================================================================\n")

# Evaluar el acuerdo del Juez 1 con los demas jueces
juez_focal <- c(3, 3, 2)  # Calificaciones del Juez 1 en 3 criterios

# Calificaciones de otros jueces (matriz: filas=criterios, columnas=jueces)
otros_jueces <- matrix(c(
  3, 3, 2, 3, 3, 2, 3, 2, 3,  # Criterio 1: J2-J10
  3, 2, 3, 3, 3, 2, 3, 3, 3,  # Criterio 2: J2-J10
  2, 3, 3, 3, 2, 3, 3, 3, 2   # Criterio 3: J2-J10
), nrow = 3, byrow = TRUE)

cat("\nCalificaciones Juez focal:", juez_focal)
cat("\nMatriz otros jueces:\n")
print(otros_jueces)

resultado_A <- coef_A(juez_focal, otros_jueces, c = 4, verbose = TRUE)
# extract_results(resultado_V, resultado_R, )
# =============================================================================
# SECCION 6: COEFICIENTE I (CONSISTENCIA INTER-ITEM)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 6: COEFICIENTE I (CONSISTENCIA INTER-ITEM)\n")
cat("================================================================\n")

# Comparar consistencia entre Item 1 e Item 2
item1_ratings <- c(3, 3, 3, 2, 3, 3, 3, 2, 3, 3)
item2_ratings <- c(3, 3, 2, 3, 3, 3, 2, 3, 3, 3)

cat("\nCalificaciones Item 1:", item1_ratings)
cat("\nCalificaciones Item 2:", item2_ratings, "\n")

resultado_I <- coef_I(item1_ratings, item2_ratings, c = 4, verbose = TRUE)

# =============================================================================
# SECCION 7: ANALISIS CON DATA.FRAME (MULTIPLES ITEMS)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 7: ANALISIS DE MULTIPLES ITEMS (DATA.FRAME)\n")
cat("================================================================\n")

# Analizar todos los items del data.frame
resultado_df <- bayes_aiken(datos_validez,
                             coefficient = "V",
                             item_col = "item",
                             criterion_col = "criterio",
                             l = 0, s = 3,
                             verbose = TRUE)

# =============================================================================
# SECCION 8: COMPARACION DE PRIORS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 8: SENSIBILIDAD A DIFERENTES PRIORS\n")
cat("================================================================\n")

cat("\nComparando el efecto de diferentes distribuciones prior:\n")

comparacion_priors <- compare_priors(ratings_item1, l = 0, s = 3,
                                      priors = list(
                                        "Uniforme Beta(1,1)" = c(1, 1),
                                        "Jeffreys Beta(0.5,0.5)" = c(0.5, 0.5),
                                        "Informativo Beta(2,2)" = c(2, 2),
                                        "Esceptico Beta(1,3)" = c(1, 3),
                                        "Optimista Beta(3,1)" = c(3, 1)
                                      ))
print(comparacion_priors)

# =============================================================================
# SECCION 9: VISUALIZACION
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 9: VISUALIZACION DE RESULTADOS\n")
cat("================================================================\n")

# Grafico de la distribucion posterior para V (individual)
cat("\n9.1 Grafico de distribucion posterior individual:\n")
plot_posterior(resultado_V$post_alpha, resultado_V$post_beta,
               coef = "V",
               prior_alpha = 1, prior_beta = 1,
               cred_level = 0.95,
               show_prior = TRUE,
               show_hdi = TRUE,
               main = "Distribucion Posterior - V de Aiken (Item 1)")

# Grafico de multiples items en un panel
cat("\n9.2 Panel de multiples items (20 items simulados):\n")

# Crear datos para 20 items
items_data <- lapply(1:20, function(i) {
  list(
    post_alpha = 25 + sample(-5:10, 1),
    post_beta = 3 + sample(0:3, 1),
    item_name = paste0("Item ", i)
  )
})

# Panel de 4 filas x 5 columnas
plot_multiple_posteriors(items_data,
                         n_cols = 5,
                         main_title = "Posterior Distributions - Items 1-20",
                         coef = "V")

# =============================================================================
# SECCION 10: EXPORTACION DE RESULTADOS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 10: EXPORTACION DE RESULTADOS\n")
cat("================================================================\n")

# Capturar salida verbose
cat("\nCapturando salida verbose...\n")
captura <- bayes_aiken_capture(ratings_item1, coefficient = "V", l = 0, s = 3)

# Mostrar texto capturado (primeras lineas)
cat("\nPrimeras 20 lineas del texto capturado:\n")
cat(paste(head(captura$text, 20), collapse = "\n"))
cat("\n...\n")

# Para exportar a archivo (descomentar para usar):
# bayes_aiken_export(resultado_multiple, "resultados_validez.txt", format = "both")

# =============================================================================
# SECCION 11: TABLA RESUMEN FINAL
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 11: TABLA RESUMEN FINAL\n")
cat("================================================================\n")

cat("\nResumen del analisis:\n")
tabla_resumen <- summary_table(resultado_multiple)
print(tabla_resumen)

# =============================================================================
# SECCION 12: INTERPRETACION DE RESULTADOS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 12: GUIA DE INTERPRETACION\n")
cat("================================================================\n")

cat("
INTERPRETACION DE RESULTADOS BAYESIANOS:

1. ESTIMACION PUNTUAL:
   - V Clasico: Estimacion tradicional de Aiken
   - V Bayesiano: Media de la distribucion posterior

2. INTERVALOS DE CREDIBILIDAD:
   - ETI (Equal-Tailed Interval): Percentiles simetricos
   - HDI (Highest Density Interval): Region de mayor densidad
   - El HDI es preferible porque contiene los valores mas probables

3. PROBABILIDADES:
   - P(V > 0.70): Probabilidad de validez adecuada
   - P(V > 0.80): Probabilidad de validez excelente

4. CRITERIOS DE DECISION:
   - Si HDI_inferior >= 0.70: Evidencia fuerte de validez adecuada
   - Si P(V > 0.70) >= 0.95: Alta confianza en validez adecuada
   - Si P(V > 0.80) >= 0.95: Alta confianza en validez excelente

5. VENTAJAS DEL ENFOQUE BAYESIANO:
   - Proporciona probabilidades directas (no p-valores)
   - Incorpora incertidumbre de forma natural
   - Permite incorporar conocimiento previo via priors
   - Los intervalos tienen interpretacion probabilistica directa
")

# =============================================================================
# SECCION 13: FUNCIONES AVANZADAS - ANALISIS BAYESIANO EXTENDIDO
# =============================================================================
# Estas funciones implementan analisis avanzados basados en las mejores
# practicas de la literatura metodologica bayesiana para validez de contenido.
#
# Referencias:
# - Penfield & Giacobbi (2004) - Intervalos de confianza para V de Aiken
# - Kruschke (2018) - ROPE y decisiones bayesianas
# - SciELO Costa Rica - Comparaciones de coeficientes V

cat("\n")
cat("================================================================\n")
cat("  SECCION 13: FUNCIONES AVANZADAS\n")
cat("================================================================\n")

# -----------------------------------------------------------------------------
# 13.1 REGLA DE DECISION PROBABILISTICA
# -----------------------------------------------------------------------------
# En lugar de usar intervalos, usamos probabilidades directas:
# "Conservar el item si P(V >= c | datos) >= 0.95"
#
# Esto es mas natural en el marco bayesiano porque:
# - Proporciona una probabilidad directa de cumplir el criterio
# - Permite reglas de decision explicitas y reproducibles
# - Es facil de comunicar a no-estadisticos

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  13.1 REGLA DE DECISION PROBABILISTICA\n")
cat("----------------------------------------------------------------\n")
cat("\nPregunta: ¿El item cumple con el criterio de validez?\n")
cat("Regla: Conservar item si P(V >= 0.70) >= 0.95\n\n")

# Aplicar regla de decision al Item 1
decision_item1 <- decision_rule(resultado_V,
                                 threshold = 0.70,
                                 min_prob = 0.95,
                                 verbose = TRUE)

# Tambien podemos usar un umbral mas estricto
cat("\n--- Con umbral mas estricto (0.80): ---\n")
decision_estricta <- decision_rule(resultado_V,
                                    threshold = 0.80,
                                    min_prob = 0.95,
                                    verbose = TRUE)

# Aplicar a multiples items (si tenemos objeto bayes_aiken con data.frame)
cat("\n--- Aplicando regla a multiples items: ---\n")
# Nota: Esto funcionaria con resultado_df si tiene la estructura correcta
# decisiones_multi <- decision_rule_multi(resultado_df, threshold = 0.70, min_prob = 0.95)
# print(decisiones_multi)

# -----------------------------------------------------------------------------
# 13.2 COMPARACION DE ITEMS (Diferencias Posteriores)
# -----------------------------------------------------------------------------
# En vez de comparar puntajes puntuales, estimamos la distribucion posterior
# de la diferencia: Delta = V_item1 - V_item2
#
# Esto nos permite responder:
# - P(Delta > 0): Probabilidad de que Item1 sea mejor que Item2
# - P(Delta > delta): Probabilidad de diferencia clinicamente relevante

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  13.2 COMPARACION DE ITEMS (Diferencias Posteriores)\n")
cat("----------------------------------------------------------------\n")
cat("\nPregunta: ¿El Item 1 tiene mejor validez que el Item 2?\n")
cat("Calculamos: P(V_item1 - V_item2 > 0)\n\n")

# Item 1: Alta validez (el que ya calculamos)
cat("Item 1 - Calificaciones:", ratings_item1, "\n")

# Item 2: Validez moderada (simulado)
ratings_item2_comp <- c(2, 3, 2, 3, 2, 2, 3, 2, 3, 2)
cat("Item 2 - Calificaciones:", ratings_item2_comp, "\n\n")

# Calcular V para Item 2
resultado_item2 <- coef_V(ratings_item2_comp, l = 0, s = 3, verbose = FALSE)

# Comparar items
comparacion_items <- compare_items(resultado_V, resultado_item2,
                                    delta = 0.10,  # Diferencia minima relevante
                                    n_samples = 10000,
                                    verbose = TRUE)

# Visualizar la distribucion de diferencias
cat("\n--- Grafico de la distribucion de diferencias: ---\n")
plot(comparacion_items, show_rope = TRUE, rope = c(-0.10, 0.10))

# -----------------------------------------------------------------------------
# 13.3 COMPARACION DE GRUPOS DE JUECES
# -----------------------------------------------------------------------------
# Si tenemos dos paneles independientes de jueces (ej. clinicos vs
# investigadores, Lima vs provincias, expertos vs experienciales),
# podemos comparar sus evaluaciones directamente.
#
# Esto responde preguntas como:
# - ¿Los clinicos son mas exigentes que los investigadores?
# - ¿Hay sesgo geografico en las evaluaciones?
# - ¿Expertos y usuarios perciben la validez de forma similar?

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  13.3 COMPARACION DE GRUPOS DE JUECES\n")
cat("----------------------------------------------------------------\n")
cat("\nEscenario: Comparar evaluaciones de clinicos vs investigadores\n")
cat("para el mismo item.\n\n")

# Grupo 1: Clinicos (5 jueces)
grupo_clinicos <- c(3, 3, 2, 3, 3)
cat("Clinicos (n=5):", grupo_clinicos, "- Media:", mean(grupo_clinicos), "\n")

# Grupo 2: Investigadores (5 jueces)
grupo_investigadores <- c(3, 2, 3, 2, 2)
cat("Investigadores (n=5):", grupo_investigadores, "- Media:", mean(grupo_investigadores), "\n\n")

# Comparar grupos
comp_grupos <- compare_judge_groups(grupo_clinicos, grupo_investigadores,
                                     l = 0, s = 3,
                                     delta = 0.05,  # Diferencia para equivalencia
                                     n_samples = 10000,
                                     verbose = TRUE)

# Ejemplo 2: Grupos con diferencia clara
cat("\n--- Ejemplo con diferencia mas marcada: ---\n")
grupo_A <- c(3, 3, 3, 3, 3, 3, 3, 2)  # Muy positivos
grupo_B <- c(2, 2, 2, 3, 2, 2, 2, 2)  # Mas criticos

comp_marcada <- compare_judge_groups(grupo_A, grupo_B,
                                      l = 0, s = 3,
                                      delta = 0.05,
                                      verbose = TRUE)

# -----------------------------------------------------------------------------
# 13.4 ANALISIS ROPE (Region de Equivalencia Practica)
# -----------------------------------------------------------------------------
# El enfoque ROPE (Region of Practical Equivalence) combina:
# - HDI (Highest Density Interval): Region con el 95% de la posterior
# - ROPE: Region donde los valores son "practicamente aceptables"
#
# Decisiones:
# - ACEPTAR: Todo el HDI esta dentro del ROPE (evidencia fuerte)
# - RECHAZAR: Todo el HDI esta fuera del ROPE (evidencia contra)
# - INDECISO: El HDI se superpone parcialmente (necesitamos mas datos)
#
# Esto es mas sofisticado que simplemente mirar si el intervalo cruza
# un umbral, porque considera toda la distribucion.

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  13.4 ANALISIS ROPE (Region de Equivalencia Practica)\n")
cat("----------------------------------------------------------------\n")
cat("\nPregunta: ¿El item tiene validez 'practicamente aceptable'?\n")
cat("ROPE definido como: [0.70, 1.00] (V >= 0.70 es aceptable)\n\n")

# Analisis ROPE para Item 1 (alta validez)
rope_item1 <- rope_analysis(resultado_V,
                             rope_lower = 0.70,
                             rope_upper = 1.00,
                             cred_level = 0.95,
                             verbose = TRUE)

# Visualizar ROPE
cat("\n--- Grafico del analisis ROPE: ---\n")
plot(rope_item1)

# Ejemplo con item de validez dudosa
cat("\n--- Ejemplo con item de validez dudosa: ---\n")
ratings_dudoso <- c(2, 2, 3, 2, 2, 3, 2, 2, 3, 2)
resultado_dudoso <- coef_V(ratings_dudoso, l = 0, s = 3, verbose = FALSE)

rope_dudoso <- rope_analysis(resultado_dudoso,
                              rope_lower = 0.70,
                              rope_upper = 1.00,
                              verbose = TRUE)

# Ejemplo con ROPE centrado (para comparacion con un valor teorico)
cat("\n--- ROPE centrado en 0.80 (validez excelente): ---\n")
rope_centrado <- rope_analysis(resultado_V,
                                rope_lower = 0.75,
                                rope_upper = 0.85,
                                verbose = TRUE)

# -----------------------------------------------------------------------------
# 13.5 PLANIFICACION DE MUESTRA (Posterior Predictivo)
# -----------------------------------------------------------------------------
# Una de las ventajas mas practicas del enfoque bayesiano es poder
# responder: "¿Cuantos jueces mas necesito?"
#
# Usando el posterior predictivo, simulamos:
# 1. Muestreamos un V "verdadero" de nuestra posterior actual
# 2. Simulamos calificaciones de m jueces adicionales
# 3. Actualizamos la posterior con los nuevos datos
# 4. Calculamos P(V >= umbral) con los datos actualizados
#
# Esto nos da una justificacion formal para el tamano del panel,
# en lugar de "usamos 8 jueces porque es lo tipico".

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  13.5 PLANIFICACION DE MUESTRA (Posterior Predictivo)\n")
cat("----------------------------------------------------------------\n")
cat("\nEscenario: Tenemos solo 5 jueces y queremos saber cuantos mas\n")
cat("necesitamos para tener 95% de confianza de que V >= 0.70\n\n")

# Situacion actual: pocos jueces con resultados mixtos
ratings_actuales <- c(3, 2, 3, 3, 2)
cat("Calificaciones actuales (n=5):", ratings_actuales, "\n")
cat("Media actual:", mean(ratings_actuales), "\n\n")

resultado_actual <- coef_V(ratings_actuales, l = 0, s = 3, verbose = FALSE)

# Planificar cuantos jueces adicionales necesitamos
plan_muestra <- sample_size_planning(resultado_actual,
                                      target_prob = 0.95,
                                      threshold = 0.70,
                                      additional_judges = 1:15,
                                      n_simulations = 1000,
                                      verbose = TRUE)

# Visualizar curva de planificacion
cat("\n--- Grafico de planificacion de muestra: ---\n")
plot(plan_muestra)

# Escenario alternativo: umbral mas estricto
cat("\n--- Planificacion con umbral mas estricto (0.80): ---\n")
plan_estricto <- sample_size_planning(resultado_actual,
                                       target_prob = 0.95,
                                       threshold = 0.80,
                                       additional_judges = 1:20,
                                       n_simulations = 1000,
                                       verbose = TRUE)

# -----------------------------------------------------------------------------
# 13.6 RESUMEN DE FUNCIONES AVANZADAS
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  RESUMEN DE FUNCIONES AVANZADAS DISPONIBLES\n")
cat("================================================================\n")

cat("
FUNCIONES AVANZADAS DE BayesAiken:

1. decision_rule(result, threshold, min_prob)
   - Regla de decision: 'Aceptar si P(V >= threshold) >= min_prob'
   - Retorna: ACEPTAR, RECHAZAR, o INDECISO

2. decision_rule_multi(results, threshold, min_prob)
   - Aplica regla de decision a multiples items
   - Retorna: data.frame con decisiones por item

3. compare_items(result1, result2, delta)
   - Compara dos items usando diferencias posteriores
   - Reporta: P(Item1 > Item2), HDI de la diferencia
   - Incluye: plot() para visualizar

4. compare_judge_groups(grupo1, grupo2, delta)
   - Compara dos paneles independientes de jueces
   - Reporta: P(Grupo1 > Grupo2), P(equivalencia)
   - Util para: detectar sesgos entre tipos de jueces

5. rope_analysis(result, rope_lower, rope_upper)
   - Analisis ROPE (Region of Practical Equivalence)
   - Decision basada en HDI + ROPE
   - Incluye: plot() para visualizar

6. sample_size_planning(result, target_prob, threshold)
   - Planificacion de tamano de muestra
   - Responde: '¿Cuantos jueces mas necesito?'
   - Incluye: plot() para visualizar curva
")

# =============================================================================
# SECCION 14: EJEMPLO INTEGRADO - FLUJO DE TRABAJO COMPLETO
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 14: FLUJO DE TRABAJO COMPLETO\n")
cat("================================================================\n")
cat("\nEjemplo de un flujo de analisis completo para un estudio de\n")
cat("validez de contenido:\n\n")

cat("PASO 1: Calcular V de Aiken bayesiano para todos los items\n")
cat("        resultado <- bayes_aiken(datos, coefficient = 'V', ...)\n\n")

cat("PASO 2: Aplicar regla de decision a cada item\n")
cat("        decisiones <- decision_rule_multi(resultado, threshold = 0.70)\n\n")

cat("PASO 3: Identificar items que necesitan revision\n")
cat("        items_revisar <- decisiones[decisiones$decision == 'INDECISO', ]\n\n")

cat("PASO 4: Para items indecisos, planificar jueces adicionales\n")
cat("        plan <- sample_size_planning(item_indeciso, target_prob = 0.95)\n\n")

cat("PASO 5: Comparar items problematicos con items buenos\n")
cat("        comp <- compare_items(item_bueno, item_problematico)\n\n")

cat("PASO 6: Si hay grupos de jueces, verificar equivalencia\n")
cat("        equiv <- compare_judge_groups(grupo1, grupo2, delta = 0.05)\n\n")

cat("PASO 7: Generar reporte final con visualizaciones\n")
cat("        plot_multiple_posteriors(resultados, n_cols = 5)\n\n")


# =============================================================================
# SECCION 15: EXTRACCION DE RESULTADOS EN DATA.FRAME
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECCION 15: EXTRACCION DE RESULTADOS (extract_results)\n")
cat("================================================================\n")

cat("\nLa funcion extract_results() permite extraer los coeficientes y sus\n")
cat("intervalos de credibilidad en un data.frame ordenado, util para:\n")
cat("- Exportar resultados a Excel/CSV\n")
cat("- Crear tablas para publicaciones\n")
cat("- Analisis adicionales\n\n")

# 15.1 Crear varios coeficientes de ejemplo
cat("--- 15.1 Crear coeficientes para varios items ---\n\n")

item1_V <- coef_V(c(3, 3, 2, 3, 3, 3, 3, 2, 3, 3), l = 0, s = 3, verbose = FALSE)
item2_V <- coef_V(c(3, 2, 3, 3, 2, 3, 2, 3, 3, 2), l = 0, s = 3, verbose = FALSE)
item3_V <- coef_V(c(2, 2, 3, 2, 3, 2, 3, 2, 2, 3), l = 0, s = 3, verbose = FALSE)
item4_V <- coef_V(c(3, 3, 3, 3, 3, 3, 3, 3, 3, 2), l = 0, s = 3, verbose = FALSE)
item5_V <- coef_V(c(2, 3, 2, 2, 3, 2, 2, 3, 2, 2), l = 0, s = 3, verbose = FALSE)

# 15.2 Extraer un solo resultado
cat("--- 15.2 Extraer un solo resultado ---\n\n")
df_single <- extract_results(item1_V)
print(df_single)

# 15.3 Extraer multiples resultados sin nombres
cat("\n--- 15.3 Extraer multiples resultados ---\n\n")
df_multi <- extract_results(item1_V, item2_V, item3_V)
print(df_multi)

# 15.4 Extraer con nombres personalizados
cat("\n--- 15.4 Con nombres personalizados ---\n\n")
df_named <- extract_results(
  Item1_Claridad = item1_V,
  Item2_Claridad = item2_V,
  Item3_Claridad = item3_V,
  Item4_Claridad = item4_V,
  Item5_Claridad = item5_V
)
print(df_named)

# 15.5 Incluyendo ambos intervalos (ETI y HDI)
cat("\n--- 15.5 Incluyendo ETI y HDI ---\n\n")
df_both <- extract_results(
  Item1 = item1_V,
  Item2 = item2_V,
  ci_type = "both"
)
print(df_both)

# 15.6 Desde una lista de resultados
cat("\n--- 15.6 Desde una lista ---\n\n")
lista_items <- list(
  Relevancia = item1_V,
  Representatividad = item2_V,
  Claridad = item3_V,
  Coherencia = item4_V,
  Suficiencia = item5_V
)
df_lista <- extract_results(lista_items)
print(df_lista)

# 15.7 Filtrar items problematicos
cat("\n--- 15.7 Filtrar items con P(V > 0.70) < 0.90 ---\n\n")
items_revisar <- df_lista[df_lista$P_70 < 0.90, ]
if (nrow(items_revisar) > 0) {
  cat("Items que requieren revision:\n")
  print(items_revisar)
} else {
  cat("Todos los items tienen P(V > 0.70) >= 0.90\n")
}

# 15.8 Exportar a CSV (ejemplo)
cat("\n--- 15.8 Exportar a CSV ---\n\n")
cat("Para exportar los resultados a CSV:\n")
cat("  write.csv(df_lista, 'resultados_validez.csv', row.names = FALSE)\n\n")

cat("\n")
cat("================================================================\n")
cat("           FIN DEL EJEMPLO REPRODUCIBLE\n")
cat("================================================================\n")
