# Leave-one-out cross-validation and permutation test for the OLS models.
# Outcome: MaxLevel.
#
# LOOCV is always applied to the FULL model (covariates + neural predictor as a
# single lm()). In each fold the full model is refit on the training data and the
# left-out subject is predicted with that model. 

# 0. Packages and shared style
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)
library(ggpubr)

source("style_REFUEL.R")   # theme_refuel(), CV palette, save_tiff(), etc.

REFUEL_SEED <- 2026

# 1. Load data (adjust the path to the data file on your machine)
df_raw <- read_excel("Hipp_integrity_predicts_Ind_dif.xlsx")

names(df_raw)

# 2. Prepare data
df_model <- df_raw %>%
  rename(
    Subject                 = SubjectID,
    MaxLevel                = Max_level,
    TBV                     = ses_pre_TBV,
    L_HC                    = ses_pre_lh_whole_HC,
    L_CA1                   = ses_pre_lh_CA1,
    L_CA23                  = ses_pre_lh_CA23,
    L_CA4DG                 = ses_pre_lh_CA4_DG,
    L_HC_tail               = ses_pre_lh_Hippocampal_tail,
    L_HC_body               = ses_pre_lh_Hippocampal_body,
    L_HC_head               = ses_pre_lh_Hippocampal_head,
    L_HC_fneurite           = fneurite_LH_HC_total,
    L_HC_anterior_fneurite  = fneurite_LH_HC_ant,
    L_HC_posterior_fneurite = fneurite_LH_HC_post
  )

cat("Total N (before dropping NA per model):", nrow(df_model), "\n")
str(df_model)

# 3. Model definitions (always covariates + predictor together)

modelos_cv <- list(
  list(nombre = "Left HC",                        predictor = "L_HC",      covariates = c("TBV", "Education")),
  list(nombre = "CA1",                            predictor = "L_CA1",     covariates = c("TBV", "Education")),
  list(nombre = "CA2/3",                          predictor = "L_CA23",    covariates = c("TBV", "Education")),
  list(nombre = "CA4-DG",                         predictor = "L_CA4DG",   covariates = c("TBV", "Education")),
  list(nombre = "Head",                           predictor = "L_HC_head", covariates = c("TBV", "Education")),
  list(nombre = "Body",                           predictor = "L_HC_body", covariates = c("TBV", "Education")),
  list(nombre = "Tail",                           predictor = "L_HC_tail", covariates = c("TBV", "Education")),
  list(nombre = "Left HC fneurite",               predictor = "L_HC_fneurite",           covariates = c("Education")),
  list(nombre = "Left HC fneurite (anterior)",    predictor = "L_HC_anterior_fneurite",  covariates = c("Education")),
  list(nombre = "Left HC fneurite (posterior)",   predictor = "L_HC_posterior_fneurite", covariates = c("Education"))
)

# 4. Helper functions

get_model_rmse <- function(model) {
  sqrt(mean(residuals(model)^2))
}

get_lm_coefficients <- function(model) {
  coef_table <- as.data.frame(summary(model)$coefficients)
  coef_table$Predictor <- rownames(coef_table)
  rownames(coef_table) <- NULL
  coef_table <- coef_table[, c("Predictor", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  names(coef_table) <- c("Predictor", "Estimate", "Std_Error", "t_value", "p_value")
  coef_table
}

# 5. Leave-one-out on the full model

loo_lm <- function(data, formula, id_col = NULL) {

  n <- nrow(data)
  outcome_name <- all.vars(formula)[1]
  observed <- data[[outcome_name]]

  predicted    <- rep(NA, n)
  residuals_cv <- rep(NA, n)

  for (i in 1:n) {
    train_data <- data[-i, ]
    test_data  <- data[i, ]

    model_i <- lm(formula, data = train_data)

    predicted[i]    <- predict(model_i, newdata = test_data)
    residuals_cv[i] <- observed[i] - predicted[i]
  }

  rmse_cv <- sqrt(mean(residuals_cv^2))
  mae_cv  <- mean(abs(residuals_cv))
  r2_cv   <- 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  cor_obs_pred <- cor(observed, predicted)

  results <- data.frame(
    observed    = observed,
    predicted   = predicted,
    residual_cv = residuals_cv,
    abs_error   = abs(residuals_cv)
  )

  if (!is.null(id_col)) {
    results <- cbind(Participant = data[[id_col]], results)
  }

  metrics <- data.frame(
    N = n,
    RMSE_CV = rmse_cv,
    MAE_CV = mae_cv,
    R2_CV = r2_cv,
    Cor_observed_predicted = cor_obs_pred
  )

  list(results = results, metrics = metrics)
}

# 6. Permutation test for the LOOCV R2 of the full model

perm_test_loocv <- function(data, formula, observed_r2cv, n_perm = 5000, seed = REFUEL_SEED) {

  set.seed(seed)

  outcome_name <- all.vars(formula)[1]
  perm_r2cv <- numeric(n_perm)
  perm_cor  <- numeric(n_perm)

  for (b in 1:n_perm) {
    data_perm <- data
    data_perm[[outcome_name]] <- sample(data_perm[[outcome_name]])

    loo_perm <- loo_lm(data = data_perm, formula = formula, id_col = NULL)

    perm_r2cv[b] <- loo_perm$metrics$R2_CV
    perm_cor[b]  <- loo_perm$metrics$Cor_observed_predicted

    if (b %% 500 == 0) cat("Permutation", b, "of", n_perm, "completed\n")
  }

  p_r2cv <- (sum(perm_r2cv >= observed_r2cv) + 1) / (n_perm + 1)

  list(
    observed_R2_CV       = observed_r2cv,
    permutation_p_R2_CV  = p_r2cv,
    permutation_R2_CV    = perm_r2cv,
    permutation_cor      = perm_cor,
    n_perm               = n_perm
  )
}

# 7. Run one model: in-sample + LOOCV + permutation

ejecutar_modelo_cv <- function(datos, modelo_cfg, outcome = "MaxLevel",
                               n_perm = 5000, seed = REFUEL_SEED) {

  vars_modelo  <- c("Subject", outcome, modelo_cfg$covariates, modelo_cfg$predictor)
  datos_modelo <- datos %>% dplyr::select(dplyr::all_of(vars_modelo)) %>% na.omit()

  formula_completa <- as.formula(
    paste(outcome, "~", paste(c(modelo_cfg$covariates, modelo_cfg$predictor), collapse = " + "))
  )

  fit_full <- lm(formula_completa, data = datos_modelo)

  insample <- data.frame(
    Model = modelo_cfg$nombre,
    N = nrow(datos_modelo),
    R2_in_sample = summary(fit_full)$r.squared,
    Adjusted_R2_in_sample = summary(fit_full)$adj.r.squared,
    RMSE_in_sample = get_model_rmse(fit_full)
  )

  loo <- loo_lm(data = datos_modelo, formula = formula_completa, id_col = "Subject")

  perm <- perm_test_loocv(
    data = datos_modelo, formula = formula_completa,
    observed_r2cv = loo$metrics$R2_CV,
    n_perm = n_perm, seed = seed
  )

  list(
    nombre   = modelo_cfg$nombre,
    formula  = formula_completa,
    fit_full = fit_full,
    insample = insample,
    loo      = loo,
    perm     = perm
  )
}

resultados_cv <- lapply(modelos_cv, ejecutar_modelo_cv, datos = df_model)
names(resultados_cv) <- sapply(modelos_cv, function(m) m$nombre)

# 8. Worst-predicted per model

worst_predictions <- lapply(resultados_cv, function(r) {
  r$loo$results %>% arrange(desc(abs_error))
})

# 9. LOOCV figures (purple/green palette, statistics reported in table/caption)
# Panel grouping: anatomical (CA1/CA2-3/CA4-DG), A-P axis (Head/Body/Tail) and
# fneurite axis (anterior/posterior).

nombres_anatomico    <- c("CA1", "CA2/3", "CA4-DG")
nombres_eje          <- c("Head", "Body", "Tail")
nombres_fneurite_eje <- c("Left HC fneurite (anterior)", "Left HC fneurite (posterior)")
nombres_individuales <- setdiff(names(resultados_cv),
                                c(nombres_anatomico, nombres_eje,
                                  nombres_fneurite_eje))

make_loocv_scatter <- function(resultado_cv, panel = FALSE, axis_limits = NULL) {
  df_plot <- resultado_cv$loo$results

  p <- ggplot(df_plot, aes(x = observed, y = predicted)) +
    geom_point(color = REFUEL_COL_CV_POINT, size = REFUEL_POINT_SIZE, alpha = REFUEL_POINT_ALPHA) +
    geom_smooth(method = "lm", se = TRUE,
                color = REFUEL_COL_CV_LINE, fill = REFUEL_COL_CV_BAND,
                alpha = 0.15, linewidth = REFUEL_SMOOTH_LWD) +
    labs(
      title    = format_fneurite_title(paste0("LOOCV — ", resultado_cv$nombre), bold = TRUE),
      subtitle = if (panel) NULL else paste(deparse(resultado_cv$formula), collapse = " "),
      x        = "Observed MaxLevel",
      y        = "Predicted MaxLevel"
    ) +
    stat_annotation("italic(R)[CV]^2",
                    resultado_cv$loo$metrics$R2_CV,
                    resultado_cv$perm$permutation_p_R2_CV) +
    theme_refuel(base_size = if (panel) 12 else 14,
                 axis_title_size = if (panel) 11 else 13,
                 axis_text_size  = if (panel) 9  else 11)

  if (!is.null(axis_limits)) {
    p <- p + coord_cartesian(xlim = axis_limits, ylim = axis_limits)
  }

  p
}

make_perm_histogram <- function(resultado_cv, panel = FALSE) {
  ggplot(data.frame(R2_CV_perm = resultado_cv$perm$permutation_R2_CV), aes(x = R2_CV_perm)) +
    geom_histogram(bins = 40, color = "white", fill = REFUEL_COL_CV_HIST, alpha = 0.75) +
    geom_vline(xintercept = resultado_cv$perm$observed_R2_CV,
               color = REFUEL_COL_CV_OBS, linewidth = 1, linetype = "dashed") +
    labs(
      title = format_fneurite_title(paste0("Permutation test — ", resultado_cv$nombre), bold = TRUE),
      x     = "Permuted LOOCV R²",
      y     = "Count"
    ) +
    theme_refuel(base_size = if (panel) 11 else 12,
                 axis_title_size = if (panel) 10 else 11,
                 axis_text_size  = if (panel) 9  else 10)
}

axis_limits_grupo <- function(nombres) {
  valores <- unlist(lapply(resultados_cv[nombres], function(r) {
    c(r$loo$results$observed, r$loo$results$predicted)
  }))
  pad_range(range(valores))
}

axis_limits_anat         <- axis_limits_grupo(nombres_anatomico)
axis_limits_eje          <- axis_limits_grupo(nombres_eje)
axis_limits_fneurite_eje <- axis_limits_grupo(nombres_fneurite_eje)

# 9a. Individual models (Left HC, whole fneurite)
for (nombre in nombres_individuales) {
  res  <- resultados_cv[[nombre]]
  slug <- gsub("[^A-Za-z0-9]+", "_", nombre)

  save_tiff(make_loocv_scatter(res), paste0("LOOCV_", slug, ".tiff"),
            width = 5, height = 4.5, dpi = 600)

  save_tiff(make_perm_histogram(res), paste0("Permutation_", slug, ".tiff"),
            width = 5, height = 4.5, dpi = 600)
}

# 9b. Anatomical panel (CA1, CA2/3, CA4-DG)
plots_loocv_anat <- lapply(resultados_cv[nombres_anatomico], make_loocv_scatter,
                           panel = TRUE, axis_limits = axis_limits_anat)
panel_loocv_anat <- ggarrange(plotlist = plots_loocv_anat, ncol = 3, nrow = 1)
panel_loocv_anat <- annotate_figure(panel_loocv_anat,
                                    top = text_grob("LOOCV — Anatomical subregions — Left hippocampus",
                                                    face = "bold", size = 13))
save_tiff(panel_loocv_anat, "LOOCV_panel_anatomico.tiff", width = 13, height = 5, dpi = 600)

plots_perm_anat <- lapply(resultados_cv[nombres_anatomico], make_perm_histogram, panel = TRUE)
panel_perm_anat <- ggarrange(plotlist = plots_perm_anat, ncol = 3, nrow = 1)
panel_perm_anat <- annotate_figure(panel_perm_anat,
                                   top = text_grob("Permutation test — Anatomical subregions — Left hippocampus",
                                                   face = "bold", size = 13))
save_tiff(panel_perm_anat, "Permutation_panel_anatomico.tiff", width = 13, height = 5, dpi = 600)

# 9c. Antero-posterior axis panel (Head, Body, Tail)
plots_loocv_eje <- lapply(resultados_cv[nombres_eje], make_loocv_scatter,
                          panel = TRUE, axis_limits = axis_limits_eje)
panel_loocv_eje <- ggarrange(plotlist = plots_loocv_eje, ncol = 3, nrow = 1)
panel_loocv_eje <- annotate_figure(panel_loocv_eje,
                                   top = text_grob("LOOCV — Antero-posterior axis — Left hippocampus",
                                                   face = "bold", size = 13))
save_tiff(panel_loocv_eje, "LOOCV_panel_eje.tiff", width = 13, height = 5, dpi = 600)

plots_perm_eje <- lapply(resultados_cv[nombres_eje], make_perm_histogram, panel = TRUE)
panel_perm_eje <- ggarrange(plotlist = plots_perm_eje, ncol = 3, nrow = 1)
panel_perm_eje <- annotate_figure(panel_perm_eje,
                                  top = text_grob("Permutation test — Antero-posterior axis — Left hippocampus",
                                                  face = "bold", size = 13))
save_tiff(panel_perm_eje, "Permutation_panel_eje.tiff", width = 13, height = 5, dpi = 600)

# 9d. fneurite antero-posterior axis panel (anterior, posterior)
plots_loocv_fneurite_eje <- lapply(resultados_cv[nombres_fneurite_eje], make_loocv_scatter,
                                   panel = TRUE, axis_limits = axis_limits_fneurite_eje)
panel_loocv_fneurite_eje <- ggarrange(plotlist = plots_loocv_fneurite_eje, ncol = 2, nrow = 1)
panel_loocv_fneurite_eje <- annotate_figure(panel_loocv_fneurite_eje,
                                            top = grid::textGrob(
                                              format_fneurite_title("LOOCV — Left HC fneurite — anterior/posterior", bold = TRUE),
                                              gp = grid::gpar(fontsize = 13)
                                            ))
save_tiff(panel_loocv_fneurite_eje, "LOOCV_panel_fneurite_eje.tiff", width = 9, height = 5, dpi = 600)

plots_perm_fneurite_eje <- lapply(resultados_cv[nombres_fneurite_eje], make_perm_histogram, panel = TRUE)
panel_perm_fneurite_eje <- ggarrange(plotlist = plots_perm_fneurite_eje, ncol = 2, nrow = 1)
panel_perm_fneurite_eje <- annotate_figure(panel_perm_fneurite_eje,
                                           top = grid::textGrob(
                                             format_fneurite_title("Permutation test — Left HC fneurite — anterior/posterior", bold = TRUE),
                                             gp = grid::gpar(fontsize = 13)
                                           ))
save_tiff(panel_perm_fneurite_eje, "Permutation_panel_fneurite_eje.tiff", width = 9, height = 5, dpi = 600)

# 10. Summary tables
insample_comparison <- bind_rows(lapply(resultados_cv, function(r) r$insample))

cv_comparison <- bind_rows(lapply(resultados_cv, function(r) {
  r$loo$metrics %>% mutate(Model = r$nombre) %>% select(Model, everything())
}))

full_model_comparison <- insample_comparison %>%
  left_join(cv_comparison, by = c("Model", "N"))

permutation_summary <- bind_rows(lapply(resultados_cv, function(r) {
  data.frame(
    Model = r$nombre,
    R2_CV = r$loo$metrics$R2_CV,
    RMSE_CV = r$loo$metrics$RMSE_CV,
    MAE_CV = r$loo$metrics$MAE_CV,
    Observed_predicted_r = r$loo$metrics$Cor_observed_predicted,
    Permutation_p_R2_CV = r$perm$permutation_p_R2_CV,
    N_permutations = r$perm$n_perm
  )
}))

print(full_model_comparison)
print(permutation_summary)

# 11. Save everything to Excel (one sheet per model where applicable)
generar_slugs <- function(nombres, prefix) {
  slugs <- character(length(nombres))
  for (i in seq_along(nombres)) {
    slugs[i] <- sheet_slug(nombres[i], prefix = prefix, existing = slugs[seq_len(i - 1)])
  }
  slugs
}

loo_sheets   <- setNames(lapply(resultados_cv, function(r) r$loo$results),
                         generar_slugs(names(resultados_cv), prefix = "LOOCV_"))
worst_sheets <- setNames(worst_predictions,
                         generar_slugs(names(worst_predictions), prefix = "Worst_"))
coef_sheets  <- setNames(lapply(resultados_cv, function(r) get_lm_coefficients(r$fit_full)),
                         generar_slugs(names(resultados_cv), prefix = "Coef_"))

write_xlsx(
  c(
    list(
      Full_model_comparison = full_model_comparison,
      CV_model_comparison   = cv_comparison,
      Permutation_summary   = permutation_summary,
      In_sample_comparison  = insample_comparison
    ),
    loo_sheets,
    worst_sheets,
    coef_sheets
  ),
  path = "LOOCV_results_left_HC_all_models.xlsx"
)

# 12. Final message
cat("\nLOOCV finished for", length(resultados_cv), "model(s).\n")
cat("Results saved to: LOOCV_results_left_HC_all_models.xlsx\n")
cat("Figures saved as TIFF (600 dpi, LZW):\n")
for (nombre in nombres_individuales) {
  slug <- gsub("[^A-Za-z0-9]+", "_", nombre)
  cat("- LOOCV_", slug, ".tiff / Permutation_", slug, ".tiff\n", sep = "")
}
cat("- LOOCV_panel_anatomico.tiff / Permutation_panel_anatomico.tiff (CA1, CA2/3, CA4-DG)\n")
cat("- LOOCV_panel_eje.tiff / Permutation_panel_eje.tiff (Head, Body, Tail)\n")
cat("- LOOCV_panel_fneurite_eje.tiff / Permutation_panel_fneurite_eje.tiff (fneurite anterior, posterior)\n")
