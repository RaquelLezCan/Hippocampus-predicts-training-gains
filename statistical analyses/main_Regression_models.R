# OLS regressions of the left hippocampus predicting MaxLevel.

library(readxl)
library(dplyr)
library(officer)
library(flextable)
library(ggplot2)
library(ggpubr)
source("style_REFUEL.R")   # theme_refuel(), palette, save_tiff(), sheet_slug()

# 0. Load and prepare data -----------------------------------------------------

datos <- read_excel("Hipp_integrity_predicts_Ind_dif.xlsx")

datos <- datos %>%
  rename(
    MaxLevel  = Max_level,
    HC_L      = ses_pre_lh_whole_HC,
    TBV       = ses_pre_TBV,
    CA1_L     = ses_pre_lh_CA1,
    CA23_L    = ses_pre_lh_CA23,
    CA4DG_L   = ses_pre_lh_CA4_DG,
    SUB_L     = ses_pre_lh_subiculum,
    PRESUB_L  = ses_pre_lh_presubiculum,
    FIMBRIA_L = ses_pre_lh_fimbria,
    HC_tail_L = ses_pre_lh_Hippocampal_tail,
    HC_body_L = ses_pre_lh_Hippocampal_body,
    HC_head_L = ses_pre_lh_Hippocampal_head,
    HC_R      = ses_pre_rh_whole_HC,
    LH_HC_fneurite           = fneurite_LH_HC_total,
    LH_HC_Anterior_fneurite  = fneurite_LH_HC_ant,
    LH_HC_Posterior_fneurite = fneurite_LH_HC_post,
    RH_HC_fneurite           = fneurite_RH_HC_total,
    LH_HC_fsoma              = fsoma_LH_HC_total,
    RH_HC_fsoma              = fsoma_RH_HC_total
  )

# 0b. Covariate screening with Spearman ----------------------------------------
# MaxLevel against sociodemographic and pre-training neuropsychological variables.
# Those significant (p < .05, uncorrected) are considered covariate candidates.
# Sex is binary, so its rho is equivalent to a Mann-Whitney; included for uniform
# criteria, to be kept in mind when writing up.

vars_sociodemo  <- c("Age", "Sex")

vars_neuropsico <- c(
  "Education",
  "ses_pre_DIGIT_forward", "ses_pre_DIGIT_backward",
  "ses_pre_CORSI_forward", "ses_pre_CORSI_backward",
  "ses_pre_REY_AVG", "ses_pre_REY_DELAY",
  "ses_pre_LETTCOMP_ACC", "ses_pre_PATTCOMP_ACC",
  "ses_pre_DIGIT_SYMBOL_ACC", "ses_pre_DOTMATRIX_ACC",
  "ses_pre_GDS",
  "ses_pre_Somatic_GAS", "ses_pre_Cognitive_GAS", "ses_pre_Affective_GAS", "ses_pre_GAS_total",
  "ses_pre_Physical_fun_SF_36", "ses_pre_Role_lim_physical_health_SF_36",
  "ses_pre_Role_lim_emoti_problems_SF_36", "ses_pre_Energy_fatigue_SF_36",
  "ses_pre_Emotional_wellbeing_SF_36", "ses_pre_Social_funct_SF_36",
  "ses_pre_Pain_SF_36", "ses_pre_General_health_SF_36"
)

vars_candidatas <- c(vars_sociodemo, vars_neuropsico)

correlacion_spearman <- function(datos, outcome, vars_candidatas) {
  resultados <- lapply(vars_candidatas, function(v) {
    x <- datos[[v]]
    y <- datos[[outcome]]

    # If the variable comes as text or factor, convert to numeric
    if (is.character(x) || is.factor(x)) {
      x <- as.numeric(as.factor(x))
    }

    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))

    data.frame(
      Variable = v,
      rho      = round(unname(ct$estimate), 3),
      S        = round(unname(ct$statistic), 1),
      p        = round(ct$p.value, 4),
      Sig      = ifelse(ct$p.value < 0.05, "*", ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, resultados)
}

tabla_corr_spearman <- correlacion_spearman(datos, "MaxLevel", vars_candidatas)

cat("\n========================================\n")
cat("Spearman correlations: MaxLevel vs. candidate covariates\n")
cat("========================================\n")
print(tabla_corr_spearman)

covariables_significativas <- tabla_corr_spearman$Variable[tabla_corr_spearman$Sig == "*"]
cat("\nVariables correlated with MaxLevel (p < .05):\n")
print(covariables_significativas)

fmt_corr_ft <- function(df) {
  flextable(df) |>
    bold(part = "header") |>
    bg(i = ~ Sig == "*", bg = "#FFF3CD") |>
    autofit()
}

# 1. Subregions and labels -----------------------------------------------------

subregiones_completo   <- c("HC_L")
subregiones_anatomicas <- c("CA1_L", "CA23_L", "CA4DG_L")
subregiones_subfields  <- c("SUB_L", "PRESUB_L", "FIMBRIA_L")
subregiones_eje        <- c("HC_head_L", "HC_body_L", "HC_tail_L")

# SANDI fneurite: left hemisphere only, without TBV
subregiones_sandi      <- c("LH_HC_fneurite")
subregiones_sandi_eje  <- c("LH_HC_Anterior_fneurite", "LH_HC_Posterior_fneurite")

# SANDI fsoma: same pattern as fneurite, without anterior/posterior axis (not
# significant in either hemisphere, so it does not enter the cross-validation)
subregiones_fsoma      <- c("LH_HC_fsoma")
subregiones_fsoma_RH   <- c("RH_HC_fsoma")

# Right hippocampus: reported to document the absence of effect
subregiones_completo_RH <- c("HC_R")
subregiones_sandi_RH    <- c("RH_HC_fneurite")

etiquetas <- c(
  HC_L                     = "Left HC",
  CA1_L                    = "CA1",
  CA23_L                   = "CA2/3",
  CA4DG_L                  = "CA4-DG",
  SUB_L                    = "Subiculum",
  PRESUB_L                 = "Presubiculum",
  FIMBRIA_L                = "Fimbria",
  HC_head_L                = "Head",
  HC_body_L                = "Body",
  HC_tail_L                = "Tail",
  LH_HC_fneurite           = "Left HC fneurite",
  LH_HC_Anterior_fneurite  = "Left HC fneurite (anterior)",
  LH_HC_Posterior_fneurite = "Left HC fneurite (posterior)",
  LH_HC_fsoma              = "Left HC fsoma",
  RH_HC_fsoma              = "Right HC fsoma",
  HC_R                     = "Right HC",
  RH_HC_fneurite           = "Right HC fneurite"
)

term_labels <- c(
  etiquetas,
  "(Intercept)"        = "Intercept",
  "TBV"                = "TBV",
  "Education"          = "Education",
  "ses_pre_GDS"           = "GDS",
  "ses_pre_DOTMATRIX_ACC" = "Dot Matrix (accuracy)"
)

# 2. OLS fitting function ------------------------------------------------------

run_ols <- function(datos, outcome, predictor, covariates) {

  formula_obj <- as.formula(
    paste(outcome, "~", paste(c(covariates, predictor), collapse = " + "))
  )

  fit      <- lm(formula_obj, data = datos)
  s        <- summary(fit)
  ci_all   <- confint(fit)
  coef_tab <- s$coefficients

  # Standardized betas for all terms, to compare effect sizes across covariates
  # with very different scales.

  vars_std  <- c(outcome, covariates, predictor)
  datos_std <- datos
  datos_std[vars_std] <- scale(datos_std[vars_std])
  fit_std   <- lm(formula_obj, data = datos_std)
  beta_std  <- unname(coef(fit_std)[predictor])
  beta_all  <- coef(fit_std)
  beta_all["(Intercept)"] <- NA_real_

  full_table <- data.frame(
    Term        = rownames(coef_tab),
    b           = round(coef_tab[, "Estimate"], 3),
    Beta        = round(unname(beta_all[rownames(coef_tab)]), 3),
    SE          = round(coef_tab[, "Std. Error"], 3),
    CI_low      = round(ci_all[, 1], 3),
    CI_high     = round(ci_all[, 2], 3),
    t           = round(coef_tab[, "t value"], 3),
    p           = round(coef_tab[, "Pr(>|t|)"], 4),
    row.names   = NULL,
    check.names = FALSE
  )
  full_table$Term <- ifelse(full_table$Term %in% names(term_labels),
                            term_labels[full_table$Term],
                            full_table$Term)
  full_table$Sig <- ifelse(full_table$p < 0.05, "*", "")

  ci <- ci_all[predictor, ]

   
  model_F    <- unname(s$fstatistic["value"])
  model_df1  <- unname(s$fstatistic["numdf"])
  model_df2  <- unname(s$fstatistic["dendf"])
  model_p    <- pf(model_F, model_df1, model_df2, lower.tail = FALSE)

  list(
    predictor  = predictor,
    label      = etiquetas[predictor],
    b          = round(coef_tab[predictor, "Estimate"], 3),
    beta_std   = round(beta_std, 3),
    se         = round(coef_tab[predictor, "Std. Error"], 3),
    ci_low     = round(ci[1], 3),
    ci_high    = round(ci[2], 3),
    t          = round(coef_tab[predictor, "t value"], 3),
    df         = fit$df.residual,
    p_raw      = round(coef_tab[predictor, "Pr(>|t|)"], 4),
    R2         = round(s$r.squared, 3),
    R2adj      = round(s$adj.r.squared, 3),
    model_F    = round(model_F, 3),
    model_df1  = model_df1,
    model_df2  = model_df2,
    model_p    = round(model_p, 4),
    full_table = full_table,
    fit        = fit
  )
}

# 3. Run the analyses ----------------------------------------------------------

correr_analisis <- function(datos, subregiones, outcome = "MaxLevel",
                            covariates = c("TBV", "Education"),
                            nombre_analisis = "") {

  resultados <- lapply(subregiones, function(pred) {
    run_ols(datos, outcome, pred, covariates)
  })

  for (i in seq_along(resultados)) {
    resultados[[i]]$sig <- ifelse(resultados[[i]]$p_raw < 0.05, "*", "")
  }

  cat("\n========================================\n")
  cat("Analysis:", nombre_analisis, "\n")
  cat("Uncorrected p (exploratory hierarchical follow-up)\n")
  cat("========================================\n")

  for (res in resultados) {
    cat(sprintf("\nPredictor: %s\n", res$predictor))
    cat(sprintf("  b = %.3f [%.3f, %.3f], t(%d) = %.3f\n",
                res$b, res$ci_low, res$ci_high, res$df, res$t))
    cat(sprintf("  p = %.4f %s\n", res$p_raw, res$sig))
    cat(sprintf("  R2 = %.3f | R2adj = %.3f\n", res$R2, res$R2adj))
    cat(sprintf("  Model F(%d, %d) = %.3f, p = %.4f\n",
                res$model_df1, res$model_df2, res$model_F, res$model_p))
  }

  resultados
}

res_completo  <- correr_analisis(datos, subregiones_completo,
                                 nombre_analisis = "Whole left HC")
res_anatomico <- correr_analisis(datos, subregiones_anatomicas,
                                 nombre_analisis = "Anatomical subregions (CA1, CA2/3, CA4-DG)")
res_subfields <- correr_analisis(datos, subregiones_subfields,
                                 nombre_analisis = "Additional subregions (subiculum, presubiculum, fimbria)")
res_eje       <- correr_analisis(datos, subregiones_eje,
                                 nombre_analisis = "Antero-posterior axis (head, body, tail)")

res_sandi     <- correr_analisis(datos, subregiones_sandi,
                                 covariates = c("Education"),
                                 nombre_analisis = "SANDI fneurite (left HC)")
res_sandi_eje <- correr_analisis(datos, subregiones_sandi_eje,
                                 covariates = c("Education"),
                                 nombre_analisis = "SANDI fneurite anterior/posterior")

res_fsoma    <- correr_analisis(datos, subregiones_fsoma,
                                covariates = c("Education"),
                                nombre_analisis = "SANDI fsoma (left HC)")
res_fsoma_RH <- correr_analisis(datos, subregiones_fsoma_RH,
                                covariates = c("Education"),
                                nombre_analisis = "SANDI fsoma (right HC)")

res_completo_RH <- correr_analisis(datos, subregiones_completo_RH,
                                   nombre_analisis = "Whole right HC (volume)")
res_sandi_RH    <- correr_analisis(datos, subregiones_sandi_RH,
                                   covariates = c("Education"),
                                   nombre_analisis = "SANDI fneurite (right HC)")

# 3b. Reference cognitive model (hierarchical step 1, no neural predictor) ------
# Three predictors: GDS (kept on theoretical grounds following Engvig et al.,
# 2012, despite not passing the screening), Education (the only significant one)
# and Dot Matrix accuracy. Education is passed as the focal predictor in run_ols(); 
# GDS and Dot Matrix remain in full_table.

res_cognitivo <- run_ols(datos, "MaxLevel", "Education",
                         covariates = c("ses_pre_GDS", "ses_pre_DOTMATRIX_ACC"))
res_cognitivo$label <- "Cognitive model (GDS + Education + Dot Matrix)"
res_cognitivo$sig   <- ifelse(res_cognitivo$p_raw < 0.05, "*", "")

cat("\n========================================\n")
cat("Reference cognitive model (step 1)\n")
cat("MaxLevel ~ GDS + Education + Dot Matrix accuracy\n")
cat("========================================\n")
cat(sprintf("\n%s\n", res_cognitivo$label))
cat(sprintf("  Model F(%d, %d) = %.3f, p = %.4f\n",
            res_cognitivo$model_df1, res_cognitivo$model_df2,
            res_cognitivo$model_F, res_cognitivo$model_p))
cat("\nFull model coefficients:\n")
print(res_cognitivo$full_table)

# 4. Common y-axis limits (residualized MaxLevel) ------------------------------
y_limits_vol   <- pad_range(range(residuals(lm(MaxLevel ~ TBV + Education, data = datos))))
y_limits_sandi <- pad_range(range(residuals(lm(MaxLevel ~ Education,       data = datos))))

# 5. Residualized scatter plots, with beta and p annotation --------------------


axis_x_label <- function(label) {
  if (grepl("fneurite", label, fixed = TRUE)) {
    format_fneurite_title(label, bold = FALSE)
  } else {
    format_fneurite_title(paste0(label, " (mm³)"), bold = FALSE)
  }
}

make_scatter <- function(datos, resultado, outcome = "MaxLevel",
                         covariates = c("TBV", "Education"), panel = FALSE,
                         y_limits = NULL) {

  pred  <- resultado$predictor
  label <- resultado$label

  form_cov <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  form_pre <- as.formula(paste(pred,    "~", paste(covariates, collapse = " + ")))

  resid_y <- residuals(lm(form_cov, data = datos))
  resid_x <- residuals(lm(form_pre, data = datos))

  df_scatter <- data.frame(x = resid_x, y = resid_y)

  p <- ggplot(df_scatter, aes(x = x, y = y)) +
    geom_point(color = REFUEL_COL_REG_POINT, size = REFUEL_POINT_SIZE, alpha = REFUEL_POINT_ALPHA) +
    geom_smooth(method = "lm", se = TRUE, color = REFUEL_COL_REG_LINE,
                fill = REFUEL_COL_REG_LINE, alpha = 0.15, linewidth = REFUEL_SMOOTH_LWD) +
    labs(
      title = format_fneurite_title(label, bold = TRUE),
      x     = axis_x_label(label),
      y     = "MaxLevel (residualized)"
    ) +
    stat_annotation("italic(beta)", resultado$beta_std, resultado$p_raw) +
    theme_refuel(base_size = if (panel) 12 else 14,
                 axis_title_size = if (panel) 11 else 13,
                 axis_text_size  = if (panel) 9  else 11)

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  p
}

p_completo <- make_scatter(datos, res_completo[[1]], y_limits = y_limits_vol)
save_tiff(p_completo, "scatter_HC_completo.tiff", width = 5, height = 4.5, dpi = 600)

plots_anat <- lapply(res_anatomico, make_scatter, datos = datos, panel = TRUE,
                     y_limits = y_limits_vol)
panel_anat <- ggarrange(plotlist = plots_anat, ncol = 3, nrow = 1)
panel_anat <- annotate_figure(panel_anat,
                              top = text_grob("Anatomical subregions — Left hippocampus",
                                              face = "bold", size = 13))
save_tiff(panel_anat, "scatter_anatomico.tiff", width = 13, height = 5, dpi = 600)

plots_subfields <- lapply(res_subfields, make_scatter, datos = datos, panel = TRUE,
                          y_limits = y_limits_vol)
panel_subfields <- ggarrange(plotlist = plots_subfields, ncol = 3, nrow = 1)
panel_subfields <- annotate_figure(panel_subfields,
                              top = text_grob("Additional subregions (subiculum, presubiculum, fimbria) — Left hippocampus",
                                              face = "bold", size = 13))
save_tiff(panel_subfields, "scatter_subfields.tiff", width = 13, height = 5, dpi = 600)

plots_eje <- lapply(res_eje, make_scatter, datos = datos, panel = TRUE,
                    y_limits = y_limits_vol)
panel_eje <- ggarrange(plotlist = plots_eje, ncol = 3, nrow = 1)
panel_eje <- annotate_figure(panel_eje,
                             top = text_grob("Antero-posterior axis — Left hippocampus",
                                             face = "bold", size = 13))
save_tiff(panel_eje, "scatter_eje.tiff", width = 13, height = 5, dpi = 600)

p_sandi <- make_scatter(datos, res_sandi[[1]], covariates = c("Education"),
                        y_limits = y_limits_sandi)
save_tiff(p_sandi, "scatter_SANDI_fneurite.tiff", width = 5, height = 4.5, dpi = 600)

p_fsoma <- make_scatter(datos, res_fsoma[[1]], covariates = c("Education"),
                        y_limits = y_limits_sandi)
save_tiff(p_fsoma, "scatter_SANDI_fsoma.tiff", width = 5, height = 4.5, dpi = 600)

plots_sandi_eje <- lapply(res_sandi_eje, make_scatter, datos = datos,
                          covariates = c("Education"), panel = TRUE,
                          y_limits = y_limits_sandi)
panel_sandi_eje <- ggarrange(plotlist = plots_sandi_eje, ncol = 2, nrow = 1)
panel_sandi_eje <- annotate_figure(panel_sandi_eje,
                                   top = grid::textGrob(
                                     format_fneurite_title("Left HC fneurite — anterior/posterior", bold = TRUE),
                                     gp = grid::gpar(fontsize = 13)
                                   ))
save_tiff(panel_sandi_eje, "scatter_SANDI_eje.tiff", width = 9, height = 5, dpi = 600)

p_completo_RH <- make_scatter(datos, res_completo_RH[[1]], y_limits = y_limits_vol)
save_tiff(p_completo_RH, "scatter_HC_completo_RH.tiff", width = 5, height = 4.5, dpi = 600)

p_sandi_RH <- make_scatter(datos, res_sandi_RH[[1]], covariates = c("Education"),
                           y_limits = y_limits_sandi)
save_tiff(p_sandi_RH, "scatter_SANDI_fneurite_RH.tiff", width = 5, height = 4.5, dpi = 600)

p_fsoma_RH <- make_scatter(datos, res_fsoma_RH[[1]], covariates = c("Education"),
                           y_limits = y_limits_sandi)
save_tiff(p_fsoma_RH, "scatter_SANDI_fsoma_RH.tiff", width = 5, height = 4.5, dpi = 600)

# 6. Tables and Word export ----------------------------------------------------

construir_tabla_resumen <- function(resultados) {
  do.call(rbind, lapply(resultados, function(r) {
    data.frame(
      Predictor  = r$label,
      b          = r$b,
      SE         = r$se,
      "IC 95%"   = paste0("[", r$ci_low, ", ", r$ci_high, "]"),
      t          = r$t,
      df         = r$df,
      p          = r$p_raw,
      Sig        = r$sig,
      R2adj      = r$R2adj,
      "Model F"  = r$model_F,
      "Model df" = paste0(r$model_df1, ", ", r$model_df2),
      "Model p"  = r$model_p,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }))
}

fmt_resumen_ft <- function(df) {
  flextable(df) |>
    bold(part = "header") |>
    bg(i = ~ Sig == "*", bg = "#FFF3CD") |>
    autofit()
}

fmt_full_ft <- function(full_table) {
  flextable(full_table[, c("Term", "b", "Beta", "SE", "CI_low", "CI_high", "t", "p", "Sig")]) |>
    set_header_labels(
      Term = "Predictor", b = "b", Beta = "β (std.)", SE = "SE",
      CI_low = "CI 95% (low)", CI_high = "CI 95% (high)",
      t = "t", p = "p", Sig = ""
    ) |>
    bold(part = "header") |>
    bg(i = ~ Sig == "*", bg = "#FFF3CD") |>
    autofit()
}

agregar_modelos_word <- function(doc, resultados, etiqueta_modelo, modelo_txt) {
  for (res in resultados) {
    doc <- doc |>
      officer::body_add_par(paste0(etiqueta_modelo, " — ", res$label), style = "heading 3") |>
      officer::body_add_par(modelo_txt) |>
      flextable::body_add_flextable(fmt_full_ft(res$full_table)) |>
      officer::body_add_par(
        sprintf("Model fit: F(%d, %d) = %.3f, p = %.4f | R² = %.3f, R²adj = %.3f",
                res$model_df1, res$model_df2, res$model_F, res$model_p,
                res$R2, res$R2adj)
      ) |>
      officer::body_add_par("")
  }
  doc
}

doc <- officer::read_docx() |>

  officer::body_add_par("Supplementary Table S0. Spearman correlations — MaxLevel vs. candidate covariates", style = "heading 2") |>
  officer::body_add_par("Sociodemographic (Age, Sex) and pre-training neuropsychological measures. Significant variables (p < .05, uncorrected) considered as covariate candidates for the regression models.") |>
  flextable::body_add_flextable(fmt_corr_ft(tabla_corr_spearman)) |>
  officer::body_add_par("") |>

  officer::body_add_par("Table 0. Cognitive model — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ GDS + Education + Dot Matrix accuracy. OLS. Step 1 of the hierarchical design (sociodemographic/cognitive/affective predictors). Row shown is the Education term, the only significant predictor and the covariate retained for the neural models (Step 2).") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(list(res_cognitivo)))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, list(res_cognitivo), "Cognitive model — full model",
                            "Model: MaxLevel ~ GDS + Education + Dot Matrix accuracy. OLS. All three coefficients shown (GDS and Dot Matrix accuracy included for theoretical/empirical completeness; only Education reached significance).")

doc <- doc |>
  officer::body_add_par("Table 1. Left hippocampus — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ TBV + Education + HC_L. OLS. ") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_completo))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_completo, "Left HC — full model",
                            "Model: MaxLevel ~ TBV + Education + HC_L. OLS.")

doc <- doc |>
  officer::body_add_par("Table 2. Anatomical subregions — Left hippocampus — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hierarchical follow-up.") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_anatomico))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_anatomico, "Anatomical subregion — full model",
                            "Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hierarchical follow-up.")

doc <- doc |>
  officer::body_add_par("Table 2b. Additional subregions (subiculum, presubiculum, fimbria) — Left hippocampus — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hierarchical follow-up.") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_subfields))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_subfields, "Additional subregion — full model",
                            "Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hierarchical follow-up.")

doc <- doc |>
  officer::body_add_par("Table 3. Antero-posterior axis — Left hippocampus — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hhierarchical follow-up.") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_eje))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_eje, "Antero-posterior subregion — full model",
                            "Model: MaxLevel ~ TBV + Education + HC_subregion. OLS. Hierarchical follow-up).")

doc <- doc |>
  officer::body_add_par("Table 4. SANDI fneurite (microstructural) — Left HC — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ Education + fneurite. OLS.") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_sandi))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_sandi, "SANDI fneurite Left HC — full model",
                            "Model: MaxLevel ~ Education + fneurite. OLS.")

doc <- doc |>
  officer::body_add_par("Table 5. SANDI fneurite — Left HC anterior/posterior — summary", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ Education + fneurite_subregion. OLS. Hierarchical follow-up to whole-HC SANDI model.") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_sandi_eje))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_sandi_eje, "SANDI fneurite subregion — full model",
                            "Model: MaxLevel ~ Education + fneurite_subregion. OLS.")

doc <- doc |>
  officer::body_add_par("Table 6. Right hippocampus volume — null result", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ TBV + Education + HC_R. OLS. Reported to document absence of effect in RH (no subregion follow-up, no LOOCV).") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_completo_RH))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_completo_RH, "Right HC volume — full model",
                            "Model: MaxLevel ~ TBV + Education + HC_R. OLS.")

doc <- doc |>
  officer::body_add_par("Table 7. Right hippocampus SANDI fneurite", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ Education + fneurite. OLS. Reported to document absence of effect in RH (no subregion follow-up, no LOOCV).") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_sandi_RH))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_sandi_RH, "Right HC SANDI fneurite — full model",
                            "Model: MaxLevel ~ Education + fneurite. OLS.")

doc <- doc |>
  officer::body_add_par("Table 8. SANDI fsoma — Left HC — null result", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ Education + fsoma. OLS. Reported to document absence of effect (no anterior/posterior follow-up, no LOOCV).") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_fsoma))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_fsoma, "SANDI fsoma Left HC — full model",
                            "Model: MaxLevel ~ Education + fsoma. OLS.")

doc <- doc |>
  officer::body_add_par("Table 9. SANDI fsoma — Right HC — null result", style = "heading 2") |>
  officer::body_add_par("Model: MaxLevel ~ Education + fsoma. OLS. TBV NOT used as covariate. Reported to document absence of effect in RH (no subregion follow-up, no LOOCV).") |>
  flextable::body_add_flextable(fmt_resumen_ft(construir_tabla_resumen(res_fsoma_RH))) |>
  officer::body_add_par("")

doc <- agregar_modelos_word(doc, res_fsoma_RH, "Right HC SANDI fsoma — full model",
                            "Model: MaxLevel ~ Education + fsoma. OLS.")

print(doc, target = "results_regressions_HC_detailled.docx")
cat("Saved: results_regressions_HC_detailled.docx\n")
cat("\nFigures saved as TIFF (600 dpi, LZW):\n")
cat("- scatter_HC_completo.tiff\n- scatter_anatomico.tiff\n- scatter_subfields.tiff\n- scatter_eje.tiff\n")
cat("- scatter_SANDI_fneurite.tiff\n- scatter_SANDI_eje.tiff\n- scatter_SANDI_fsoma.tiff\n")
cat("- scatter_HC_completo_RH.tiff\n- scatter_SANDI_fneurite_RH.tiff\n- scatter_SANDI_fsoma_RH.tiff\n")
