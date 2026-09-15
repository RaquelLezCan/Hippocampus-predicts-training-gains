# Direct comparison of associations using dependent partial correlations, for the
# volume-and-SANDI manuscript (outcome: MaxLevel).
#
# Only the two confirmatory contrasts are reported: (1) left vs. right whole-HC
# volume, and (2) fneurite vs. fsoma in the left HC. The subdivisions (A-P axis,
# subfields) are exploratory and are not entered into this test.
#
# Method: partial correlation between MaxLevel and each measure (residualizing on
# Education, plus TBV for volume), Williams (1959) test for dependent correlations
# sharing a variable, effect size as the difference of partial correlations (delta r),
# 95% CI from 5,000 bootstrap resamples, and permutation p from 5,000 permutations.
# Partial correlations are compared instead of entering left+right in a single OLS
# because the left-right collinearity (r = .88) flips the partial betas.

library(readxl)

# 0. Data
ruta <- "Hipp_integrity_predicts_Ind_dif.xlsx"
raw  <- as.data.frame(read_excel(ruta))
d <- data.frame(
  MaxLevel   = raw[["Max_level"]],
  Education  = raw[["Education"]],
  TBV        = raw[["ses_pre_TBV"]],
  HC_L       = raw[["ses_pre_lh_whole_HC"]],
  HC_R       = raw[["ses_pre_rh_whole_HC"]],
  fneurite_L = raw[["fneurite_LH_HC_total"]],
  fsoma_L    = raw[["fsoma_LH_HC_total"]]
)
cat("n =", nrow(d), "| complete cases =", sum(complete.cases(d)), "\n")

# 1. Functions
resid_on <- function(v, covars, data)
  residuals(lm(as.formula(paste(v, "~", paste(covars, collapse = "+"))), data = data))

# Williams's t comparing r(Y,x1) and r(Y,x2), dependent and sharing Y
williams <- function(r12, r13, r23, n) {
  detR <- 1 - r12^2 - r13^2 - r23^2 + 2*r12*r13*r23
  df <- n - 3
  t  <- (r12 - r13) * sqrt( ((n-1)*(1+r23)) /
        ( 2*((n-1)/(n-3))*detR + ((r12+r13)^2/4)*(1-r23)^3 ) )
  c(t = t, df = df, p = 2*pt(-abs(t), df))
}

compare_pcor <- function(data, Y, x1, x2, covars, label,
                         nboot = 5000, nperm = 5000, seed = 1) {
  set.seed(seed)
  data <- data[complete.cases(data[, c(Y, x1, x2, covars)]), , drop = FALSE]
  n <- nrow(data); k <- length(covars); n_adj <- n - k   # n adjusted for covariates
  ry  <- resid_on(Y,  covars, data)
  r1v <- resid_on(x1, covars, data)
  r2v <- resid_on(x2, covars, data)
  r1  <- cor(ry, r1v); r2 <- cor(ry, r2v); r12 <- cor(r1v, r2v)
  w   <- williams(r1, r2, r12, n_adj)
  # 95% CI of the difference by bootstrap
  bo <- replicate(nboot, {
    i <- sample(n, n, TRUE); dd <- data[i, ]
    ry_ <- resid_on(Y, covars, dd); a <- resid_on(x1, covars, dd); b <- resid_on(x2, covars, dd)
    cor(ry_, a) - cor(ry_, b) })
  ci <- quantile(bo, c(.025, .975), na.rm = TRUE)
  # permutation p (permutes the residual of Y)
  obs <- r1 - r2
  pe  <- replicate(nperm, { rys <- sample(ry); cor(rys, r1v) - cor(rys, r2v) })
  pperm <- (sum(abs(pe) >= abs(obs)) + 1) / (nperm + 1)
  data.frame(comparison = label,
             r_partial_1 = round(r1, 3), r_partial_2 = round(r2, 3),
             diff_r = round(r1 - r2, 3), CI_low = round(ci[1], 3), CI_high = round(ci[2], 3),
             t = round(w["t"], 2), df = w["df"], p = round(w["p"], 4),
             p_perm = round(pperm, 4), r_between = round(r12, 3),
             row.names = NULL)
}

# 2. The two confirmatory contrasts
res <- rbind(
  compare_pcor(d, "MaxLevel", "HC_L", "HC_R", c("Education","TBV"),
               "Volume: Left vs Right"),
  compare_pcor(d, "MaxLevel", "fneurite_L", "fsoma_L", "Education",
               "Left HC: fneurite vs fsoma")
)
# Holm across the two pre-specified contrasts
res$p_holm <- round(p.adjust(res$p, "holm"), 4)

options(width = 200)
cat("\n==== Direct comparison of associations (partial correlations) ====\n\n")
print(res, row.names = FALSE)
cat("\ndiff_r = r_partial_1 - r_partial_2. p = Williams; p_perm = permutation (5,000).\n",
    "r_between = correlation between the two predictors (dependency).\n", sep = "")
