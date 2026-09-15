# Shared visual style for the REFUEL project.
# Loaded at the top of each script with source("style_REFUEL.R") so that figures
# across analyses keep consistent colors and theme.

# Palette for the regression models (whole HC, subregions, A-P axis)
REFUEL_COL_REG_POINT <- "#2C5F8A"   # points
REFUEL_COL_REG_LINE  <- "#C84B31"   # fit line and CI band

# Palette for the cross-validation figures (kept distinct from the one above)
REFUEL_COL_CV_POINT <- "#6A4C93"
REFUEL_COL_CV_LINE  <- "#3B8F5C"
REFUEL_COL_CV_BAND  <- "#3B8F5C"
REFUEL_COL_CV_HIST  <- "#6A4C93"
REFUEL_COL_CV_OBS   <- "#3B8F5C"

# Point and fit-line geometry
REFUEL_POINT_SIZE  <- 3
REFUEL_POINT_ALPHA <- 0.75
REFUEL_SMOOTH_LWD  <- 1.1

# Base theme. base_size 14 for single plots; lower it (e.g. 12) for row panels
# of several subplots so the axes do not crowd.
theme_refuel <- function(base_size = 14, axis_title_size = 13, axis_text_size = 11,
                         show_title = TRUE) {
  tema <- theme_classic(base_size = base_size) +
    theme(
      axis.title    = element_text(size = axis_title_size),
      axis.text     = element_text(size = axis_text_size),
      plot.subtitle = element_text(size = base_size - 3, color = "gray30")
    )
  if (show_title) {
    tema <- tema + theme(plot.title = element_text(face = "bold", size = base_size - 1))
  } else {
    tema <- tema + theme(plot.title = element_blank())
  }
  tema
}

# Statistic annotation (beta or R2_CV, plus p) in the top-right corner. Symbols
# in italics, numbers in regular font. APA format: no leading zero (.44 not 0.44)
# and a typographic minus sign.
format_apa_num <- function(x, digits = 2) {
  s <- sprintf(paste0("%.", digits, "f"), x)
  s <- sub("^-", "−", s)
  sub("(−?)0\\.", "\\1.", s)
}

format_p_apa <- function(p) {
  if (p < .001) return("< .001")
  digits <- if (p < .01) 3 else 2
  paste0("= ", format_apa_num(p, digits = digits))
}

# label_expr: plotmath fragment for the left-hand symbol (e.g. "italic(beta)").
stat_annotation <- function(label_expr, value, p_value, value_digits = 2) {
  value_str <- format_apa_num(value, digits = value_digits)
  p_str     <- format_p_apa(p_value)
  expr_texto <- sprintf('paste(%s, " = %s,  ", italic(p), " %s")',
                        label_expr, value_str, p_str)
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
           label = expr_texto, parse = TRUE, size = 3.5, color = "gray20")
}

# p formatting for tables and console (not for the figures)
format_p <- function(p) {
  ifelse(p < .001, "< .001", paste0("= ", sprintf("%.3f", p)))
}

# Export figures as high-resolution TIFF
save_tiff <- function(plot, filename, width = 5, height = 4.5, dpi = 600) {
  ggsave(
    filename    = filename,
    plot        = plot,
    device      = "tiff",
    compression = "lzw",
    width       = width,
    height      = height,
    dpi         = dpi,
    bg          = "white"
  )
  cat("Saved (TIFF,", dpi, "dpi):", filename, "\n")
}

# Add a symmetric margin to a range, to set common axis limits across subplots
# that share the y-axis variable.
pad_range <- function(r, factor = 0.08) {
  span <- diff(r)
  r + c(-1, 1) * span * factor
}

# Format "fneurite" in figure text (f in italics, "neurite" as subscript). If the
# text does not contain "fneurite" it is returned unchanged, or bold if bold = TRUE.
# bold = TRUE is used in panel titles built with grid::textGrob(), where
# gpar(fontface) does not affect plotmath expressions.
format_fneurite_title <- function(texto, bold = FALSE) {

  if (!grepl("fneurite", texto, fixed = TRUE)) {
    if (bold) {
      return(parse(text = sprintf('bold("%s")', texto))[[1]])
    }
    return(texto)
  }

  partes  <- strsplit(texto, "fneurite", fixed = TRUE)[[1]]
  antes   <- partes[1]
  despues <- if (length(partes) > 1) partes[2] else ""

  if (bold) {
    f_expr       <- "bolditalic(f)"
    neurite_expr <- "bold(neurite)"
    antes_expr   <- sprintf('bold("%s")', antes)
    despues_expr <- sprintf('bold("%s")', despues)
  } else {
    f_expr       <- "italic(f)"
    neurite_expr <- "neurite"
    antes_expr   <- sprintf('"%s"', antes)
    despues_expr <- sprintf('"%s"', despues)
  }

  expr_texto <- sprintf('paste(%s, %s[%s], %s)', antes_expr, f_expr, neurite_expr, despues_expr)
  parse(text = expr_texto)[[1]]
}

# Excel sheet names: at most 31 characters, no odd characters, unique within a
# given set (adds a numeric suffix if truncation to 31 chars collides).
sheet_slug <- function(x, prefix = "", existing = character(0)) {
  s <- paste0(prefix, gsub("[^A-Za-z0-9]+", "_", x))
  s <- substr(s, 1, 31)

  if (s %in% existing) {
    base <- substr(s, 1, 28)
    i <- 2
    repeat {
      candidato <- paste0(base, "_", i)
      if (!(candidato %in% existing)) {
        s <- candidato
        break
      }
      i <- i + 1
    }
  }
  s
}
