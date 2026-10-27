#' ============================================================
#' nhst-sampling-variability-lm-b.R   (VERSION B: one predictor)
#'
#' Shiny-free helpers behind VERSION B of the linear-model sampling-variability
#' app (nhst-sampling-variability-lm-b-app/app.R). Same idea as the correlation
#' app (nhst-sampling-variability.R), but the statistic is the slope b1 (and R^2)
#' of the simple linear model
#'
#'     y = b0 + b1 * x + error
#'
#' where x and y are two of the project's tests (reading experience = ART3,
#' vocabulary = StuVoc1), in RAW score units.
#'
#'   column 1: the observed ("focal") dataset: scatterplot + geom_smooth(lm)
#'   column 2: four more datasets from the same population (same slope, same N)
#'   column 3: two piles of dots, 100 studies each: H0 true (slope = 0, upper)
#'             and HA true (slope = the slider, lower); dots blue = not
#'             significant, orange = significant
#'
#' PHASE A: one continuous predictor (reading experience or vocabulary).
#' PHASE B: gender (M vs F) as a binary predictor, same file. The slope is then
#' the M - F mean difference (F = reference, coded 0; M = 1).
#'
#' This file is Version B only: it is the current Version C helper file
#' (nhst-sampling-variability-lm.R) with the Phase C features removed (age, adjusting
#' for other predictors, added-variable plots, added R^2, VIF). Function names are the
#' same as in Version C, so DO NOT source both files into one R session. With no
#' other predictors ticked, Version C gives identical results (checked by a parity
#' script; see nhst-sampling-variability-lm-b-walkthrough.qmd).
#'
#' Why this can reuse the correlation machinery. In simple regression the t-test
#' of b1 is the same test as the t-test of r (same t, df and p), and R^2 = r^2.
#' So the population is still simulated with simulate_pair() at a population
#' correlation r = b1 * sd(x) / sd(y); what is new is the estimation of b1 and
#' the conversion between slope and correlation.
#'
#' Built on, and modifying none of:
#'   psych-battery-dgp.R: simulate_pair(), get_calibration(), default_*()
#'   nhst-fisher-neyman-pearson-demo.R: simulate_many_studies(),
#'     plot_scatter_grid(), plot_dotpile_significance(), critical_r()
#'   nhst-sampling-variability.R: expected_power(), format_p(), r_size_label(),
#'     sv_n_choices, sv_r_limits, sv_point_color
#'
#' Seeds: focal = master seed, four-grid = master + 20000 (+ study number),
#' piles >= 1e6 (see simulate_lm_piles()), so none can coincide.
#'
#' Known approximations (documented in the walkthrough):
#'  * The "true slope" is NOMINAL: it is converted to a population r using the
#'    model-implied SDs, but the realised observed r runs a few per cent below
#'    nominal (linear disattenuation applied to a nonlinear item model).
#'  * The dashed "critical slope" lines use population SDs; each study's
#'    exact critical slope is t_crit * its own SE, which varies from study to
#'    study. Dot colours always use each study's exact p.
#'
#' Console usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("nhst-fisher-neyman-pearson-demo.R")
#'   source("nhst-sampling-variability.R")
#'   source("nhst-sampling-variability-lm-b.R")
#'   fg <- simulate_lm_focal_and_grid(n = 30, beta = .2, x = "ART3", y = "StuVoc1", seed = 1)
#'   pl <- simulate_lm_piles(n = 30, beta = .2, x = "ART3", y = "StuVoc1", pile_seed = 1)
#'   plot_lm_focal(fg)
#'   plot_lm_grid(fg)
#'   plot_lm_pile(pl, "h0", "beta", fg$focal_stats)
#'   plot_lm_pile(pl, "ha", "r2", fg$focal_stats)
#'   ```
#' ============================================================

library(ggplot2)

# ---- variables --------------------------------------------------------------------

#' The two continuous predictors/outcomes (the two tests), as
#' label = column name (the labels are what the app's drop-down shows).
lm_var_choices <- c("Reading experience" = "ART3", "Vocabulary" = "StuVoc1")

#' Predictor choices: the two tests and gender (binary: F = 0 reference, M = 1).
lm_x_choices <- c(lm_var_choices, "Gender (M vs F)" = "gender")

#' TRUE for the binary predictor.
lm_is_binary <- function(v) identical(v, "gender")

#' Plain-language label for a column name.
lm_var_label <- function(v) names(lm_x_choices)[match(v, lm_x_choices)]

#' Score range of each test (50 and 90 items), fixed so panels and redraws are
#' comparable.
lm_var_limits <- list(ART3 = c(0, 90), StuVoc1 = c(0, 50))

#' The other variable of the pair (there are exactly two tests).
lm_other_var <- function(v) setdiff(unname(lm_var_choices), v)

#' Population SD of a test's observed score, from the calibrated item model
#' (the same machinery simulate_pair() uses), so slope <-> r conversions need
#' no simulation.
lm_pop_sd <- function(v) {
  # Gender: two exactly equal groups coded 0/1 -> population SD 0.5. With this,
  # r = b * sd_x / sd_y gives r = (M - F difference) / (2 * sd_y), the equal-groups
  # point-biserial correlation, so every slope <-> r function below works unchanged.
  if (lm_is_binary(v)) return(0.5)
  sp  <- default_test_specs()[[v]]
  rel <- default_reliability_vec()[[v]]
  get_calibration(sp$mu, rel, sp$n_items, sp$guessing)$moments$sd_total
}

# ---- slope <-> correlation ------------------------------------------------------------

#' Nominal population correlation for a population slope of y on x:
#' r = b1 * sd(x) / sd(y).
lm_beta_to_r <- function(beta, x, y) beta * lm_pop_sd(x) / lm_pop_sd(y)

#' Inverse: b1 = r * sd(y) / sd(x), in units of y per unit of x.
lm_r_to_beta <- function(r, x, y) r * lm_pop_sd(y) / lm_pop_sd(x)

#' Slider step in raw units: 0.01 for the tests, 0.1 score points for the
#' (much larger) gender mean difference.
lm_slope_step <- function(x) if (lm_is_binary(x)) 0.1 else 0.01

#' Largest |slope| the slider offers: the slope that corresponds to the
#' population-r cap `sv_r_limits[2]` (.85; attenuation ceiling), rounded down
#' to the slider's step.
lm_slope_limit <- function(x, y) {
  st <- lm_slope_step(x)
  round(floor(lm_r_to_beta(sv_r_limits[2], x, y) / st + 1e-9) * st, 2)
}

#' The slope that corresponds to r = 1 (x-axis half-width of the full-range
#' slope pile).
lm_slope_scale <- function(x, y) lm_r_to_beta(1, x, y)

# ---- estimation --------------------------------------------------------------------------

#' Slope, SE, t, df, exact two-sided p and R^2 of the simple regression of y on x
#' (closed form; identical to summary(lm(y ~ x)), checked in the tests).
#'
#' @return one-row data frame: `beta`, `se` (raw units), `beta_std`, `se_std`
#'   (standardised slope = r, and its SE, from the sample SDs), `t`, `df`, `p`,
#'   `r2`, `significant`.
slope_stats <- function(x, y, alpha = 0.05) {
  n   <- length(x)
  sxx <- sum((x - mean(x))^2)
  sxy <- sum((x - mean(x)) * (y - mean(y)))
  syy <- sum((y - mean(y))^2)
  b   <- sxy / sxx
  r2  <- sxy^2 / (sxx * syy)
  se  <- sqrt(syy * (1 - r2) / (n - 2) / sxx)
  t   <- b / se
  p   <- 2 * stats::pt(-abs(t), n - 2)
  data.frame(beta = b, se = se, beta_std = sxy / sqrt(sxx * syy), se_std = se * sqrt(sxx / syy),
             t = t, df = n - 2, p = p, r2 = r2, significant = p < alpha)
}

#' Name used for the slope in labels: "b" (raw units) or "beta" (standardised;
#' spelled out, not a Greek letter, so it renders in any font, including webR).
lm_sym <- function(units) if (identical(units, "std")) "beta" else "b"

#' "b = 0.21 (SE 0.09), t(28) = 2.34, p = .027, R^2 = .16" (this order); with
#' `two_lines = TRUE` the line breaks after t(df) (for the narrow focal panel).
#' `units` = "raw" or "std".
lm_stat_label <- function(s, with_p = TRUE, two_lines = FALSE, units = "raw") {
  std <- identical(units, "std")
  r2  <- lm_r2_text(s)
  est <- sprintf("%s = %.2f (SE %.2f), t(%d) = %.2f", lm_sym(units), if (std) s$beta_std else s$beta,
                 if (std) s$se_std else s$se, as.integer(s$df), s$t)
  tail <- if (with_p) paste(format_p(s$p), r2, sep = ", ") else r2
  paste0(est, if (two_lines) ",\n" else ", ", tail)
}

#' R^2 text: "R^2 = .16" (leading zero dropped).
lm_r2_text <- function(s) sub("= 0\\.", "= .", sprintf("R\u00b2 = %.2f", s$r2))

#' Wording for a slope: direction and the Cohen-benchmark size of its
#' standardised version (= r).
lm_beta_label <- function(beta, x, y) r_size_label(lm_beta_to_r(beta, x, y))

# ---- simulation ----------------------------------------------------------------------------

#' Observed test scores from latent standardised abilities `z`: the same item
#' model simulate_pair() applies (logistic item curve, guessing floor, binomial
#' item noise), written out so a latent z built any way can be scored.
lm_item_scores <- function(z, var, reliability = default_reliability_vec()[[var]]) {
  n   <- length(z)
  sp  <- default_test_specs()[[var]]
  cal <- get_calibration(sp$mu, reliability, sp$n_items, sp$guessing)
  p   <- sp$guessing + (1 - sp$guessing) * plogis(cal$theta_mu + cal$theta_sigma * z)
  if (identical(var, "ART3")) {
    rbinom(n, sp$split[["authors"]], p) + rbinom(n, sp$split[["foils"]], p)
  } else {
    rbinom(n, sp$n_items, p)
  }
}

#' Simulate `n` people with gender (balanced groups) and one test score.
#'
#' Gender is coded F = 0 (reference), M = 1 and the groups are exactly equal
#' (n / 2 each; for odd n the extra person is F), so group sizes add no noise.
#' The outcome's latent ability is z = b * g + sqrt(1 - b^2) * e, with g = -1 / +1
#' (unit population variance) and e ~ N(0, 1); gender has no measurement error,
#' so the latent b is the target observed correlation divided by sqrt(reliability).
#' As with simulate_pair(), the TOTAL variance of the outcome is held fixed and
#' only the correlation changes.
#'
#' @param r  target observed-score (point-biserial) correlation between gender and `y`.
#' @return a data frame with columns `gender` (0 / 1) and `y`.
lm_simulate_gender <- function(n, r, y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  b <- disattenuate(r, 1, default_reliability_vec()[[y]])
  g <- c(rep(0, ceiling(n / 2)), rep(1, floor(n / 2)))
  z <- b * (g - 0.5) / 0.5 + sqrt(1 - b^2) * rnorm(n)
  out <- data.frame(gender = g)
  out[[y]] <- lm_item_scores(z, y)
  out
}

#' One simulated study for any (x, y): a test pair via simulate_pair(), or
#' gender + a test via lm_simulate_gender(). Returns the columns `x` and `y`.
lm_simulate <- function(n, r, x, y, seed = NULL) {
  if (lm_is_binary(x)) return(lm_simulate_gender(n, r, y, seed))
  simulate_pair(n, x, y, true_corr = r, seed = seed)[c(x, y)]
}

#' `n_studies` independent studies (study s uses seed + s, exactly as
#' simulate_many_studies() does, so test pairs are unchanged), as one long data
#' frame with a `study` factor.
lm_simulate_many <- function(n, n_studies, r, x, y, seed) {
  d <- do.call(rbind, lapply(seq_len(n_studies), function(s) {
    cbind(lm_simulate(n, r, x, y, seed + s), study = sprintf("Study %d", s))
  }))
  d$study <- factor(d$study, levels = sprintf("Study %d", seq_len(n_studies)))
  d
}

#' Simulate the focal dataset and the four comparison datasets.
#'
#' @param n,beta,seed  sample size, population slope (raw units; for gender, the M - F
#'   difference), seed.
#' @param x,y          column names (predictor, outcome).
#' @return a list: `n`, `beta`, `r_pop` (the standardised slope = r), `x`, `y`, `focal`
#'   (data frame), `focal_stats` (one row, see slope_stats()), `grid` (long data frame with
#'   a `study` factor), `grid_stats` (one row per study; the strip labels are built in
#'   plot_lm_grid(), so they can follow the units).
simulate_lm_focal_and_grid <- function(n, beta, x, y, seed, alpha = 0.05) {
  r_pop <- lm_beta_to_r(beta, x, y)
  focal <- lm_simulate(n, r_pop, x, y, seed)
  focal_stats <- slope_stats(focal[[x]], focal[[y]], alpha)

  grid <- lm_simulate_many(n, 4, r_pop, x, y, seed + 20000)
  grid_stats <- do.call(rbind, lapply(levels(grid$study), function(s) {
    g <- grid[grid$study == s, ]
    cbind(study = s, slope_stats(g[[x]], g[[y]], alpha))
  }))
  list(n = n, beta = beta, r_pop = r_pop, x = x, y = y, focal = focal,
       focal_stats = focal_stats, grid = grid, grid_stats = grid_stats)
}

#' Simulate the two 100-study piles: H0 true (population slope 0) and HA true (population
#' slope `beta`), both at sample size `n`.
#'
#' Seed design: study `i` of a world uses
#' `pile_seed + offset + i*1e6 + n*1e3`, with offset 0 (H0) or 12345 (HA). The
#' offset can never be cancelled by the other terms, so no seed is shared
#' between worlds even when beta = 0 (the two piles are then independent draws
#' from the same population). Always >= 1e6, so never equal to a focal/grid seed.
#'
#' @return a list: `n`, `beta`, `r_pop`, `x`, `y`, `alpha`, `h0`, `ha` (one row per study:
#'   slope_stats()), `crit_r`, `crit_beta`, `crit_r2`, `crit_by_world` (approximate critical
#'   values, the same in both worlds, for the slope in raw and standardised units and for
#'   R^2), `true_r2` (the HA world's population R^2), `r2_label`, `expected_h0`, `expected_ha`.
simulate_lm_piles <- function(n, beta, x, y, pile_seed, reps = 100, alpha = 0.05) {
  r_pop <- lm_beta_to_r(beta, x, y)
  one_world <- function(r, offset) {
    do.call(rbind, lapply(seq_len(reps), function(i) {
      d <- lm_simulate(n, r, x, y, pile_seed + offset + i * 1e6 + n * 1e3)
      slope_stats(d[[x]], d[[y]], alpha)
    }))
  }
  crit_r <- critical_r(n, alpha)
  cw <- list(beta = lm_r_to_beta(crit_r, x, y), beta_std = crit_r, r2 = crit_r^2)
  list(n = n, beta = beta, r_pop = r_pop, x = x, y = y, alpha = alpha,
       h0 = one_world(0, 0), ha = one_world(r_pop, 12345),
       crit_r = crit_r, crit_beta = cw$beta, crit_r2 = cw$r2,
       crit_by_world = list(h0 = cw, ha = cw), true_r2 = r_pop^2, r2_label = "R\u00b2",
       expected_h0 = alpha, expected_ha = expected_power(n, r_pop, alpha))
}


# ---- plots ---------------------------------------------------------------------------------

#' Colours (colour-blind safe: no red/blue pairing). Points keep the project's
#' blue; the fitted line, its band and the "our observed" marker are near-black,
#' which cannot be confused with the Okabe-Ito blue/orange significance dots.
lm_ink_color <- "#222222"

#' Colour of the fitted line, its confidence band and the group means / CIs:
#' Okabe-Ito bluish green, distinct from the blue points for all common forms of
#' colour-blindness. Each is drawn over a thicker white "halo" so it stays visible
#' against a dense cloud of points (large N).
lm_interval_color <- "#009E73"

#' Point transparency in the focal panel (a little lower than 1 so the fitted line
#' and intervals show through at large N).
lm_focal_alpha <- 0.7

#' The regression line + confidence band drawn on the scatterplots.
lm_smooth <- function(linewidth = 1) {
  list(
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "white",
                linewidth = linewidth + 1.8),
    geom_smooth(method = "lm", formula = y ~ x, color = lm_interval_color, fill = lm_interval_color,
                alpha = 0.15, linewidth = linewidth)
  )
}

#' A text label on a white background (no border) anchored in the headroom of a
#' pile, so dashed/dotted lines never run through the text. ggplot2 4.0 renamed the
#' border argument, so this checks the version (webR may ship an older ggplot2).
lm_label <- function(x, text, vjust, hjust, color) {
  border <- if (utils::packageVersion("ggplot2") >= "4.0.0") list(linewidth = 0) else list(label.size = 0)
  do.call(annotate, c(list("label", x = x, y = Inf, label = text, vjust = vjust, hjust = hjust,
                           size = 3.3, color = color, fill = "white", alpha = 0.9,
                           label.padding = grid::unit(0.12, "lines")), border))
}

#' Group means with 95% t confidence intervals, by `by` columns (e.g. "group", or
#' c("study", "group")).
lm_group_means <- function(d, y, by) {
  keys <- unique(d[by])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    sel <- Reduce(`&`, lapply(by, function(b) d[[b]] == keys[[b]][i]))
    v <- d[[y]][sel]
    m <- mean(v)
    h <- stats::qt(0.975, length(v) - 1) * stats::sd(v) / sqrt(length(v))
    cbind(keys[i, , drop = FALSE], mean = m, lo = m - h, hi = m + h)
  }))
}

#' Layers for a two-group "dot plot": jittered points, the two group means with
#' 95% CIs, and a line joining the means (its slope IS b, the M - F difference).
lm_group_layers <- function(means, y, point_alpha = lm_focal_alpha, point_size = 3) {
  list(
    geom_point(color = sv_point_color, alpha = point_alpha, size = point_size,
               position = position_jitter(width = 0.15, height = 0, seed = 1)),
    geom_errorbar(data = means, aes(y = mean, ymin = lo, ymax = hi), width = 0.16,
                  color = "white", alpha = 0.6, linewidth = 2.4),
    geom_line(data = means, aes(y = mean, group = 1), color = "white", alpha = 0.6, linewidth = 2.8),
    geom_errorbar(data = means, aes(y = mean, ymin = lo, ymax = hi), width = 0.12,
                  color = lm_interval_color, linewidth = 0.9),
    geom_line(data = means, aes(y = mean, group = 1), color = lm_interval_color, linewidth = 1.1),
    geom_point(data = means, aes(y = mean), color = "white", alpha = 0.6, size = 5, shape = 18),
    geom_point(data = means, aes(y = mean), color = lm_interval_color, size = 3.8, shape = 18),
    scale_y_continuous(limits = lm_var_limits[[y]])
  )
}

#' Gender version of column 1: jittered points by group, group means with 95%
#' CIs and the line joining them; subtitle as in plot_lm_focal().
plot_lm_focal_binary <- function(fg, units = "raw") {
  d <- fg$focal
  d$group <- factor(d[[fg$x]], levels = c(0, 1), labels = c("F", "M"))
  ggplot(d, aes(group, .data[[fg$y]])) +
    lm_group_layers(lm_group_means(d, fg$y, "group"), fg$y) +
    labs(x = "Gender (F = reference, M)", y = lm_var_label(fg$y),
         subtitle = lm_stat_label(fg$focal_stats, two_lines = TRUE, units = units)) +
    theme_bw(base_size = 14) +
    theme(plot.subtitle = element_text(face = "bold", size = 12))
}

#' Strip labels for the four-dataset grid, on three lines: study number; slope
#' (in `units`) and t; R^2.
lm_strip_labels <- function(gs, units = "raw") {
  est <- if (identical(units, "std")) gs$beta_std else gs$beta
  sprintf("%s\n%s = %.2f, t = %.2f\n%s", gs$study, lm_sym(units), est, gs$t,
          vapply(seq_len(nrow(gs)), function(i) lm_r2_text(gs[i, ]), character(1)))
}

#' Gender version of column 2: the four comparison datasets as a 2 x 2 grid.
plot_lm_grid_binary <- function(fg, units = "raw") {
  gs  <- fg$grid_stats
  labels <- lm_strip_labels(gs, units)
  d <- fg$grid
  d$study <- factor(d$study, levels = gs$study, labels = labels)
  d$group <- factor(d[[fg$x]], levels = c(0, 1), labels = c("F", "M"))
  ggplot(d, aes(group, .data[[fg$y]])) +
    lm_group_layers(lm_group_means(d, fg$y, c("study", "group")), fg$y,
                    point_alpha = 0.6, point_size = 1.8) +
    facet_wrap(~study, ncol = 2) +
    labs(x = "Gender (F = reference, M)", y = lm_var_label(fg$y)) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 8, hjust = 0))
}

#' Column 1: the focal scatterplot with the fitted line, annotated with the
#' slope (raw `b` or standardised `beta`), SE, t, p and R^2. The axes stay in raw
#' score units whichever `units` is chosen.
#'
#' @param fg     a simulate_lm_focal_and_grid() result.
#' @param units  "raw" or "std" (label only).
plot_lm_focal <- function(fg, units = "raw") {
  if (lm_is_binary(fg$x)) return(plot_lm_focal_binary(fg, units))
  ggplot(fg$focal, aes(.data[[fg$x]], .data[[fg$y]])) +
    geom_point(color = sv_point_color, alpha = lm_focal_alpha, size = 3) +
    lm_smooth() +
    scale_x_continuous(limits = lm_var_limits[[fg$x]]) +
    scale_y_continuous(limits = lm_var_limits[[fg$y]]) +
    labs(x = lm_var_label(fg$x), y = lm_var_label(fg$y),
         subtitle = lm_stat_label(fg$focal_stats, two_lines = TRUE, units = units)) +
    theme_bw(base_size = 14) +
    theme(plot.subtitle = element_text(face = "bold", size = 12))
}

#' Column 2: the four comparison datasets as a 2 x 2 grid, each with its own
#' fitted line, each strip showing that study's slope (in `units`) and R^2.
#'
#' @param fg     a simulate_lm_focal_and_grid() result.
#' @param units  "raw" or "std" (strip labels only).
plot_lm_grid <- function(fg, units = "raw") {
  if (lm_is_binary(fg$x)) return(plot_lm_grid_binary(fg, units))
  gs  <- fg$grid_stats
  labels <- lm_strip_labels(gs, units)
  d <- fg$grid
  d$study <- factor(d$study, levels = gs$study, labels = labels)
  plot_scatter_grid(d, fg$x, fg$y, ncol = 2) +
    lm_smooth(linewidth = 0.8) +
    scale_x_continuous(limits = lm_var_limits[[fg$x]]) +
    scale_y_continuous(limits = lm_var_limits[[fg$y]]) +
    labs(x = lm_var_label(fg$x), y = lm_var_label(fg$y)) +
    theme(strip.text = element_text(face = "bold", size = 8, hjust = 0))
}

#' x-axis range shared by both piles. Full scale: slope pile +/- the slope that
#' corresponds to r = 1 (raw: sd_y/sd_x; standardised: 1); R^2 pile 0 to 1. With
#' `zoom = TRUE`: the range of everything drawn, padded by 5% (R^2 floored at 0,
#' slope within the full range).
lm_pile_xlim <- function(piles, stat, observed, zoom = FALSE, units = "raw") {
  std <- identical(units, "std")
  if (stat == "beta") {
    scale <- if (std) 1 else lm_slope_scale(piles$x, piles$y)
    col   <- if (std) "beta_std" else "beta"
    crit  <- unlist(lapply(piles$crit_by_world, function(w) if (std) w$beta_std else w$beta))
    if (!zoom) return(c(-scale, scale))
    rng <- range(c(piles$h0[[col]], piles$ha[[col]], observed[[col]], crit, -crit))
    pad <- 0.05 * diff(rng)
    c(max(-scale, rng[1] - pad), min(scale, rng[2] + pad))
  } else {
    if (!zoom) return(c(0, 1))
    rng <- range(c(piles$h0$r2, piles$ha$r2, observed$r2, unlist(lapply(piles$crit_by_world, `[[`, "r2"))))
    pad <- 0.05 * diff(rng)
    c(max(0, rng[1] - pad), min(1, rng[2] + pad))
  }
}

#' Fraction of the panel height the tallest stack may fill (calibrated so that no
#' overflow warning appears across N, zoom, statistic and panel sizes; see the
#' walkthrough).
lm_dot_fill <- 0.5

#' Multiplier on the drawn dot diameter (ggdist `dotsize`; 1 = dots exactly one bin wide, so
#' neighbours just touch). Above 1 the dots overlap slightly: bigger and easier to see.
lm_dot_size <- 1.25

#' Assumed size (pixels, width x height) of one pile's plotting panel in the app; sets the
#' shared dot size. Matches the layout of nhst-sampling-variability-lm-b-app/app.R.
lm_pile_panel_px <- c(430, 150)

#' A dot diameter (in data units) shared by BOTH piles. ggdist picks each pile's
#' binwidth separately so that its own tallest stack fits, which makes the more
#' concentrated pile's dots smaller (the H0 pile, usually); with one shared
#' binwidth the dots mean the same thing in both. Chooses the largest binwidth
#' for which the tallest stack in either pile stays below `fill` of the panel
#' height, given the panel's shape (`panel_px` = width and height in pixels).
#'
#' @param values_list  list of the numeric vectors that will be plotted.
#' @param xlim         the x range of the panel.
#' @param panel_px     c(width, height) of the plotting panel in pixels.
lm_shared_binwidth <- function(values_list, xlim, panel_px = lm_pile_panel_px, fill = lm_dot_fill) {
  span   <- diff(xlim)
  aspect <- panel_px[1] / panel_px[2]
  bws    <- exp(seq(log(span / 3000), log(span / 12), length.out = 300))
  # Most dots within +/- 0.75 binwidth of any one value: a slightly pessimistic
  # stand-in for ggdist's stack height (it never underestimates it in tests).
  tallest <- function(bw) {
    max(vapply(values_list, function(v) {
      v <- sort(v)
      max(findInterval(v + 0.75 * bw, v) - findInterval(v - 0.75 * bw, v, left.open = TRUE))
    }, numeric(1)))
  }
  ok <- vapply(bws, function(bw) tallest(bw) * bw * aspect / span <= fill, logical(1))
  if (any(ok)) max(bws[ok]) else min(bws)
}

#' Column 3: one pile of 100 dots (one per simulated study), blue = not
#' significant, orange = significant, with the world's population value
#' (dotted), the critical value(s) (dashed) and the focal dataset's observed
#' value (solid black). Built on plot_dotpile_significance(); the two text
#' labels are drawn here, on white backgrounds and on separate rows, so lines
#' never run through them.
#'
#' @param piles     a simulate_lm_piles() result.
#' @param world     "h0" or "ha".
#' @param stat      "beta" (the slope) or "r2".
#' @param observed  the focal dataset's slope_stats() row.
#' @param zoom      see lm_pile_xlim().
#' @param units     "raw" or "std": what the slope pile shows (ignored for R^2).
#' @param panel_px  c(width, height) of the plotting panel in pixels (sets the shared dot size).
#' @param dotsize   multiplier on the drawn dot diameter (see lm_dot_size).
plot_lm_pile <- function(piles, world = c("h0", "ha"), stat = c("beta", "r2"), observed, zoom = FALSE,
                         units = c("raw", "std"), panel_px = lm_pile_panel_px, dotsize = lm_dot_size) {
  world <- match.arg(world)
  stat  <- match.arg(stat)
  units <- match.arg(units)
  std   <- units == "std"
  bin   <- lm_is_binary(piles$x)
  sweep <- piles[[world]]
  xlim  <- lm_pile_xlim(piles, stat, observed, zoom, units)
  expected <- if (world == "h0") piles$expected_h0 else piles$expected_ha
  k <- sum(sweep$significant)

  if (stat == "beta") {
    sym      <- lm_sym(units)
    col      <- if (std) "beta_std" else "beta"
    values   <- sweep[[col]]
    pop      <- if (std) piles$r_pop else piles$beta
    true_val <- if (world == "h0") 0 else pop
    line_nm  <- if (world == "h0") sprintf("%s = 0 (H0)", sym) else sprintf("true %s = %.2f", sym, true_val)
    cw       <- piles$crit_by_world[[world]]
    crit_abs <- if (std) cw$beta_std else cw$beta
    crit     <- c(-crit_abs, crit_abs)
    obs_val  <- observed[[col]]
    xlab     <- if (std) "Standardised slope (beta) in each simulated study" else
      if (bin) "Estimated M - F difference b per study" else
        "Estimated slope b in each simulated study"
    obs_lab  <- paste("our observed", sym)
    what     <- if (std) "population standardised beta" else
      if (bin) "M - F difference" else "population slope"
    title <- if (world == "h0") sprintf("If H0 is true (%s = 0)", what) else
      sprintf("If HA is true (%s = %.2f)", what, true_val)
  } else {
    values   <- sweep$r2
    rl       <- piles$r2_label
    true_val <- if (world == "h0") 0 else piles$true_r2
    line_nm  <- if (world == "h0") sprintf("%s = 0 (H0)", rl) else sprintf("true %s = %.2f", rl, true_val)
    crit     <- piles$crit_by_world[[world]]$r2
    obs_val  <- observed$r2
    xlab     <- sprintf("%s in each simulated study", rl)
    obs_lab  <- paste("our observed", rl)
    title <- if (world == "h0") sprintf("If H0 is true (population %s = 0)", rl) else
      sprintf("If HA is true (population %s = %.2f)", rl, true_val)
  }
  subtitle <- sprintf("%d of %d studies significant (about %.0f expected)", k, nrow(sweep),
                      round(expected * nrow(sweep)))

  side <- function(v) if (v > mean(xlim)) 1.02 else -0.02
  col_name <- if (stat == "r2") "r2" else if (std) "beta_std" else "beta"
  bw <- lm_shared_binwidth(list(piles$h0[[col_name]], piles$ha[[col_name]]), xlim, panel_px)
  p <- plot_dotpile_significance(values, xlab = xlab, significant = sweep$significant, binwidth = bw) +
    geom_vline(xintercept = true_val, color = "gray40", linetype = "dotted") +
    geom_vline(xintercept = crit, color = "gray40", linetype = "dashed") +
    geom_vline(xintercept = obs_val, color = lm_ink_color, linewidth = 1) +
    lm_label(true_val, line_nm, vjust = 1.1, hjust = side(true_val), color = "gray30") +
    lm_label(obs_val, obs_lab, vjust = 2.7, hjust = side(obs_val), color = lm_ink_color) +
    scale_y_continuous(limits = c(0, 1.4), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim) +
    labs(title = title, subtitle = subtitle) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
  # The dot layer is the shared function's first layer; scale its dot diameter
  p$layers[[1]]$geom_params$dotsize <- dotsize
  p
}
