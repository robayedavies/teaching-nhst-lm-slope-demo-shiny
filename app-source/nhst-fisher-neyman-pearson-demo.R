#' ============================================================
#' nhst-fisher-neyman-pearson-demo.R
#'
#' Demonstration figures contrasting Fisher's significance-testing
#' framework with Neyman & Pearson's hypothesis-testing framework, for the
#' hypothesis-testing-significance.qmd lecture. Built entirely on top of
#' psych-battery-dgp.R's `simulate_pair()`/`run_sweep()`, using one focal
#' pair throughout:
#'
#'   StuVoc1 (vocabulary, 50 items)  --  ART3 (print exposure, 90 items)
#'   Vermeiren-default observed r = .52 (default_r_obs_matrix())
#'
#' The lecture simulates this pair under two worlds:
#'   H0 true: true_corr = 0   (vocabulary and print exposure are unrelated)
#'   H1 true: true_corr = .5  (close to the Vermeiren-default .52)
#'
#' A single sample size, N = 30, is used for every *fixed-N* demonstration
#' (both the Fisher "suppose H0 is true" demos and the Neyman-Pearson
#' H0-true/H1-true demos) -- chosen because it is close to the sample size
#' needed for power = .80 at r = .5, alpha = .05 two-sided
#' (pwr::pwr.r.test(r = .5, sig.level = .05, power = .8) gives n \u2248 28.2;
#' the empirical power of this script's own DGP at n = 30 is \u2248 .82, see
#' the walkthrough). Using the *same* N for both worlds is deliberate: it
#' keeps "what changed?" unambiguous (only the true correlation), and it's
#' also statistically apt, since the empirical false-positive rate under
#' H0 doesn't depend on N -- alpha is a property of the decision rule, not
#' of the sample. Only the "funnel"/corridor demonstrations, which are
#' specifically about the effect of varying N, use N = 50/100/1000.
#'
#' Significance is encoded by colour alone throughout (no point shape),
#' using the Okabe-Ito colour-blind-safe pair already used by
#' person-grid.R: "#0072B2" (blue) = not significant, "#E69F00" (orange)
#' = significant.
#'
#' Requires: psych-battery-dgp.R sourced first (uses simulate_pair(),
#' fit_focal_cor(), run_sweep(), default_r_obs_matrix()); ggplot2; ggdist.
#'
#' Quarto usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("nhst-fisher-neyman-pearson-demo.R")
#'
#'   # Fisher section: one grid of 6 studies, H0 true, n = 30
#'   studies <- simulate_many_studies(n = 30, n_studies = 6, true_corr = 0, seed = 11)
#'   plot_scatter_grid(studies, "ART3", "StuVoc1")
#'
#'   # ... and the same 6 studies' r / t values as dot-piles
#'   fits <- fit_many_studies(studies, "ART3", "StuVoc1")
#'   plot_dotpile_stat(fits$estimate, xlab = "r")
#'   plot_dotpile_stat(t_from_r(fits$estimate, n = 30), xlab = "t",
#'                      vline = critical_t(30), vline_label = "critical t (alpha = .05)")
#'
#'   # Neyman-Pearson section: 1000 studies, H1 true, n = 30
#'   sweep <- run_sweep(n_grid = 30, true_corr_grid = .5, reps = 1000,
#'                       var_x = "StuVoc1", var_y = "ART3", seed_base = 22)
#'   plot_dotpile_significance(sweep$estimate, xlab = "r", significant = sweep$significant,
#'                              vlines = c(`critical r` = critical_r(30), `true r` = .5))
#'
#'   # Funnel plots: same effect, varying N
#'   plot_funnel("StuVoc1", "ART3", n_grid = c(50, 100, 1000), true_corr = .5, reps = 300, seed_base = 33)
#'   plot_funnel_multi_effect("StuVoc1", "ART3", n_grid = c(50, 100, 1000),
#'                             true_corr_grid = c(.25, .5, .75), reps = 300, seed_base = 44)
#'   ```
#' ============================================================

library(ggplot2)

# ---- significance palette (Okabe-Ito, colour-blind-safe; colour-only, no shape) ----

sig_colors <- c("not significant" = "#0072B2", "significant" = "#E69F00")

#' Turn a logical "significant" vector into the labelled factor used by
#' every colour-by-significance plot in this script, so the legend
#' wording and factor level order are identical everywhere.
significance_label <- function(significant) {
  factor(significant, levels = c(FALSE, TRUE), labels = c("not significant", "significant"))
}

# ---- r / t / critical-value helpers ----------------------------------------

#' The t-statistic for a Pearson correlation, given r and sample size n.
#' t = r * sqrt(n - 2) / sqrt(1 - r^2); see fit_focal_cor()/cor.test() for
#' the version computed directly from data -- this is the same formula,
#' applied to an already-computed r (e.g. from run_sweep()'s output).
t_from_r <- function(r, n) {
  r * sqrt(n - 2) / sqrt(1 - r^2)
}

#' The critical two-sided t-value for a Pearson correlation test with
#' n - 2 degrees of freedom.
critical_t <- function(n, alpha = 0.05) {
  stopifnot("`n` must be at least 3" = n >= 3)
  qt(1 - alpha / 2, df = n - 2)
}

#' The critical two-sided |r| boundary at sample size n -- the observed r
#' at which the t-test for that correlation exactly reaches the critical
#' t. Inverts t_from_r() algebraically: r_crit = t_crit / sqrt(t_crit^2 + n - 2).
critical_r <- function(n, alpha = 0.05) {
  tcrit <- critical_t(n, alpha)
  tcrit / sqrt(tcrit^2 + n - 2)
}

# ---- simulating small grids of studies (for scatter grids / small dot-piles) ----

#' Simulate `n_studies` independent studies of `n` participants each, all
#' drawn from the same (var_x, var_y) population correlation -- the small-
#' n_studies analogue of run_sweep(), returning the *raw* per-person data
#' (not just fitted estimates), for scatter-grid figures.
#'
#' @param n            participants per study.
#' @param n_studies    number of independent studies.
#' @param true_corr    population *observed*-score correlation to simulate
#'                      from (see simulate_pair()'s corr_scale argument).
#' @param var_x, var_y the focal pair of tests.
#' @param reliability  reliability override applied to both tests; NULL
#'                      (default) uses each test's Vermeiren-default value.
#' @param corr_scale   "observed" (default) or "latent"; passed through to
#'                      simulate_pair().
#' @param seed         base seed; study s uses seed + s, so each study is
#'                      independently drawn but the whole grid is
#'                      reproducible.
#' @return a long data frame: `study` (factor, "Study 1" .. "Study n_studies"),
#'   plus the `var_x`/`var_y` score columns from simulate_pair().
simulate_many_studies <- function(n, n_studies, true_corr = 0, var_x = "StuVoc1", var_y = "ART3",
                                   reliability = NULL, corr_scale = c("observed", "latent"), seed = NULL) {
  corr_scale <- match.arg(corr_scale)
  panels <- vector("list", n_studies)
  for (s in seq_len(n_studies)) {
    seed_s <- if (is.null(seed)) NULL else seed + s
    d <- simulate_pair(n, var_x, var_y, true_corr = true_corr, corr_scale = corr_scale,
                        reliability_x = reliability, reliability_y = reliability, seed = seed_s)
    d$study <- sprintf("Study %d", s)
    panels[[s]] <- d
  }
  out <- do.call(rbind, panels)
  out$study <- factor(out$study, levels = sprintf("Study %d", seq_len(n_studies)))
  out
}

#' Fit a Pearson correlation test within each `study` of a
#' simulate_many_studies()-style data frame, returning one row per study.
#'
#' @param data           a data frame with a `study` column (as returned
#'                         by simulate_many_studies()).
#' @param var_x, var_y   the focal pair of tests.
#' @param alpha          significance threshold.
#' @return a data frame: `study`, `estimate`, `ci_low`, `ci_high`,
#'   `p_value`, `significant`.
fit_many_studies <- function(data, var_x, var_y, alpha = 0.05) {
  studies <- levels(droplevels(data$study))
  rows <- lapply(studies, function(s) {
    fit <- fit_focal_cor(data[data$study == s, , drop = FALSE], var_x, var_y)
    cbind(study = s, fit, significant = fit$p_value < alpha)
  })
  out <- do.call(rbind, rows)
  out$study <- factor(out$study, levels = studies)
  out
}

# ---- scatter grid -----------------------------------------------------------

#' A facet_wrap(~study) grid of scatterplots, one panel per simulated
#' study -- shared theme/colour with the rest of the project's raw-data
#' scatterplots (plain blue points, no significance encoding: this is
#' "what the data look like", before any test is applied).
#'
#' @param data            a data frame with a `study` column (as returned
#'                          by simulate_many_studies()).
#' @param var_x, var_y    the focal pair of tests (column names in `data`).
#' @param ncol            panels per row.
#' @return a ggplot object.
plot_scatter_grid <- function(data, var_x, var_y, ncol = 3) {
  ggplot(data, aes(.data[[var_x]], .data[[var_y]])) +
    geom_point(color = "#3B6EA5", alpha = 0.6) +
    facet_wrap(~study, ncol = ncol) +
    labs(x = var_x, y = var_y) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
}

# ---- dot-piles: a statistic's distribution, with and without significance ----

#' A sensible default dot diameter (in data units) for ggdist::stat_dots(),
#' used whenever a dot-pile function's `binwidth` argument is left NULL.
#'
#' ggdist's own auto-binwidth (`binwidth = NA`) is tuned for the "many
#' points, one panel" case and works well there (verified up to n = 1000
#' for this script's r/t-value ranges), but for a *handful* of points
#' (e.g. the six-study demonstrations) it tends to pick a binwidth coarse
#' enough that distinct-but-nearby values get treated as ties and nudged
#' sideways into overlapping "petal" clusters rather than cleanly stacked
#' points. Below a small-n threshold, this instead scales the binwidth
#' directly off the data's own range and n, which keeps genuinely distinct
#' values visually separate while still stacking true (or near-)ties.
#' Above the threshold, ggdist's own auto-binwidth (`NA`) is used
#' unchanged.
#'
#' @param values  the numeric vector that will be plotted.
#' @return a single numeric binwidth, or `NA` (ggdist's own "auto").
default_dotpile_binwidth <- function(values) {
  n <- length(values)
  if (n <= 20) diff(range(values)) / (3 * n) else NA
}

#' Manually stack a small number of values into dot-pile-style (x, y)
#' positions, grouping each value with its neighbour into the same column
#' whenever they're within `binwidth` of each other (greedy, after
#' sorting). This exists only so that individual points can be labelled
#' exactly: `ggdist::stat_dots()`'s own stacking is computed at draw time
#' in grid units and isn't recoverable from `ggplot_build()`, so exact
#' per-point label placement (`plot_dotpile_labeled()`) isn't otherwise
#' possible. Intended for the same small-n regime as
#' `default_dotpile_binwidth()`'s n <= 20 branch -- not a general
#' replacement for ggdist's own layout.
#'
#' @param values     numeric vector to position.
#' @param binwidth   grouping tolerance; NULL (default) picks
#'                     `default_dotpile_binwidth()`.
#' @return a data frame: `x` (sorted values), `y` (stack level within its
#'   column, starting at 0), `order` (the index into the original
#'   `values` that each row corresponds to, for attaching per-point
#'   labels/colours after sorting).
stack_dot_positions <- function(values, binwidth = NULL) {
  if (is.null(binwidth) || is.na(binwidth)) binwidth <- default_dotpile_binwidth(values)
  if (is.na(binwidth)) binwidth <- diff(range(values)) / 30
  ord <- order(values)
  x <- values[ord]
  y <- integer(length(x))
  level <- 0
  for (i in seq_along(x)) {
    if (i == 1 || (x[i] - x[i - 1]) > binwidth) level <- 0 else level <- level + 1
    y[i] <- level
  }
  data.frame(x = x, y = y, order = ord)
}

#' A dot-pile with each point labelled (e.g. by study number), built on
#' `stack_dot_positions()` rather than `ggdist::stat_dots()` so that
#' labels land exactly on their point -- intended for small `n` (e.g. the
#' six-study demonstrations), where identifying which point came from
#' which study is pedagogically useful.
#'
#' @param values        numeric vector (r or t values, one per study).
#' @param labels         vector the same length as `values` (e.g. study
#'                         numbers) to print inside each dot.
#' @param xlab          x-axis label, e.g. "r" or "t".
#' @param vline         optional single reference value to draw as a
#'                        dotted vertical line.
#' @param vline_label    optional short label placed near the top of
#'                        `vline`.
#' @param binwidth       grouping tolerance; NULL (default) picks
#'                         `default_dotpile_binwidth()`.
#' @param dot_size, text_size  point size and label text size.
#' @return a ggplot object.
plot_dotpile_labeled <- function(values, labels, xlab, vline = NULL, vline_label = NULL,
                                  binwidth = NULL, dot_size = 11, text_size = 3) {
  pos <- stack_dot_positions(values, binwidth)
  pos$label <- labels[pos$order]
  p <- ggplot(pos, aes(x = x, y = y)) +
    geom_point(size = dot_size, color = "#0072B2") +
    geom_text(aes(label = label), color = "white", size = text_size, fontface = "bold") +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(limits = c(-0.6, max(pos$y) + 0.6)) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
  if (!is.null(vline)) {
    full_range <- range(c(values, vline))
    hjust <- if (vline > mean(full_range)) 1.05 else -0.05
    p <- p + geom_vline(xintercept = vline, color = "gray40", linetype = "dotted")
    if (!is.null(vline_label)) {
      p <- p + annotate("text", x = vline, y = Inf, label = vline_label,
                         vjust = 1.5, hjust = hjust, size = 3, color = "gray40")
    }
  }
  p
}

#' A single-colour "pile of dots" (ggdist::stat_dots) showing the
#' distribution of a statistic (r or t) across simulated studies, with no
#' significance encoding -- used for Fisher-style demonstrations, which
#' never classify an individual result as "significant"/"not significant"
#' (that distinction belongs to Neyman-Pearson).
#'
#' @param values        numeric vector (r or t values, one per study).
#' @param xlab          x-axis label, e.g. "r" or "t".
#' @param vline         optional single reference value (e.g. a critical
#'                        t, or one highlighted "observed" value) to draw
#'                        as a dotted vertical line.
#' @param vline_label    optional short label placed near the top of
#'                        `vline`, annotating what it represents.
#' @param binwidth       dot diameter in data units; NULL (default) picks
#'                        `default_dotpile_binwidth()`.
#' @return a ggplot object.
plot_dotpile_stat <- function(values, xlab, vline = NULL, vline_label = NULL, binwidth = NULL) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plot_dotpile_stat(). Install with install.packages('ggdist').")
  }
  if (is.null(binwidth)) binwidth <- default_dotpile_binwidth(values)
  df <- data.frame(value = values)
  p <- ggplot(df, aes(x = value)) +
    ggdist::stat_dots(color = "#0072B2", fill = "#0072B2", binwidth = binwidth) +
    scale_x_continuous(name = xlab) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, color = "gray40", linetype = "dotted")
    if (!is.null(vline_label)) {
      # Right-align (and tuck inside the panel) labels whose line sits in
      # the right half of the plotted range, so the label text can't run
      # off the right edge; left-align (as before) otherwise.
      full_range <- range(c(values, vline))
      hjust <- if (vline > mean(full_range)) 1.05 else -0.05
      p <- p + annotate("text", x = vline, y = Inf, label = vline_label,
                         vjust = 1.5, hjust = hjust, size = 3, color = "gray40")
    }
  }
  p
}

#' As plot_dotpile_stat(), but coloured by significance (Okabe-Ito pair,
#' colour only -- no shape), with one or more labelled reference lines
#' (e.g. the critical r boundary and/or the true population r) -- used
#' for every Neyman-Pearson-section dot-pile.
#'
#' @param values        numeric vector (r values, one per study).
#' @param xlab          x-axis label.
#' @param significant   logical vector, same length as `values`.
#' @param vlines        optional named numeric vector of reference values,
#'                        e.g. c(`critical r` = .31, `true r` = .5); names
#'                        are used as labels.
#' @param binwidth       dot diameter in data units; NULL (default) picks
#'                        `default_dotpile_binwidth()`.
#' @return a ggplot object.
plot_dotpile_significance <- function(values, xlab, significant, vlines = NULL, binwidth = NULL) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for plot_dotpile_significance(). Install with install.packages('ggdist').")
  }
  if (is.null(binwidth)) binwidth <- default_dotpile_binwidth(values)
  df <- data.frame(value = values, significant_label = significance_label(significant))
  p <- ggplot(df, aes(x = value, color = significant_label, fill = significant_label)) +
    ggdist::stat_dots(binwidth = binwidth) +
    scale_color_manual(values = sig_colors, name = NULL) +
    scale_fill_manual(values = sig_colors, name = NULL) +
    scale_x_continuous(name = xlab) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.9),
          axis.ticks.x = element_line(color = "black"))
  if (!is.null(vlines)) {
    # Order by position, so staggered vjust steps always run low-x-to-
    # high-x; right-align labels on the right half of the range (so their
    # text tucks back into the panel instead of running off the edge),
    # and stagger each label's vjust so two reference lines close enough
    # in x that their labels would otherwise collide (e.g. critical r and
    # true r sitting a few hundredths of r apart) are stacked instead.
    ord <- order(vlines)
    vline_df <- data.frame(x = unname(vlines)[ord], label = names(vlines)[ord])
    full_range <- range(c(values, vlines))
    vline_df$hjust <- ifelse(vline_df$x > mean(full_range), 1.05, -0.05)
    vline_df$vjust <- 1.5 + 1.5 * (seq_len(nrow(vline_df)) - 1)
    p <- p + geom_vline(data = vline_df, aes(xintercept = x), color = "gray40", linetype = "dotted",
                         inherit.aes = FALSE)
    for (i in seq_len(nrow(vline_df))) {
      p <- p + annotate("text", x = vline_df$x[i], y = Inf, label = vline_df$label[i],
                         vjust = vline_df$vjust[i], hjust = vline_df$hjust[i],
                         size = 3, color = "gray40")
    }
  }
  p
}

# ---- funnel / corridor plots: one effect, many n -----------------------------

#' A "funnel"/corridor plot for one fixed population correlation, one
#' panel per sample size -- each panel shows run_sweep()'s replicate
#' studies at that n as jittered points (x = observed r), coloured by
#' significance, with that panel's own critical-r boundary and a vline at
#' the true population r. Panels are faceted (not stacked on one y-axis,
#' as nhst-power.R's plot_corridor() does) because the critical-r boundary
#' depends on n and needs its own position in each panel.
#'
#' @param var_x, var_y     the focal pair of tests.
#' @param n_grid            sample sizes to compare (one facet panel each).
#' @param true_corr         the single population *observed*-score
#'                           correlation to simulate from (0 for the
#'                           H0-true funnel, e.g. .5 for an H1-true funnel).
#' @param reps              independent studies simulated per sample size.
#' @param reliability       reliability override applied to both tests;
#'                           NULL (default) uses Vermeiren-default values.
#' @param alpha             significance threshold.
#' @param seed_base         base seed, passed through to run_sweep().
#' @return a ggplot object.
plot_funnel <- function(var_x, var_y, n_grid, true_corr, reps = 300, reliability = NULL,
                         alpha = 0.05, seed_base = 1234) {
  raw <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr, reliability_grid = reliability,
                    reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  raw$significant_label <- significance_label(raw$significant)
  crit_df <- data.frame(n = n_grid, crit = vapply(n_grid, critical_r, numeric(1), alpha = alpha))

  # Conventional funnel-plot orientation: precision (n) increasing upward,
  # so the narrow, high-precision studies sit at the top and the wide,
  # more variable small-n studies sit at the bottom.
  n_levels <- sort(unique(n_grid), decreasing = TRUE)
  raw$n <- factor(raw$n, levels = n_levels)
  crit_df$n <- factor(crit_df$n, levels = n_levels)

  ggplot(raw, aes(x = estimate, y = 0, color = significant_label)) +
    geom_vline(xintercept = true_corr, color = "gray40", linetype = "dotted") +
    geom_vline(data = crit_df, aes(xintercept = crit), color = "gray40", linetype = "dashed",
               inherit.aes = FALSE) +
    geom_vline(data = crit_df, aes(xintercept = -crit), color = "gray40", linetype = "dashed",
               inherit.aes = FALSE) +
    geom_jitter(height = 0.3, alpha = 0.6, size = 2.2) +
    scale_color_manual(values = sig_colors, name = NULL) +
    facet_wrap(~ n, ncol = 1, labeller = label_both) +
    labs(x = sprintf("Observed r(%s, %s)", var_x, var_y),
         subtitle = sprintf(
           "%d simulated studies per n; true (population) correlation = %.2f (dotted); dashed lines = +/- critical r at alpha = %.2f",
           reps, true_corr, alpha)) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold"))
}

#' As plot_funnel(), but for several population correlations at once --
#' faceted on a 2-D grid of n x true_corr (rows x columns), so true_corr
#' doesn't have to be double-encoded through colour (colour is already
#' spent on significance).
#'
#' @inheritParams plot_funnel
#' @param true_corr_grid  population *observed*-score correlations to
#'                         compare, e.g. c(.25, .5, .75).
#' @return a ggplot object.
plot_funnel_multi_effect <- function(var_x, var_y, n_grid, true_corr_grid, reps = 300,
                                      reliability = NULL, alpha = 0.05, seed_base = 1234) {
  raw <- run_sweep(n_grid = n_grid, true_corr_grid = true_corr_grid, reliability_grid = reliability,
                    reps = reps, var_x = var_x, var_y = var_y, alpha = alpha, seed_base = seed_base)
  raw$significant_label <- significance_label(raw$significant)

  crit_df <- expand.grid(n = n_grid, true_corr = true_corr_grid)
  crit_df$crit <- vapply(crit_df$n, critical_r, numeric(1), alpha = alpha)
  true_df <- data.frame(true_corr = true_corr_grid)

  # Conventional funnel-plot orientation: precision (n) increasing upward.
  n_levels <- sort(unique(n_grid), decreasing = TRUE)
  raw$n <- factor(raw$n, levels = n_levels)
  crit_df$n <- factor(crit_df$n, levels = n_levels)

  ggplot(raw, aes(x = estimate, y = 0, color = significant_label)) +
    geom_vline(data = true_df, aes(xintercept = true_corr), color = "gray40", linetype = "dotted",
               inherit.aes = FALSE) +
    geom_vline(data = crit_df, aes(xintercept = crit), color = "gray40", linetype = "dashed",
               inherit.aes = FALSE) +
    geom_vline(data = crit_df, aes(xintercept = -crit), color = "gray40", linetype = "dashed",
               inherit.aes = FALSE) +
    geom_jitter(height = 0.3, alpha = 0.5, size = 0.9) +
    scale_color_manual(values = sig_colors, name = NULL) +
    facet_grid(n ~ true_corr,
               labeller = labeller(n = label_both, true_corr = function(x) paste0("true r: ", x))) +
    labs(x = sprintf("Observed r(%s, %s)", var_x, var_y),
         subtitle = sprintf(
           "%d simulated studies per cell; dotted line = true (population) correlation; dashed lines = +/- critical r at alpha = %.2f",
           reps, alpha)) +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold"))
}
