#' ============================================================
#' nhst-sampling-variability.R
#'
#' Shiny-free helpers behind the sampling-variability comparison app
#' (nhst-sampling-variability-app/app.R). The app lets a student compare
#' ONE observed (simulated) dataset's r and t against TWO representations
#' of sampling variability:
#'
#'   column 1: the observed ("focal") dataset, large scatterplot, r / t / p
#'   column 2: six more datasets from the same population (same r, same N)
#'   column 3: two piles of dots, 100 studies each, one dot per study:
#'             upper pile = H0 true (r = 0), lower pile = HA true (r = the
#'             population r that generated the focal dataset); dots blue =
#'             not significant, orange = significant at alpha
#'
#' Everything is built on the project's existing machinery -- nothing here
#' modifies it:
#'   psych-battery-dgp.R: simulate_pair(), run_sweep(), default_* helpers
#'   nhst-fisher-neyman-pearson-demo.R: simulate_many_studies(),
#'     plot_scatter_grid(), plot_dotpile_significance(), t_from_r(),
#'     critical_r(), sig_colors
#'
#' Focal pair, as throughout the lecture: ART3 (reading experience, x-axis)
#' and StuVoc1 (vocabulary, y-axis).
#'
#' Requires (sourced first): psych-battery-dgp.R,
#' nhst-fisher-neyman-pearson-demo.R; packages ggplot2, ggdist, MASS.
#'
#' Console usage:
#'   ```{r}
#'   source("psych-battery-dgp.R")
#'   source("nhst-fisher-neyman-pearson-demo.R")
#'   source("nhst-sampling-variability.R")
#'   fg <- simulate_focal_and_grid(n = 30, r = .3, seed = 1)
#'   pl <- simulate_piles(n = 30, r = .3, pile_seed = 1)
#'   plot_focal_scatter(fg)
#'   plot_comparison_grid(fg)
#'   plot_pile(pl, world = "h0", observed_r = fg$focal_stats$r)
#'   plot_pile(pl, world = "ha", observed_r = fg$focal_stats$r)
#'   ```
#' ============================================================

library(ggplot2)

# ---- constants --------------------------------------------------------------------

#' Focal pair: which simulate_pair() test goes on which axis.
sv_var_x <- "ART3"      # reading experience (x-axis)
sv_var_y <- "StuVoc1"   # vocabulary (y-axis)

#' Sample sizes offered by the app: includes the lecture's 30 / 100 / 1000.
sv_n_choices <- c(10, 20, 30, 50, 100, 200, 500, 1000)

#' Population-r range offered by the app. The upper bound is set by the
#' attenuation ceiling: with reliabilities .88 (StuVoc1) and .90 (ART3) the
#' largest achievable observed r is sqrt(.88 * .90), about .89, and beyond
#' it simulate_pair() caps and warns (see disattenuate()).
sv_r_limits <- c(-0.85, 0.85)

#' Scatterplot axis limits: the full score ranges of the two tests (50 and
#' 90 items), fixed so that redraws and panels are comparable.
sv_x_limits <- c(0, 90)
sv_y_limits <- c(0, 50)

#' Colours: points as in the lecture's scatterplots; the observed-r line as
#' in the lecture's "our observed sample" marker.
sv_point_color    <- "#3B6EA5"
sv_observed_color <- "#B23A48"

# ---- labels and statistics -------------------------------------------------------

#' A plain-language size/sign label for a population r, using the Cohen
#' benchmark *nearest* to |r| (cut points at .05, .20, .40), the same
#' convention as cohen_category() in student-datasets.R.
#'
#' @param r  a single population correlation.
#' @return e.g. "positive, medium", or "no correlation" when r = 0.
r_size_label <- function(r) {
  if (abs(r) < 1e-9) return("no correlation")
  size <- as.character(cut(abs(r), breaks = c(-Inf, 0.05, 0.20, 0.40, Inf), right = FALSE,
                           labels = c("negligible", "small", "medium", "large")))
  paste0(if (r > 0) "positive, " else "negative, ", size)
}

#' APA-style exact p: "p = .003", or "p < .001" below .001.
format_p <- function(p) {
  ifelse(p < 0.001, "p < .001", sub("= 0\\.", "= .", sprintf("p = %.3f", p)))
}

#' r, t, df and exact two-sided p for a correlation r from n people.
#'
#' @param r  observed correlation(s).
#' @param n  sample size.
#' @param alpha significance threshold (two-sided).
#' @return a data frame with `r`, `t`, `df`, `p`, `significant`.
study_stats <- function(r, n, alpha = 0.05) {
  df <- n - 2
  t <- t_from_r(r, n)
  p <- 2 * stats::pt(-abs(t), df)
  data.frame(r = r, t = t, df = df, p = p, significant = p < alpha)
}

#' Annotation text: r and t (and optionally df and exact p).
stat_label <- function(stats_row, with_p = FALSE) {
  base <- sprintf("r = %.2f, t(%d) = %.2f", stats_row$r, as.integer(stats_row$df), stats_row$t)
  if (with_p) paste0(base, ", ", format_p(stats_row$p)) else base
}

#' The power of the two-sided test of r = 0 at sample size n when the
#' population correlation is r (equal to `pwr::pwr.r.test()`'s calculation,
#' written in base R to avoid the dependency). At r = 0 this returns alpha.
expected_power <- function(n, r, alpha = 0.05) {
  stopifnot("`n` must be at least 4" = n >= 4)
  if (abs(r) < 1e-9) return(alpha)
  tcrit <- stats::qt(alpha / 2, df = n - 2, lower.tail = FALSE)
  rc  <- sqrt(tcrit^2 / (tcrit^2 + n - 2))
  zr  <- atanh(r) + r / (2 * (n - 1))
  zrc <- atanh(rc)
  stats::pnorm((zr - zrc) * sqrt(n - 3)) + stats::pnorm((-zr - zrc) * sqrt(n - 3))
}

# ---- simulation -------------------------------------------------------------------

#' Simulate the focal dataset and the six comparison datasets.
#'
#' All seven are drawn from the same population (correlation `r`, size `n`).
#' The focal dataset uses `seed`; the six use `seed + 20000` (plus 1..6
#' inside simulate_many_studies()). Seeds stay below 1e6, so they can never
#' coincide with a run_sweep() seed (always at least 1e6) -- see
#' simulate_piles().
#'
#' @param n,r,seed  sample size, population observed-score r, seed.
#' @param alpha     significance threshold.
#' @return a list: `n`, `r`, `focal` (data frame, columns `ART3`, `StuVoc1`),
#'   `focal_stats` (one row), `grid` (long data frame with a `study` factor
#'   labelled with each study's r and t), `grid_stats` (one row per study).
simulate_focal_and_grid <- function(n, r, seed, alpha = 0.05) {
  focal <- simulate_pair(n, sv_var_x, sv_var_y, true_corr = r, seed = seed)[c(sv_var_x, sv_var_y)]
  focal_stats <- study_stats(stats::cor(focal[[sv_var_x]], focal[[sv_var_y]]), n, alpha)

  grid <- simulate_many_studies(n = n, n_studies = 6, true_corr = r,
                                var_x = sv_var_x, var_y = sv_var_y, seed = seed + 20000)
  r_grid <- vapply(levels(grid$study), function(s) {
    g <- grid[grid$study == s, ]
    stats::cor(g[[sv_var_x]], g[[sv_var_y]])
  }, numeric(1))
  grid_stats <- cbind(study = names(r_grid), study_stats(unname(r_grid), n, alpha))
  # Put each study's r and t in its facet strip, so the annotation can never
  # sit on top of the points.
  strip_labels <- sprintf("%s\nr = %.2f, t = %.2f", grid_stats$study, grid_stats$r, grid_stats$t)
  grid$study <- factor(grid$study, levels = levels(grid$study), labels = strip_labels)
  list(n = n, r = r, focal = focal, focal_stats = focal_stats, grid = grid, grid_stats = grid_stats)
}

#' Simulate the two 100-study piles: H0 true (population r = 0) and HA true
#' (population r = `r`), both at sample size `n`, using run_sweep().
#'
#' Seed design: run_sweep() derives each replicate's seed from
#' `seed_base + rep*1e6 + n*1e3 + round(reliability*100) + round((true_corr + 1)*1000)`.
#' If both piles shared a `seed_base`, they would be IDENTICAL at r = 0
#' (same seeds, same population). The HA pile therefore uses
#' `pile_seed + 12345`: the true-corr term differs between the two piles by
#' at most 850, so the offset can never cancel and no seed can coincide.
#'
#' @param n,r        sample size, HA population r.
#' @param pile_seed  base seed for the piles (H0 uses it as is).
#' @param reps       studies per pile (default 100).
#' @param alpha      significance threshold.
#' @return a list: `n`, `r`, `alpha`, `h0` and `ha` (run_sweep() output, one
#'   row per study), `crit_r`, `expected_h0`, `expected_ha`.
simulate_piles <- function(n, r, pile_seed, reps = 100, alpha = 0.05) {
  h0 <- run_sweep(n_grid = n, true_corr_grid = 0, reps = reps, var_x = sv_var_x, var_y = sv_var_y,
                  alpha = alpha, seed_base = pile_seed)
  ha <- run_sweep(n_grid = n, true_corr_grid = r, reps = reps, var_x = sv_var_x, var_y = sv_var_y,
                  alpha = alpha, seed_base = pile_seed + 12345)
  list(n = n, r = r, alpha = alpha, h0 = h0, ha = ha, crit_r = critical_r(n, alpha),
       expected_h0 = alpha, expected_ha = expected_power(n, r, alpha))
}

# ---- plots ---------------------------------------------------------------------------

#' Column 1: the large scatterplot of the focal dataset, annotated with r,
#' t(df) and the exact p. Styled like the lecture's
#' `fig-scatterplot-experience-vocab`; axes fixed to the tests' score ranges.
#'
#' @param fg  a simulate_focal_and_grid() result.
plot_focal_scatter <- function(fg) {
  ggplot(fg$focal, aes(.data[[sv_var_x]], .data[[sv_var_y]])) +
    geom_point(color = sv_point_color, alpha = 0.8, size = 3) +
    scale_x_continuous(limits = sv_x_limits) +
    scale_y_continuous(limits = sv_y_limits) +
    labs(x = "Reading experience", y = "Vocabulary",
         subtitle = stat_label(fg$focal_stats, with_p = TRUE)) +
    theme_bw(base_size = 14) +
    theme(plot.subtitle = element_text(face = "bold", size = 15))
}

#' Column 2: the six comparison datasets as a 2 x 3 grid, each strip
#' annotated with that study's r and t. Built on plot_scatter_grid().
#'
#' @param fg  a simulate_focal_and_grid() result.
plot_comparison_grid <- function(fg) {
  plot_scatter_grid(fg$grid, sv_var_x, sv_var_y, ncol = 3) +
    scale_x_continuous(limits = sv_x_limits, breaks = c(0, 45, 90)) +
    scale_y_continuous(limits = sv_y_limits) +
    labs(x = "Reading experience", y = "Vocabulary") +
    theme(strip.text = element_text(face = "bold", size = 8.5))
}

#' The x-axis range shared by both piles. Default: the full r scale
#' (-1 to 1), so the narrowing of the piles with sample size is visible.
#' With `zoom = TRUE`: the range of everything drawn (both piles, the
#' observed r, the critical-r lines), padded by 5%, still shared by both
#' piles.
pile_xlim <- function(piles, observed_r, zoom = FALSE) {
  if (!zoom) return(c(-1, 1))
  rng <- range(c(piles$h0$estimate, piles$ha$estimate, observed_r, piles$crit_r, -piles$crit_r))
  pad <- 0.05 * diff(rng)
  c(max(-1, rng[1] - pad), min(1, rng[2] + pad))
}

#' Column 3: one pile of 100 dots (one dot per simulated study; x = observed
#' r), blue = not significant, orange = significant, with the world's
#' population r (dotted), the two critical-r boundaries (dashed), and the
#' focal dataset's observed r (solid). Built on plot_dotpile_significance().
#'
#' @param piles       a simulate_piles() result.
#' @param world       "h0" (population r = 0) or "ha" (population r = piles$r).
#' @param observed_r  the focal dataset's observed r.
#' @param zoom        see pile_xlim().
plot_pile <- function(piles, world = c("h0", "ha"), observed_r, zoom = FALSE) {
  world <- match.arg(world)
  sweep <- piles[[world]]
  true_r <- if (world == "h0") 0 else piles$r
  line_name <- if (world == "h0") "r = 0 (H0)" else sprintf("true r = %.2f", true_r)
  expected <- if (world == "h0") piles$expected_h0 else piles$expected_ha
  k <- sum(sweep$significant)
  xlim <- pile_xlim(piles, observed_r, zoom)

  title <- if (world == "h0") {
    "If H0 is true (population r = 0)"
  } else {
    sprintf("If HA is true (population r = %.2f)", true_r)
  }
  subtitle <- sprintf("%d of %d studies significant (about %.0f expected)", k, nrow(sweep),
                      round(expected * nrow(sweep)))

  obs_hjust <- if (observed_r > mean(xlim)) 1.05 else -0.05
  plot_dotpile_significance(sweep$estimate, xlab = "Observed r in each simulated study",
                            significant = sweep$significant,
                            vlines = stats::setNames(true_r, line_name)) +
    geom_vline(xintercept = c(-piles$crit_r, piles$crit_r), color = "gray40", linetype = "dashed") +
    geom_vline(xintercept = observed_r, color = sv_observed_color, linewidth = 1) +
    annotate("text", x = observed_r, y = Inf, label = "our observed r",
             vjust = 3.2, hjust = obs_hjust, size = 3.5, color = sv_observed_color) +
    # Headroom above the dots, so the line labels never sit on top of them
    scale_y_continuous(limits = c(0, 1.4), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim) +
    labs(title = title, subtitle = subtitle) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
}
