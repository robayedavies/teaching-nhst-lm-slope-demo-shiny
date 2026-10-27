# Sampling-variability app, linear-model VERSION B (one predictor: a test, or gender)
#
# Compares ONE observed (simulated) dataset's slope b and R^2 with two representations of
# sampling variability: four more datasets from the same population, and two piles of
# 100 simulated studies (H0 true / HA true). All logic lives in
# ../nhst-sampling-variability-lm-b.R; this file is only the interface.
#
# Version B = the Version C app (nhst-sampling-variability-lm-app/) without Phase C: no age,
# no "Also adjust for" boxes, no added-variable plots. See
# nhst-sampling-variability-lm-b-walkthrough.qmd.
#
# Run from the project root with:
#   shiny::runApp("nhst-sampling-variability-lm-b-app")
#
# Requires: shiny, bslib, ggplot2, ggdist, MASS.

library(shiny)
library(bslib)
library(ggplot2)

# Helper scripts: looked for next to this file first (for a deployed bundle), then in the
# project root (the folder above this one).
src <- function(f) source(if (file.exists(f)) f else file.path("..", f))
src("psych-battery-dgp.R")
src("nhst-fisher-neyman-pearson-demo.R")
src("nhst-sampling-variability.R")
src("nhst-sampling-variability-lm-b.R")

random_seed <- function() sample.int(10000, 1)

# Fewer value labels under the slider: the default labelling crowds the largest values at
# the right-hand end. data-grid-num is the ionRangeSlider option for the number of grid
# intervals (4 intervals = 5 labels: -max, -max/2, 0, max/2, max). shiny sets its own value
# first, so that one is removed (adding would just append a second value to the attribute).
few_ticks <- function(slider, n = 4) {
  htmltools::tagQuery(slider)$find("input.js-range-slider")$
    removeAttrs("data-grid-num")$addAttrs(`data-grid-num` = n)$allTags()
}

x_default <- "ART3"
y_default <- "StuVoc1"
limit_default <- lm_slope_limit(x_default, y_default)

# Slider label for the current predictor / units.
slider_label <- function(x, units) {
  if (units == "std") return("Population standardised beta (= r)")
  if (lm_is_binary(x)) "Population mean difference b (M - F, points)" else "Population slope b (raw units)"
}

# Choices offered for the outcome: with a test as predictor the outcome is the other test (one
# choice, so the box shows it but cannot change it); with gender either test can be the outcome.
y_choices <- function(x) if (lm_is_binary(x)) lm_var_choices else lm_var_choices[lm_var_choices == lm_other_var(x)]

# The outcome: the other test, unless the predictor is gender (then the Outcome selector decides).
outcome_of <- function(x, y_input) if (lm_is_binary(x)) y_input else lm_other_var(x)

# Largest slider value in standardised units, rounded down to the slider's step.
std_slider_limit <- function() floor(sv_r_limits[2] / 0.05 + 1e-9) * 0.05


# ---- User interface ------------------------------------------------------------------

ui <- page_fillable(
  title = "Sampling variability: the slope and R-squared of a linear model",
  theme = bs_theme(version = 5),
  padding = 12,
  gap = 12,

  # Compact control card: tighter spacing between inputs and around the card's contents
  tags$head(tags$style(HTML("
    .ctl .card-body { padding: 0.4rem 0.75rem; gap: 0; }
    .ctl .shiny-input-container { margin-bottom: 0.2rem; }
    .ctl .form-label, .ctl .control-label { margin-bottom: 0.1rem; }
    .ctl .form-check { margin-bottom: 0; min-height: 1.2rem; }
    .ctl .irs { margin-top: 0.1rem; }
  "))),

  card(
    fill = FALSE,
    class = "py-0 ctl",
    layout_columns(
      col_widths = c(2, 3, 3, 2, 2),
      selectizeInput("n", "N per study", choices = sv_n_choices, selected = 30,
                     width = "100%", options = list(dropdownParent = "body")),
      selectInput("x_var", "Predictor (X)", choices = lm_x_choices, selected = x_default,
                  width = "100%"),
      selectInput("y_var", "Outcome (Y)", choices = y_choices(x_default), selected = y_default,
                  width = "100%"),
      numericInput("seed", "Seed", value = random_seed(), min = 1, max = 10000, step = 1,
                   width = "100%"),
      div(
        class = "d-flex flex-column justify-content-center h-100 gap-1",
        actionButton("new_samples", "Draw new samples", class = "btn-primary btn-sm w-100"),
        actionButton("new_piles", "Redraw the piles", class = "btn-outline-secondary btn-sm w-100"),
        checkboxInput("zoom", "Zoom to piles", value = FALSE, width = "100%")
      )
    ),
    layout_columns(
      class = "mt-3",
      col_widths = c(8, 4),
      # Slider package: label and units buttons on one line, the slider below
      div(
        div(class = "d-flex justify-content-between align-items-baseline gap-3",
            tags$label(class = "control-label fw-bold", `for` = "beta", textOutput("slider_title", inline = TRUE)),
            radioButtons("units", NULL, choices = c("Raw units" = "raw", "Standardised" = "std"),
                         selected = "raw", inline = TRUE)),
        few_ticks(sliderInput("beta", NULL, min = -limit_default, max = limit_default,
                              value = 0.2, step = 0.01, width = "100%"))
      ),
      # Under the seed box: model and statistic
      div(
        div(class = "small text-muted mb-1", uiOutput("model_note")),
        radioButtons("stat", NULL, choices = c("Slope b" = "beta", "R\u00b2" = "r2"),
                     selected = "beta", inline = TRUE)
      )
    )
  ),

  # Fills the rest of the window; never shorter than min_height (the page scrolls below that)
  layout_columns(
    min_height = "560px",
    col_widths = c(4, 3, 5),
    card(
            card_header(class = "py-1", "1. Our observed dataset"),
      plotOutput("focal"),
      card_footer(class = "small py-1", uiOutput("focal_note"))
    ),
    card(
            card_header(class = "py-1", "2. Four more datasets, same population"),
      plotOutput("grid"),
      card_footer(class = "small py-1", uiOutput("grid_note"))
    ),
    card(
            card_header(class = "py-1", "3. Where does our estimate fall? 100 studies per world"),
      plotOutput("pile_h0"),
      plotOutput("pile_ha"),
      card_footer(class = "small py-1", uiOutput("pile_note"))
    )
  )
)


# ---- Server --------------------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(pile_seed = random_seed())

  observeEvent(input$new_samples, updateNumericInput(session, "seed", value = random_seed()))
  observeEvent(input$new_piles, rv$pile_seed <- random_seed())

  x_val <- reactive({ req(input$x_var); input$x_var })
  y_val <- reactive(outcome_of(x_val(), input$y_var))
  n_val <- reactive(as.integer(input$n))

  # The outcome box is always shown. With a test as predictor it holds only the other test;
  # with gender it offers both tests (keeping the current choice).
  observeEvent(input$x_var, {
    ch <- y_choices(input$x_var)
    updateSelectInput(session, "y_var", choices = ch,
                      selected = if (input$y_var %in% ch && lm_is_binary(input$x_var)) input$y_var else
                        if (lm_is_binary(input$x_var)) ch[[1]] else lm_other_var(input$x_var))
  }, ignoreInit = TRUE)

  output$slider_title <- renderText(slider_label(x_val(), input$units))

  # Changing the predictor, the outcome (gender only) or the units re-ranges the slider
  # (units and limits change). Keep the same standardised slope rather than the same number.
  prev <- reactiveValues(x = x_default, y = y_default, units = "raw")
  observeEvent(list(input$x_var, input$y_var, input$units), {
    x <- input$x_var; y <- outcome_of(x, input$y_var)
    req(x != y)
    if (identical(prev$x, x) && identical(prev$y, y) && identical(prev$units, input$units)) return()
    r_keep <- if (prev$units == "std") input$beta else lm_beta_to_r(input$beta, prev$x, prev$y)
    if (input$units == "std") {
      lim <- std_slider_limit()
      updateSliderInput(session, "beta",
                        min = -lim, max = lim, step = 0.05,
                        value = max(-lim, min(lim, round(r_keep * 20) / 20)))
    } else {
      lim <- lm_slope_limit(x, y); st <- lm_slope_step(x)
      v <- round(round(lm_r_to_beta(r_keep, x, y) / st) * st, 2)
      updateSliderInput(session, "beta", min = -lim, max = lim,
                        step = st, value = max(-lim, min(lim, v)))
    }
    prev$x <- x; prev$y <- y; prev$units <- input$units
  }, ignoreInit = TRUE)

  seed_val <- reactive({
    s <- input$seed
    validate(need(!is.null(s) && !is.na(s) && s >= 1 && s <= 10000 && s == round(s),
                  "Seed must be a whole number from 1 to 10000."))
    as.integer(s)
  })

  # The slider's value is in the chosen units; internally everything uses the raw slope.
  # While the slider is being re-ranged its value can briefly lie outside the new range
  # (or in the old units); wait for the update rather than simulate an impossible r.
  beta_val <- reactive({
    v <- input$beta
    req(!is.null(v), input$units %in% c("raw", "std"))
    lim <- if (input$units == "std") std_slider_limit() else lm_slope_limit(x_val(), y_val())
    req(abs(v) <= lim + 1e-9)
    if (input$units == "std") lm_r_to_beta(v, x_val(), y_val()) else v
  })

  # cache = "session": the same app.R also runs as a static shinylive page
  fg <- reactive(simulate_lm_focal_and_grid(n_val(), beta_val(), x_val(), y_val(), seed_val())) |>
    bindCache(n_val(), beta_val(), x_val(), y_val(), seed_val(), cache = "session")

  piles <- reactive(simulate_lm_piles(n_val(), beta_val(), x_val(), y_val(), rv$pile_seed)) |>
    bindCache(n_val(), beta_val(), x_val(), y_val(), rv$pile_seed, cache = "session")

  pile_for_app <- function(world) {
    # Skip the very first draw, before the browser has given the plot a size
    h <- session$clientData[[paste0("output_pile_", world, "_height")]]
    req(!is.null(h), h > 30)
    w <- session$clientData[[paste0("output_pile_", world, "_width")]]
    # Plotting-panel size = plot size minus the title / subtitle / axis space; sets the dot size
    plot_lm_pile(piles(), world, input$stat, fg()$focal_stats, zoom = input$zoom, units = input$units,
                 panel_px = c(max(w - 50, 200), max(h - 100, 50))) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 11),
            plot.margin = margin(5.5, 10, 14, 5.5))
  }

  output$focal <- renderPlot(plot_lm_focal(fg(), input$units), res = 96)
  output$grid  <- renderPlot(plot_lm_grid(fg(), input$units), res = 96)
  output$pile_h0 <- renderPlot(pile_for_app("h0"), res = 96)
  output$pile_ha <- renderPlot(pile_for_app("ha"), res = 96)

  # Model formula
  output$model_note <- renderUI(span(lm_var_label(y_val()), " ~ ", lm_var_label(x_val())))

  output$focal_note <- renderUI({
    r <- lm_beta_to_r(beta_val(), x_val(), y_val())
    if (lm_is_binary(x_val())) {
      helpText(sprintf(paste("HA world: population mean difference b = %.1f points (M - F; F is the",
                             "reference group, coded 0), i.e. Cohen's d = %.2f, point-biserial r = standardised",
                             "beta = %.2f (%s), population R\u00b2 = %.2f; N = %d (equal groups).",
                             "Same t and p as the equal-variance two-sample t test."),
                       beta_val(), 2 * r / sqrt(1 - r^2), r, r_size_label(r), r^2, n_val()))
    } else {
      helpText(sprintf(paste("HA world: population slope b = %.2f (raw units), i.e. standardised beta = r = %.2f",
                             "(%s) and population R\u00b2 = %.2f; N = %d. Same t and p as the test of r."),
                       beta_val(), r, r_size_label(r), r^2, n_val()))
    }
  })

  output$grid_note <- renderUI({
    helpText("Same population slope, same N, independent samples: each gives a different b and R\u00b2.")
  })

  output$pile_note <- renderUI({
    base <- tagList(
      span(style = "color:#0072B2; font-weight:bold", "Blue"), ": not significant; ",
      span(style = "color:#E69F00; font-weight:bold", "orange"), ": significant. ",
      "Dashed: critical value (approximate: each study's own SE differs; colours use its exact p). ",
      span(style = "font-weight:bold", "Solid black line"), ": our observed value (dotted: true value).")
    extra <- if (input$stat == "r2") {
      paste("R\u00b2 cannot be negative, so even under H0 it averages about 1/(N-1), not 0;",
            "the significance test is two-sided in b but one line in R\u00b2.")
    } else if (abs(beta_val()) < 1e-9) {
      "At b = 0, H0 and HA describe the same population, so both piles are draws from it."
    } else NULL
    helpText(base, extra)
  })
}

shinyApp(ui, server)
