# Sampling variability of a slope and R-squared (linear model, one predictor)

An interactive teaching demonstration. It fits a simple linear model to ONE
simulated dataset and compares the slope b (or R-squared) with two
representations of sampling variability: four more datasets from the same
population, and two piles of 100 simulated studies (H0 true and HA true). The
predictor is a test (reading experience or vocabulary) or gender.

This is "Version B" of the linear-model app from the `teaching-critical-thinking`
project (one predictor; no adjustment for other predictors).

It is a Shiny app that runs **entirely in the browser** (shinylive / webR), so
this repository is a static website: there is no R server.

Live page:
`https://robayedavies.github.io/teaching-nhst-lm-slope-demo-shiny/`

## Using the page

- The first visit downloads about 60-90 MB (the R engine and packages): allow
  several seconds to a minute. Later visits load from the browser's cache.
- Open it before a lecture and keep the tab open.
- Use **Chrome or Edge**. Safari is not supported (the sister correlation app
  fails there with a shinylive/webR error, posit-dev/r-shinylive issue #204; this app
  has not been tested in Safari and is expected to fail the same way).
- Controls: N per study, Predictor (X) and Outcome (Y), a seed, "Draw new samples",
  "Redraw the piles", "Zoom to piles", Raw/Standardised units, the population slope
  slider, and a Slope b / R-squared toggle for the piles.

## What is in this folder

| Path | What it is |
|---|---|
| `index.html`, `app.json`, `shinylive/`, `edit/`, `shinylive-sw.js` | The exported site (generated; do not edit by hand) |
| `app-source/` | `app.R` plus the four helper scripts it loads. **Edit these, then rebuild** |
| `rebuild-site.R` | Regenerates the site from `app-source/` |
| `.nojekyll` | Tells GitHub Pages not to run Jekyll over the files |

`app-source/` holds copies of files from the `teaching-critical-thinking` project
(`nhst-sampling-variability-lm-b-app/app.R`, `nhst-sampling-variability-lm-b.R`,
`psych-battery-dgp.R`, `nhst-fisher-neyman-pearson-demo.R`,
`nhst-sampling-variability.R`), copied on 2026-10-27 (project commit `3baa5c2`).
The copies are not linked to the originals: re-copy them if the originals change.

## Updating the site

1. Replace the files in `app-source/` with the new versions.
2. From this folder run `Rscript rebuild-site.R` (needs the R package `shinylive`).
3. Test locally: `httpuv::runStaticServer(".")` in R, then open the address it prints.
   A page like this cannot be opened as a plain file.
4. Commit and push. After a push the browser may keep showing the old version: hard
   refresh (Cmd+Shift+R) or use a private window.

## Notes

- Everything in `app-source/` ends up publicly readable on the published page (the
  source is inside `app.json`). Do not put anything private there.
- The app needs `cache = "session"` in its `bindCache()` calls, because shinylive does
  not set up shiny's default app-level cache.
- Exporting prints a warning that the installed `distributional` package (0.8.1) differs
  from the WebAssembly version (0.9.0). The app ran with no errors in testing.

## Browser test log

| Date | Browser | Result |
|---|---|---|
| 2026-10-27 | Chrome (headless), local copy of this site, fresh profile | Works. All four plots after 7.0 s; N = 1000 settled in 2.5 s; "Redraw the piles" 1.5 s; "Draw new samples" 2.1 s (seed and focal plot change); gender predictor works; no output errors |
| 2026-10-27 | Chrome (user), live GitHub Pages page | Works; looks as intended |
| 2026-10-27 | Chrome (headless), live GitHub Pages page, fresh profile | Works. All four plots appear after 33.2 s on a first visit (varies with the network; 7.0 s on localhost), no reload, no output errors. After the first render the page settled in 0.9 s. N = 1000 settled in 2.3 s; "Redraw the piles" 1.4 s (piles changed); "Draw new samples" 2.0 s (seed and focal plot changed); back to N = 30 2.2 s; gender predictor works (slider re-ranged, mean-difference footer). Only the harmless messages: a `favicon.ico` 404, the webR PostMessage warning and the "LazyFiles on gzip" warnings |
