# Rebuild the static shinylive site from app-source/.
#
# Run from this folder (the repo root):   Rscript rebuild-site.R
# Needs the 'shinylive' package. The first run downloads the shinylive web assets.
#
# What it does: exports app-source/ into a temporary folder, deletes the previously
# exported site files from this folder, and copies the new ones in. app-source/,
# README.md, .nojekyll and .git are left alone.

stopifnot("run this from the folder that contains app-source/" = file.exists("app-source/app.R"))

tmp <- tempfile("site")
shinylive::export("app-source", tmp)

old <- c("index.html", "app.json", "edit", "shinylive", "shinylive-sw.js")
unlink(old, recursive = TRUE)

ok <- file.copy(list.files(tmp, full.names = TRUE, all.files = TRUE, no.. = TRUE), ".", recursive = TRUE)
stopifnot("copying the new site failed" = all(ok))
unlink(tmp, recursive = TRUE)
message("Site rebuilt. Commit and push to publish.")
