# Chunked plan vs one whole-tile gdalwarp, both with the exact transformer (-et 0).
# Run from the repo root after R/run.R (reads ref/plan_warp.rds).
suppressMessages(library(gdalraster)); source("R/plan.R")
fc <- jsonlite::fromJSON("fixtures/items.json", simplifyVector = FALSE); pop <- population(fc$features)
tile <- geobox("EPSG:4326", c(-180 + 1195000 * 0.0003, 0.0003, 0, -90 + 245000 * 0.0003, 0, 0.0003), c(5000L, 5000L))
plan <- readRDS("ref/plan_warp.rds")
gs <- sort(unique(plan$group))
for (rs in c("near","bilinear")) for (g in seq_along(gs)) {
  r <- execute_plan(plan[plan$group == gs[g], ], tile, rs, c("-et", 0))[,,1]
  w <- gdal_whole(pop, tile, rs, c("-et", 0))[, , g]
  x <- compare(r, w)
  cat(sprintf("%-8s %s  exact transformer both sides: valid %d / %d, differing %d, max %g\n",
              rs, gs[g], x$valid_a, x$valid_b, round((1 - x$equal) * x$valid_a), x$max_abs_diff))
}
