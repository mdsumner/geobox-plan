# Chunked plan vs one whole-tile gdalwarp, both with the exact transformer (-et 0).
# Run from the repo root after R/run.R (reads ref/plan_warp.rds).
suppressMessages(library(gdalraster)); source("R/plan.R")
fc <- jsonlite::fromJSON("fixtures/items.json", simplifyVector = FALSE); pop <- population(fc$features)
tile <- geobox("EPSG:4326", c(-180 + 1195000 * 0.0003, 0.0003, 0, -90 + 245000 * 0.0003, 0, 0.0003), c(5000L, 5000L))
plan <- readRDS("ref/plan_warp.rds")
whole <- function(srcs, rs) { dst <- create("MEM","",5000,5000,1,"UInt16",return_obj=TRUE)
  dst$setGeoTransform(tile$gt); dst$setProjection(srs_to_wkt("EPSG:4326")); dst$setNoDataValue(1,0); dst$fillRaster(1,0,0)
  warp(srcs, dst, "", cl_arg=c("-r",rs,"-srcnodata",0,"-dstnodata",0,"-et",0), quiet=TRUE)
  m <- matrix(dst$read(1,0,0,5000,5000,5000,5000),5000,byrow=TRUE); m[is.na(m)] <- 0L; m }
gs <- sort(unique(plan$group))
for (rs in c("near","bilinear")) for (g in seq_along(gs)) {
  r <- execute_plan(plan[plan$group == gs[g], ], tile, rs, c("-et", 0))[,,1]
  w <- whole(rev(pop$href[solar_day(pop$datetime, pop$lon_centre) == gs[g]]), rs)
  d <- abs(as.numeric(r) - as.numeric(w))
  cat(sprintf("%-8s %s  exact transformer both sides: differing %d, max %g\n", rs, gs[g], sum(d > 0), max(d)))
}
