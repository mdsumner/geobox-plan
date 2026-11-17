# Side by side: the R plan table vs odc-stac on the same GeoBoxes.
# Run from the repo root after python/odc_reference.py.
suppressMessages(library(gdalraster))
source("R/plan.R")

fc <- jsonlite::fromJSON("fixtures/items.json", simplifyVector = FALSE)
pop <- population(fc$features)
meta <- jsonlite::fromJSON("ref/meta.json")

# Case 1: WGS84GRID30[49, 239]. Integer arithmetic gives the same origin
# odc-geo does (178.49999999999994); written out here as the literal grid.
tile <- geobox("EPSG:4326", c(-180 + 1195000 * 0.0003, 0.0003, 0, -90 + 245000 * 0.0003, 0, 0.0003),
               c(5000L, 5000L))
# Case 2: 1024 x 1024 on scene C's own grid, on its 512 blocks.
cgt <- unlist(pop[pop$item == "C", paste0("sgt", 1:6)])
aligned <- geobox("EPSG:32760", c(cgt[1] + 1024 * cgt[2], cgt[2], 0, cgt[4] + 2048 * cgt[6], 0, cgt[6]),
                  c(1024L, 1024L))

read_ref <- function(f) {
  ds <- new(GDALRaster, f); on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize(); nb <- ds$getRasterCount()
  simplify2array(lapply(seq_len(nb), function(b) {
    m <- matrix(ds$read(b, 0, 0, nx, ny, nx, ny), ny, byrow = TRUE); m[is.na(m)] <- 0L; m }))
}

# GDAL's own answer: one gdalwarp per group over the whole GeoBox, sources
# listed in reverse fuse order (gdalwarp paints later sources over earlier).
gdal_whole <- function(pop, gbox, resampling) {
  pop$group <- solar_day(pop$datetime, pop$lon_centre)
  gs <- sort(unique(pop$group))
  simplify2array(lapply(gs, function(g) {
    s <- rev(pop$href[pop$group == g])
    dst <- create("MEM", "", gbox$dim[1], gbox$dim[2], 1, "UInt16", return_obj = TRUE)
    dst$setGeoTransform(gbox$gt); dst$setProjection(srs_to_wkt(gbox$crs))
    dst$setNoDataValue(1, 0); dst$fillRaster(1, 0, 0)
    warp(s, dst, "", cl_arg = c("-r", resampling, "-srcnodata", 0, "-dstnodata", 0), quiet = TRUE)
    m <- matrix(dst$read(1, 0, 0, gbox$dim[1], gbox$dim[2], gbox$dim[1], gbox$dim[2]),
                gbox$dim[2], byrow = TRUE)
    m[is.na(m)] <- 0L
    dst$close(); m
  }))
}

compare <- function(a, b) {
  va <- a != 0; vb <- b != 0; both <- va & vb
  d <- abs(as.numeric(a[both]) - as.numeric(b[both]))
  data.frame(valid_a = sum(va), valid_b = sum(vb), only_a = sum(va & !vb), only_b = sum(vb & !va),
             equal = mean(d == 0), max_abs_diff = if (length(d)) max(d) else NA)
}

results <- list()
for (case in c("warp", "aligned")) {
  gbox <- if (case == "warp") tile else aligned
  t0 <- proc.time()[["elapsed"]]
  plan <- geobox_plan(pop, gbox, chunk = c(2048L, 2048L), pad = 2L)
  t_plan <- proc.time()[["elapsed"]] - t0
  cat(sprintf("\n== %s: plan %d rows in %.2f s\n", case, nrow(plan), t_plan))
  print(plan[, c("chunk", "group", "fuse_order", "item", "engine",
                 "src_xoff", "src_yoff", "src_xsize", "src_ysize")], row.names = FALSE)
  saveRDS(plan, sprintf("ref/plan_%s.rds", case))
  write.csv(plan, sprintf("ref/plan_%s.csv", case), row.names = FALSE)
  for (rs in c("nearest", "bilinear")) {
    gr <- if (rs == "nearest") "near" else "bilinear"
    t0 <- proc.time()[["elapsed"]]
    r <- execute_plan(plan, gbox, gr)
    t_exec <- proc.time()[["elapsed"]] - t0
    odc <- read_ref(sprintf("ref/%s_%s.tif", case, rs))
    gw <- gdal_whole(pop, gbox, gr)
    for (g in seq_len(dim(r)[3])) {
      x <- rbind(cbind(case = case, resampling = rs, group = dimnames(r)[[3]][g],
                       odc_time = meta[[paste0(case, "_", rs)]]$times[g], a = "R plan", b = "odc-stac",
                       compare(r[, , g], odc[, , g])),
                 cbind(case = case, resampling = rs, group = dimnames(r)[[3]][g],
                       odc_time = meta[[paste0(case, "_", rs)]]$times[g], a = "R plan", b = "gdalwarp whole",
                       compare(r[, , g], gw[, , g])))
      results[[length(results) + 1L]] <- x
    }
    cat(sprintf("   executed %s in %.1f s\n", rs, t_exec))
  }
}
res <- do.call(rbind, results)
print(res, row.names = FALSE)
write.csv(res, "ref/comparison.csv", row.names = FALSE)
