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
print(pop[, c("item", "src_crs", "grid_from", "dtype", "block_x", "compression", "lon_centre")], row.names = FALSE)
aligned <- geobox("EPSG:32760", c(cgt[1] + 1024 * cgt[2], cgt[2], 0, cgt[4] + 2048 * cgt[6], 0, cgt[6]),
                  c(1024L, 1024L))

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
    odc <- read_stack(sprintf("ref/%s_%s.tif", case, rs))
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
