# HLS L30 B04: the R plan table vs odc-stac vs whole-tile gdalwarp.
# Run from the repo root after hls/search.py and hls/odc_reference.py, in a
# shell with the Earthdata GDAL settings (hls/README.md).
#   Rscript hls/run_hls.R            # uses hls/items.json, ref_hls/
#   HLS_ITEMS=... HLS_OUT=... Rscript hls/run_hls.R
suppressMessages(library(gdalraster))
source("R/plan.R")
# LP DAAC returns the odd 502/503: let GDAL retry rather than fail the run
if (!nzchar(Sys.getenv("GDAL_HTTP_MAX_RETRY"))) Sys.setenv(GDAL_HTTP_MAX_RETRY = "5")
if (!nzchar(Sys.getenv("GDAL_HTTP_RETRY_DELAY"))) Sys.setenv(GDAL_HTTP_RETRY_DELAY = "2")

items_file <- Sys.getenv("HLS_ITEMS", "hls/items.json")
out <- Sys.getenv("HLS_OUT", "ref_hls")
fc <- jsonlite::fromJSON(items_file, simplifyVector = FALSE)
meta <- jsonlite::fromJSON(file.path(out, "meta.json"))
vsi <- function(h) if (grepl("^https?://", h)) paste0("/vsicurl/", h) else h

t0 <- proc.time()[["elapsed"]]
pop <- population(fc$features, asset = "B04", header = TRUE, href_fn = vsi)
cat(sprintf("population: %d items, headers read in %.1f s\n", nrow(pop), proc.time()[["elapsed"]] - t0))
pop$solar_day <- solar_day(pop$datetime, pop$lon_centre)
print(pop[, c("item", "src_crs", "grid_from", "dtype", "nodata", "block_x", "compression",
              "predictor", "lon_centre", "solar_day")], row.names = FALSE)

# GeoBoxes: exactly the grids odc-stac used (from its meta.json)
gb_from_meta <- function(m) {
  a <- m$affine   # a, b, c, d, e, f
  geobox(m$crs, c(a[3], a[1], a[2], a[6], a[4], a[5]), c(m$shape[3], m$shape[2]))
}
cases <- intersect(c("warp", "aligned", "aligned_s"), sub("_nearest$", "", names(meta)))

results <- list()
for (case in cases) {
  gbox <- gb_from_meta(meta[[paste0(case, "_nearest")]])
  gbox$crs <- short_crs(gbox$crs)
  t0 <- proc.time()[["elapsed"]]
  plan <- geobox_plan(pop, gbox, chunk = c(2048L, 2048L), pad = 2L)
  t_plan <- proc.time()[["elapsed"]] - t0
  cat(sprintf("\n== %s (%s): plan %d rows in %.2f s; engines: %s\n", case, gbox$crs, nrow(plan), t_plan,
              paste(names(table(plan$engine)), table(plan$engine), sep = "=", collapse = " ")))
  print(plan[, c("chunk", "group", "fuse_order", "item", "engine",
                 "src_xoff", "src_yoff", "src_xsize", "src_ysize")], row.names = FALSE)
  write.csv(plan, file.path(out, sprintf("plan_%s.csv", case)), row.names = FALSE)
  for (rs in c("nearest", "bilinear")) {
    gr <- if (rs == "nearest") "near" else "bilinear"
    t0 <- proc.time()[["elapsed"]]
    r <- execute_plan(plan, gbox, gr)
    t_exec <- proc.time()[["elapsed"]] - t0
    refs <- list(`odc-stac` = sprintf("%s_%s.tif", case, rs),
                 `odc-stac one chunk` = sprintf("%s_%s_onechunk.tif", case, rs))
    refs <- refs[file.exists(file.path(out, unlist(refs)))]
    others <- lapply(refs, function(f) read_stack(file.path(out, f)))
    others[["gdalwarp whole"]] <- gdal_whole(pop, gbox, gr)
    groups <- dimnames(r)[[3]]
    for (nm in names(others)) {
      o <- others[[nm]]
      if (length(dim(o)) == 2) o <- array(o, c(dim(o), 1))
      if (dim(o)[3] != length(groups)) {
        cat(sprintf("   %s has %d groups, R plan has %d: not compared\n", nm, dim(o)[3], length(groups)))
        next
      }
      for (g in seq_along(groups)) {
        results[[length(results) + 1L]] <- cbind(case = case, resampling = rs, group = groups[g],
                                                 b = nm, compare(r[, , g], o[, , g]))
      }
    }
    cat(sprintf("   executed %s in %.1f s\n", rs, t_exec))
  }
}
res <- do.call(rbind, results)
print(res, row.names = FALSE)
write.csv(res, file.path(out, "comparison.csv"), row.names = FALSE)
