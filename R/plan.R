# A GeoBox plan table: population x target chunks -> one row per read.
#
# The plan is a data frame. Every row says: to fill this target chunk, for
# this group (time step), read this window of this asset, with this engine.
# Building it reads no pixels. Executing it is a loop over rows.
#
# Dependencies: gdalraster (transforms, warp, reads), jsonlite.
# Grid arithmetic is written out here; in a package it is vaster's.

# ---- target grid ------------------------------------------------------------

#' A GeoBox: CRS, GDAL geotransform, dim = c(ncol, nrow).
geobox <- function(crs, gt, dim) {
  structure(list(crs = crs, gt = as.numeric(gt), dim = as.integer(dim)),
            class = "geobox")
}

#' Split a GeoBox into chunks. Chunk origins are computed from integer pixel
#' offsets of the parent, so a chunk's grid is exactly its parent's grid.
geobox_chunks <- function(gbox, chunk = c(2048L, 2048L)) {
  nx <- ceiling(gbox$dim[1] / chunk[1]); ny <- ceiling(gbox$dim[2] / chunk[2])
  g <- expand.grid(ci = seq_len(nx) - 1L, cj = seq_len(ny) - 1L)
  g$xoff <- g$ci * chunk[1]; g$yoff <- g$cj * chunk[2]
  g$ncol <- pmin(chunk[1], gbox$dim[1] - g$xoff)
  g$nrow <- pmin(chunk[2], gbox$dim[2] - g$yoff)
  g$chunk <- sprintf("%d_%d", g$cj, g$ci)
  g[, c("chunk", "ci", "cj", "xoff", "yoff", "ncol", "nrow")]
}

chunk_gt <- function(gbox, xoff, yoff) {
  gt <- gbox$gt
  c(gt[1] + xoff * gt[2] + yoff * gt[3], gt[2], gt[3],
    gt[4] + xoff * gt[5] + yoff * gt[6], gt[5], gt[6])
}

# ---- population -------------------------------------------------------------

#' One row per (item, asset) from STAC items carrying proj:*. Accepts
#' proj:epsg or proj:code. The longitude centre is antimeridian-aware.
population <- function(features, asset = "red") {
  rows <- lapply(features, function(f) {
    p <- f$properties; a <- f$assets[[asset]]
    crs <- if (!is.null(p[["proj:code"]])) p[["proj:code"]] else paste0("EPSG:", p[["proj:epsg"]])
    tr <- unlist(p[["proj:transform"]])          # affine order a,b,c,d,e,f
    gt <- c(tr[3], tr[1], tr[2], tr[6], tr[4], tr[5])
    bb <- unlist(f$bbox)
    xmax <- if (bb[3] < bb[1]) bb[3] + 360 else bb[3]
    lon <- (bb[1] + xmax) / 2; lon <- ((lon + 180) %% 360) - 180
    nd <- a[["raster:bands"]][[1]]$nodata
    data.frame(item = f$id, datetime = as.POSIXct(p$datetime, "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"),
               asset = asset, href = a$href, src_crs = crs,
               sgt1 = gt[1], sgt2 = gt[2], sgt3 = gt[3], sgt4 = gt[4], sgt5 = gt[5], sgt6 = gt[6],
               src_ncol = as.integer(unlist(p[["proj:shape"]])[2]),
               src_nrow = as.integer(unlist(p[["proj:shape"]])[1]),
               nodata = if (is.null(nd)) NA_real_ else nd, lon_centre = lon)
  })
  do.call(rbind, rows)
}

#' odc-stac's "solar_day": UTC time shifted by longitude / 15 hours.
solar_day <- function(datetime, lon) as.Date(datetime + lon / 15 * 3600, tz = "UTC")

# ---- footprint: target into source, never source into target ---------------

# Pixel-corner sample of a chunk: dense edges plus an interior lattice. The
# points go from the target CRS into the source CRS. Projecting the target
# into the source keeps the antimeridian out of it: a UTM scene that
# straddles 180 is continuous in its own CRS, and its outline only breaks
# when it is projected into longitude.
chunk_points <- function(gt, ncol, nrow, n = 17L) {
  i <- unique(c(seq(0, ncol, length.out = n), 0, ncol))
  j <- unique(c(seq(0, nrow, length.out = n), 0, nrow))
  p <- expand.grid(i = i, j = j)
  cbind(gt[1] + p$i * gt[2] + p$j * gt[3], gt[4] + p$i * gt[5] + p$j * gt[6])
}

src_window <- function(xy_src, sgt, src_ncol, src_nrow, pad) {
  ok <- is.finite(xy_src[, 1]) & is.finite(xy_src[, 2])
  if (!any(ok)) return(NULL)
  inv <- gdalraster::inv_geotransform(sgt)
  px <- inv[1] + xy_src[ok, 1] * inv[2] + xy_src[ok, 2] * inv[3]
  py <- inv[4] + xy_src[ok, 1] * inv[5] + xy_src[ok, 2] * inv[6]
  x0 <- max(0, floor(min(px)) - pad); x1 <- min(src_ncol, ceiling(max(px)) + pad)
  y0 <- max(0, floor(min(py)) - pad); y1 <- min(src_nrow, ceiling(max(py)) + pad)
  if (x1 <= x0 || y1 <= y0) return(NULL)
  c(x0, y0, x1 - x0, y1 - y0)
}

# ---- engine choice ----------------------------------------------------------

# "bytes":  same CRS, same pixel size, integer pixel offset, and the window
#           falls on whole source blocks: chunks are fetched and decoded as is.
# "window": same grid, not on blocks: a plain windowed read, no resampling.
# "warp":   everything else.
choose_engine <- function(dst_crs, cgt, ncol, nrow, src, block = 512L, tol = 1e-9) {
  if (!gdalraster::srs_is_same(dst_crs, src$src_crs)) return(list(engine = "warp"))
  sgt <- c(src$sgt1, src$sgt2, src$sgt3, src$sgt4, src$sgt5, src$sgt6)
  if (any(abs(cgt[c(2, 3, 5, 6)] - sgt[c(2, 3, 5, 6)]) > tol)) return(list(engine = "warp"))
  ox <- (cgt[1] - sgt[1]) / sgt[2]; oy <- (cgt[4] - sgt[4]) / sgt[6]
  if (abs(ox - round(ox)) > 1e-6 || abs(oy - round(oy)) > 1e-6) return(list(engine = "warp"))
  ox <- round(ox); oy <- round(oy)
  onblocks <- ox %% block == 0 && oy %% block == 0 &&
    (ncol %% block == 0 || ox + ncol >= src$src_ncol) &&
    (nrow %% block == 0 || oy + nrow >= src$src_nrow)
  list(engine = if (onblocks) "bytes" else "window", ox = ox, oy = oy)
}

# ---- the plan ---------------------------------------------------------------

#' Build the plan table. pad is the resampling kernel's reach in source
#' pixels (1 for nearest, 2 for bilinear/cubic).
geobox_plan <- function(pop, gbox, chunk = c(2048L, 2048L), pad = 2L, block = 512L) {
  ch <- geobox_chunks(gbox, chunk)
  pop$group <- solar_day(pop$datetime, pop$lon_centre)
  pop$fuse_order <- ave(seq_len(nrow(pop)), pop$group, FUN = seq_along)
  out <- list()
  for (k in seq_len(nrow(ch))) {
    cgt <- chunk_gt(gbox, ch$xoff[k], ch$yoff[k])
    pts <- chunk_points(cgt, ch$ncol[k], ch$nrow[k])
    for (s in seq_len(nrow(pop))) {
      src <- pop[s, ]
      sgt <- c(src$sgt1, src$sgt2, src$sgt3, src$sgt4, src$sgt5, src$sgt6)
      xy <- suppressWarnings(gdalraster::transform_xy(pts, gbox$crs, src$src_crs))
      w <- src_window(xy, sgt, src$src_ncol, src$src_nrow, pad)
      if (is.null(w)) next
      e <- choose_engine(gbox$crs, cgt, ch$ncol[k], ch$nrow[k], src, block)
      if (e$engine != "warp") {          # exact window, no pad needed
        w <- c(e$ox, e$oy, min(ch$ncol[k], src$src_ncol - e$ox), min(ch$nrow[k], src$src_nrow - e$oy))
        if (w[3] <= 0 || w[4] <= 0) next
      }
      out[[length(out) + 1L]] <- data.frame(
        chunk = ch$chunk[k], dst_xoff = ch$xoff[k], dst_yoff = ch$yoff[k],
        dst_ncol = ch$ncol[k], dst_nrow = ch$nrow[k],
        group = src$group, fuse_order = src$fuse_order, item = src$item,
        asset = src$asset, href = src$href, src_crs = src$src_crs,
        engine = e$engine, src_xoff = w[1], src_yoff = w[2], src_xsize = w[3], src_ysize = w[4],
        nodata = src$nodata)
    }
  }
  p <- do.call(rbind, out)
  p[order(p$group, p$chunk, p$fuse_order), ]
}

# ---- execution --------------------------------------------------------------

# Byte engine: block offsets and sizes come from the TIFF tags (GDAL exposes
# them as BLOCK_OFFSET_x_y / BLOCK_SIZE_x_y). Here they are read from GDAL
# metadata for brevity; in the real thing they are a byte-ref table
# (tabularberg refs_* rows, a kerchunk Parquet, a get-refs harvest) and GDAL
# is not involved at read time. DEFLATE, no predictor, uint16 only.
read_bytes_window <- function(href, xoff, yoff, xsize, ysize, block = 512L) {
  ds <- new(gdalraster::GDALRaster, href); on.exit(ds$close())
  out <- matrix(0L, ysize, xsize)
  for (bj in seq(yoff %/% block, (yoff + ysize - 1) %/% block)) {
    for (bi in seq(xoff %/% block, (xoff + xsize - 1) %/% block)) {
      off <- as.numeric(ds$getMetadataItem(1, sprintf("BLOCK_OFFSET_%d_%d", bi, bj), "TIFF"))
      sz <- as.numeric(ds$getMetadataItem(1, sprintf("BLOCK_SIZE_%d_%d", bi, bj), "TIFF"))
      con <- file(href, "rb"); seek(con, off); raw <- readBin(con, "raw", sz); close(con)
      v <- readBin(memDecompress(raw, type = "gzip"), "integer", n = block * block,
                   size = 2, signed = FALSE, endian = "little")
      m <- matrix(v, block, block, byrow = TRUE)
      r <- (bj * block - yoff) + seq_len(block); c <- (bi * block - xoff) + seq_len(block)
      keep_r <- r >= 1 & r <= ysize; keep_c <- c >= 1 & c <= xsize
      out[r[keep_r], c[keep_c]] <- m[keep_r, keep_c]
    }
  }
  out
}

# Warp engine: the plan's source window becomes a VRT (gdal_translate
# -srcwin), so the plan, not the warper's own boundary walk, decides which
# source pixels are in play. Rows are painted in reverse fuse order, so with
# nodata honoured the first valid source wins, as odc-stac's default fuser.
execute_chunk <- function(rows, gbox, resampling = "near", warp_args = character()) {
  r1 <- rows[1, ]
  cgt <- chunk_gt(gbox, r1$dst_xoff, r1$dst_yoff)
  acc <- matrix(0L, r1$dst_nrow, r1$dst_ncol)
  rows <- rows[order(rows$fuse_order), ]
  for (s in seq_len(nrow(rows))) {
    r <- rows[s, ]
    if (r$engine == "warp") {
      vrt <- sprintf("/vsimem/src_%s_%s.vrt", r$item, r$chunk)
      gdalraster::translate(r$href, vrt, c("-of", "VRT", "-srcwin", r$src_xoff, r$src_yoff,
                                           r$src_xsize, r$src_ysize), quiet = TRUE)
      dst <- gdalraster::create("MEM", "", r$dst_ncol, r$dst_nrow, 1, "UInt16", return_obj = TRUE)
      dst$setGeoTransform(cgt); dst$setProjection(gdalraster::srs_to_wkt(gbox$crs))
      dst$setNoDataValue(1, 0); dst$fillRaster(1, 0, 0)
      gdalraster::warp(vrt, dst, "", cl_arg = c("-r", resampling, "-srcnodata", r$nodata,
                                                 "-dstnodata", 0, warp_args), quiet = TRUE)
      m <- matrix(dst$read(1, 0, 0, r$dst_ncol, r$dst_nrow, r$dst_ncol, r$dst_nrow),
                  r$dst_nrow, byrow = TRUE)
      m[is.na(m)] <- 0L                # gdalraster reads nodata as NA
      dst$close(); gdalraster::vsi_unlink(vrt)
    } else if (r$engine == "bytes") {
      m <- read_bytes_window(r$href, r$src_xoff, r$src_yoff, r$src_xsize, r$src_ysize)
      m <- m[seq_len(r$dst_nrow), seq_len(r$dst_ncol), drop = FALSE]
    } else {
      ds <- new(gdalraster::GDALRaster, r$href)
      m <- matrix(ds$read(1, r$src_xoff, r$src_yoff, r$src_xsize, r$src_ysize,
                          r$src_xsize, r$src_ysize), r$src_ysize, byrow = TRUE)
      m[is.na(m)] <- 0L
      ds$close()
    }
    fill <- acc == 0 & m != 0          # first valid wins
    acc[fill] <- m[fill]
  }
  acc
}

#' Run a plan: returns an array [y, x, group] on the GeoBox.
execute_plan <- function(plan, gbox, resampling = "near", warp_args = character()) {
  groups <- sort(unique(plan$group))
  out <- array(0L, c(gbox$dim[2], gbox$dim[1], length(groups)),
               dimnames = list(NULL, NULL, as.character(groups)))
  for (g in seq_along(groups)) {
    pg <- plan[plan$group == groups[g], ]
    for (ck in unique(pg$chunk)) {
      rows <- pg[pg$chunk == ck, ]
      m <- execute_chunk(rows, gbox, resampling, warp_args)
      yy <- rows$dst_yoff[1] + seq_len(rows$dst_nrow[1])
      xx <- rows$dst_xoff[1] + seq_len(rows$dst_ncol[1])
      out[yy, xx, g] <- m
    }
  }
  out
}
