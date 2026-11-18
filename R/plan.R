# A GeoBox plan table: population x target chunks -> one row per read.
#
# The plan is a data frame. Every row says: to fill this target chunk, for
# this group (time step), read this window of this asset, with this engine.
# Building it reads no pixels (at most a file header, when an item carries no
# proj:* metadata). Executing it is a loop over rows.
#
# Dependencies: gdalraster (transforms, warp, reads, VSIFile), jsonlite.
# Grid arithmetic is written out here; in a package it is vaster's.
#
# Inside execution nodata is NA, whatever the source's nodata value is.

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

# Centre longitude as a circular mean of the footprint's vertices, so a
# granule at the antimeridian centres near 180, not near 0. A bbox that
# claims [-180, 180] for such a granule would otherwise put it at 0 and
# shift its solar day.
lon_centre <- function(f) {
  lon <- NULL
  if (!is.null(f$geometry)) lon <- unlist(f$geometry$coordinates)   # x, y, x, y, ...
  if (!is.null(lon) && length(lon) >= 2) {
    lon <- lon[seq(1, length(lon), by = 2)]
  } else {
    bb <- unlist(f$bbox); lon <- c(bb[1], bb[3])
  }
  r <- lon * pi / 180
  atan2(mean(sin(r)), mean(cos(r))) * 180 / pi
}

# "EPSG:nnnn" when the WKT is a known code, else the WKT itself.
short_crs <- function(wkt) {
  e <- tryCatch(gdalraster::srs_find_epsg(wkt), error = function(e) NA)
  if (length(e) && !is.na(e[1]) && nzchar(e[1])) as.character(e[1]) else wkt
}

# Grid and encoding from the file header: used when an item has no proj:*.
# Costs one header read per asset.
header_grid <- function(href) {
  ds <- new(gdalraster::GDALRaster, href); on.exit(ds$close())
  ims <- ds$getMetadata(0, "IMAGE_STRUCTURE")
  val <- function(k) { x <- sub(paste0("^", k, "="), "", grep(paste0("^", k, "="), ims, value = TRUE)); if (length(x)) x else NA_character_ }
  list(crs = short_crs(ds$getProjection()), gt = ds$getGeoTransform(),
       ncol = ds$getRasterXSize(), nrow = ds$getRasterYSize(),
       nodata = ds$getNoDataValue(1), dtype = ds$getDataTypeName(1),
       block = ds$getBlockSize(1), compression = val("COMPRESSION"), predictor = val("PREDICTOR"))
}

#' One row per (item, asset) from STAC items. Grid from proj:* (item or asset
#' level, proj:epsg or proj:code) when present, else from the file header.
#' Encoding (dtype, block size, compression) always from the header when
#' header = TRUE, since the byte engine needs it.
population <- function(features, asset = "red", header = TRUE, href_fn = identity) {
  rows <- lapply(features, function(f) {
    p <- f$properties; a <- f$assets[[asset]]
    href <- href_fn(a$href)
    pget <- function(k) if (!is.null(a[[k]])) a[[k]] else p[[k]]
    crs <- if (!is.null(pget("proj:code"))) pget("proj:code") else if (!is.null(pget("proj:epsg"))) paste0("EPSG:", pget("proj:epsg")) else NA
    tr <- unlist(pget("proj:transform")); shp <- unlist(pget("proj:shape"))
    h <- if (header || is.na(crs) || is.null(tr) || is.null(shp)) header_grid(href) else NULL
    if (!is.null(tr) && !is.null(shp) && !is.na(crs)) {
      gt <- c(tr[3], tr[1], tr[2], tr[6], tr[4], tr[5]); ncol <- shp[2]; nrow <- shp[1]
      grid_from <- "proj"
    } else {
      crs <- h$crs; gt <- h$gt; ncol <- h$ncol; nrow <- h$nrow; grid_from <- "header"
    }
    nd <- a[["raster:bands"]][[1]]$nodata
    if (is.null(nd) && !is.null(h)) nd <- h$nodata
    dt <- if (!is.null(p$datetime)) p$datetime else p$start_datetime
    data.frame(item = f$id, datetime = as.POSIXct(sub("Z$", "", dt), "%Y-%m-%dT%H:%M:%OS", tz = "UTC"),
               asset = asset, href = href, src_crs = crs, grid_from = grid_from,
               sgt1 = gt[1], sgt2 = gt[2], sgt3 = gt[3], sgt4 = gt[4], sgt5 = gt[5], sgt6 = gt[6],
               src_ncol = as.integer(ncol), src_nrow = as.integer(nrow),
               nodata = if (is.null(nd)) NA_real_ else as.numeric(nd),
               dtype = if (is.null(h)) NA_character_ else h$dtype,
               block_x = if (is.null(h)) NA_integer_ else as.integer(h$block[1]),
               block_y = if (is.null(h)) NA_integer_ else as.integer(h$block[2]),
               compression = if (is.null(h)) NA_character_ else h$compression,
               predictor = if (is.null(h)) NA_character_ else h$predictor,
               lon_centre = lon_centre(f))
  })
  do.call(rbind, rows)
}

#' odc-stac's "solar_day": UTC time shifted by longitude / 15 hours.
solar_day <- function(datetime, lon) as.Date(datetime + lon / 15 * 3600, tz = "UTC")

#' Group a population by acquisition pass, labelled by solar day.
#'
#' Per-item solar_day splits a pass that straddles the antimeridian: a
#' granule centred at -179.96 shifts by -12 h and one at +179.98 by +12 h, so
#' one overpass lands on two days. Items are first joined into passes (start
#' times within `tol` seconds of the previous item, in time order), and each
#' pass gets one solar day from the circular mean longitude of its members.
pass_groups <- function(datetime, lon, tol = 600) {
  o <- order(datetime)
  t <- as.numeric(datetime[o])
  pass <- integer(length(t)); pass[o] <- cumsum(c(TRUE, diff(t) > tol))
  r <- lon * pi / 180
  plon <- tapply(seq_along(lon), pass, function(i) atan2(mean(sin(r[i])), mean(cos(r[i]))) * 180 / pi)
  ptime <- tapply(seq_along(lon), pass, function(i) min(datetime[i]))
  ptime <- as.POSIXct(ptime, origin = "1970-01-01", tz = "UTC")
  solar_day(ptime, plon)[as.character(pass)]
}

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

# Same grid, possibly in a different CRS spelling: if the two CRSs are the
# same projection up to a false easting/northing (UTM 60N with negative
# northings, as HLS uses, against UTM 60S), the target origin moved into the
# source CRS lands on a source pixel corner. The test is done by moving the
# chunk's corners, not by comparing CRS strings.
same_grid_offset <- function(dst_crs, cgt, ncol, nrow, src, tol = 1e-6) {
  sgt <- c(src$sgt1, src$sgt2, src$sgt3, src$sgt4, src$sgt5, src$sgt6)
  if (any(abs(cgt[c(3, 5)]) > 0) || any(abs(sgt[c(3, 5)]) > 0)) return(NULL)
  same <- gdalraster::srs_is_same(dst_crs, src$src_crs)
  if (!same) {
    # corners and an interior point moved into the source CRS must differ
    # from the target by one constant shift, i.e. a pure offset
    xy <- rbind(c(cgt[1], cgt[4]), c(cgt[1] + ncol * cgt[2], cgt[4] + nrow * cgt[6]),
                c(cgt[1] + ncol * cgt[2] / 2, cgt[4]))
    s <- suppressWarnings(gdalraster::transform_xy(xy, dst_crs, src$src_crs))
    if (any(!is.finite(s))) return(NULL)
    d <- s - xy
    if (max(abs(sweep(d, 2, d[1, ]))) > 1e-3) return(NULL)
    o <- c(cgt[1] + d[1, 1], cgt[4] + d[1, 2])
  } else {
    o <- c(cgt[1], cgt[4])
  }
  if (abs(cgt[2] - sgt[2]) > 1e-9 || abs(cgt[6] - sgt[6]) > 1e-9) return(NULL)
  ox <- (o[1] - sgt[1]) / sgt[2]; oy <- (o[2] - sgt[4]) / sgt[6]
  if (abs(ox - round(ox)) > tol || abs(oy - round(oy)) > tol) return(NULL)
  c(round(ox), round(oy))
}

# "bytes":  same grid (see above), window on whole source blocks, and an
#           encoding the R decoder handles (DEFLATE, predictor 1 or 2, 16-bit).
# "window": same grid otherwise: a plain windowed read, no resampling.
# "warp":   everything else.
choose_engine <- function(dst_crs, cgt, ncol, nrow, src) {
  o <- same_grid_offset(dst_crs, cgt, ncol, nrow, src)
  if (is.null(o)) return(list(engine = "warp"))
  bx <- src$block_x; by <- src$block_y
  decodable <- isTRUE(src$compression == "DEFLATE") && (is.na(src$predictor) || src$predictor %in% c("1", "2")) &&
    isTRUE(src$dtype %in% c("UInt16", "Int16"))
  onblocks <- !is.na(bx) && o[1] %% bx == 0 && o[2] %% by == 0 &&
    (ncol %% bx == 0 || o[1] + ncol >= src$src_ncol) &&
    (nrow %% by == 0 || o[2] + nrow >= src$src_nrow)
  list(engine = if (onblocks && decodable) "bytes" else "window", ox = o[1], oy = o[2])
}

# ---- the plan ---------------------------------------------------------------

#' Build the plan table. pad is the resampling kernel's reach in source
#' pixels (1 for nearest, 2 for bilinear/cubic).
geobox_plan <- function(pop, gbox, chunk = c(2048L, 2048L), pad = 2L) {
  ch <- geobox_chunks(gbox, chunk)
  pop$group <- pass_groups(pop$datetime, pop$lon_centre)
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
      e <- choose_engine(gbox$crs, cgt, ch$ncol[k], ch$nrow[k], src)
      if (e$engine != "warp") {          # exact window, no pad needed
        x0 <- max(0, e$ox); y0 <- max(0, e$oy)
        x1 <- min(src$src_ncol, e$ox + ch$ncol[k]); y1 <- min(src$src_nrow, e$oy + ch$nrow[k])
        if (x1 <= x0 || y1 <= y0) next
        w <- c(x0, y0, x1 - x0, y1 - y0)
        dx <- x0 - e$ox; dy <- y0 - e$oy  # where the window lands in the chunk
      } else { dx <- NA; dy <- NA }
      out[[length(out) + 1L]] <- data.frame(
        chunk = ch$chunk[k], dst_xoff = ch$xoff[k], dst_yoff = ch$yoff[k],
        dst_ncol = ch$ncol[k], dst_nrow = ch$nrow[k],
        group = src$group, fuse_order = src$fuse_order, item = src$item,
        asset = src$asset, href = src$href, src_crs = src$src_crs,
        engine = e$engine, src_xoff = w[1], src_yoff = w[2], src_xsize = w[3], src_ysize = w[4],
        in_chunk_x = dx, in_chunk_y = dy,
        nodata = src$nodata, dtype = src$dtype, block_x = src$block_x, block_y = src$block_y,
        compression = src$compression, predictor = src$predictor)
    }
  }
  p <- do.call(rbind, out)
  p <- p[order(p$group, p$chunk, p$fuse_order), ]
  # every group in the population, including ones with no rows here, so an
  # empty time step comes back as all nodata (as odc-stac returns it)
  attr(p, "groups") <- sort(unique(pop$group))
  p
}

# ---- execution --------------------------------------------------------------

# Byte engine: block offsets and sizes come from the TIFF tags (GDAL exposes
# them as BLOCK_OFFSET_x_y / BLOCK_SIZE_x_y). Here they are read from GDAL
# metadata for brevity; in the real thing they are a byte-ref table
# (tabularberg refs_* rows, a kerchunk Parquet, a get-refs harvest) and GDAL
# is not involved at read time. The bytes are fetched with VSIFile, so local
# paths and /vsicurl/ both work. DEFLATE, predictor none or 2 (horizontal
# differencing, which HLS uses), 16-bit integers.
read_bytes_window <- function(r) {
  bx <- r$block_x; by <- r$block_y
  signed <- identical(r$dtype, "Int16")
  ds <- new(gdalraster::GDALRaster, r$href); on.exit(ds$close())
  f <- new(gdalraster::VSIFile, r$href, "r"); on.exit(f$close(), add = TRUE)
  out <- matrix(NA_real_, r$src_ysize, r$src_xsize)
  for (bj in seq(r$src_yoff %/% by, (r$src_yoff + r$src_ysize - 1) %/% by)) {
    for (bi in seq(r$src_xoff %/% bx, (r$src_xoff + r$src_xsize - 1) %/% bx)) {
      off <- as.numeric(ds$getMetadataItem(1, sprintf("BLOCK_OFFSET_%d_%d", bi, bj), "TIFF"))
      sz <- as.numeric(ds$getMetadataItem(1, sprintf("BLOCK_SIZE_%d_%d", bi, bj), "TIFF"))
      if (!is.finite(off) || !is.finite(sz) || sz == 0) next   # sparse block: nodata
      f$seek(off, gdalraster::SEEK_SET); raw <- f$read(sz)
      v <- readBin(memDecompress(raw, type = "gzip"), "integer", n = bx * by,
                   size = 2, signed = FALSE, endian = "little")
      m <- matrix(as.numeric(v), by, bx, byrow = TRUE)
      if (identical(r$predictor, "2")) {   # horizontal differencing, per row, mod 2^16
        m <- t(apply(m, 1, cumsum)) %% 65536
      }
      if (signed) m[m > 32767] <- m[m > 32767] - 65536
      rr <- (bj * by - r$src_yoff) + seq_len(by); cc <- (bi * bx - r$src_xoff) + seq_len(bx)
      kr <- rr >= 1 & rr <= r$src_ysize; kc <- cc >= 1 & cc <= r$src_xsize
      out[rr[kr], cc[kc]] <- m[kr, kc]
    }
  }
  if (!is.na(r$nodata)) out[out == r$nodata] <- NA
  out
}

# Warp engine: the plan's source window becomes a VRT (gdal_translate
# -srcwin), so the plan, not the warper's own boundary walk, decides which
# source pixels are in play. Rows are taken in fuse order and the first valid
# value wins, as odc-stac's default fuser.
execute_chunk <- function(rows, gbox, resampling = "near", warp_args = character()) {
  r1 <- rows[1, ]
  cgt <- chunk_gt(gbox, r1$dst_xoff, r1$dst_yoff)
  acc <- matrix(NA_real_, r1$dst_nrow, r1$dst_ncol)
  rows <- rows[order(rows$fuse_order), ]
  for (s in seq_len(nrow(rows))) {
    r <- rows[s, ]
    dtype <- if (is.na(r$dtype)) "Float64" else r$dtype
    if (r$engine == "warp") {
      vrt <- sprintf("/vsimem/src_%s_%s.vrt", r$item, r$chunk)
      gdalraster::translate(r$href, vrt, c("-of", "VRT", "-srcwin", r$src_xoff, r$src_yoff,
                                           r$src_xsize, r$src_ysize), quiet = TRUE)
      nd <- if (is.na(r$nodata)) 0 else r$nodata
      dst <- gdalraster::create("MEM", "", r$dst_ncol, r$dst_nrow, 1, dtype, return_obj = TRUE)
      dst$setGeoTransform(cgt); dst$setProjection(gdalraster::srs_to_wkt(gbox$crs))
      dst$setNoDataValue(1, nd); dst$fillRaster(1, nd, 0)
      gdalraster::warp(vrt, dst, "", cl_arg = c("-r", resampling, "-srcnodata", nd,
                                                 "-dstnodata", nd, warp_args), quiet = TRUE)
      m <- matrix(dst$read(1, 0, 0, r$dst_ncol, r$dst_nrow, r$dst_ncol, r$dst_nrow),
                  r$dst_nrow, byrow = TRUE)            # gdalraster reads nodata as NA
      dst$close(); gdalraster::vsi_unlink(vrt)
    } else {
      if (r$engine == "bytes") {
        w <- read_bytes_window(r)
      } else {
        ds <- new(gdalraster::GDALRaster, r$href)
        w <- matrix(ds$read(1, r$src_xoff, r$src_yoff, r$src_xsize, r$src_ysize,
                            r$src_xsize, r$src_ysize), r$src_ysize, byrow = TRUE)
        ds$close()
      }
      m <- matrix(NA_real_, r$dst_nrow, r$dst_ncol)
      m[r$in_chunk_y + seq_len(nrow(w)), r$in_chunk_x + seq_len(ncol(w))] <- w
    }
    fill <- is.na(acc) & !is.na(m)    # first valid wins
    acc[fill] <- m[fill]
  }
  acc
}

#' Run a plan: returns an array [y, x, group] on the GeoBox, NA for nodata.
execute_plan <- function(plan, gbox, resampling = "near", warp_args = character()) {
  groups <- attr(plan, "groups")
  if (is.null(groups)) groups <- sort(unique(plan$group))
  out <- array(NA_real_, c(gbox$dim[2], gbox$dim[1], length(groups)),
               dimnames = list(NULL, NULL, as.character(groups)))
  for (g in seq_along(groups)) {
    pg <- plan[plan$group == groups[g], ]
    if (!nrow(pg)) next
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

# ---- shared helpers for the run scripts -------------------------------------

#' Read a multi-band file as [y, x, band], NA for nodata.
read_stack <- function(f) {
  ds <- new(gdalraster::GDALRaster, f); on.exit(ds$close())
  nx <- ds$getRasterXSize(); ny <- ds$getRasterYSize(); nb <- ds$getRasterCount()
  simplify2array(lapply(seq_len(nb), function(b) matrix(ds$read(b, 0, 0, nx, ny, nx, ny), ny, byrow = TRUE)))
}

#' GDAL's own answer: one gdalwarp per group over the whole GeoBox, sources
#' in reverse fuse order (gdalwarp paints later sources over earlier ones).
gdal_whole <- function(pop, gbox, resampling, warp_args = character()) {
  pop$group <- pass_groups(pop$datetime, pop$lon_centre)
  gs <- sort(unique(pop$group))
  dtype <- if (is.na(pop$dtype[1])) "Float64" else pop$dtype[1]
  nd <- if (is.na(pop$nodata[1])) 0 else pop$nodata[1]
  simplify2array(lapply(gs, function(g) {
    s <- rev(pop$href[pop$group == g])
    dst <- gdalraster::create("MEM", "", gbox$dim[1], gbox$dim[2], 1, dtype, return_obj = TRUE)
    dst$setGeoTransform(gbox$gt); dst$setProjection(gdalraster::srs_to_wkt(gbox$crs))
    dst$setNoDataValue(1, nd); dst$fillRaster(1, nd, 0)
    gdalraster::warp(s, dst, "", cl_arg = c("-r", resampling, "-srcnodata", nd, "-dstnodata", nd, warp_args), quiet = TRUE)
    m <- matrix(dst$read(1, 0, 0, gbox$dim[1], gbox$dim[2], gbox$dim[1], gbox$dim[2]),
                gbox$dim[2], byrow = TRUE)
    dst$close(); m
  }))
}

compare <- function(a, b) {
  va <- !is.na(a); vb <- !is.na(b); both <- va & vb
  d <- abs(a[both] - b[both])
  data.frame(valid_a = sum(va), valid_b = sum(vb), only_a = sum(va & !vb), only_b = sum(vb & !va),
             equal = mean(d == 0), max_abs_diff = if (length(d)) max(d) else NA)
}
