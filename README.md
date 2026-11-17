# GeoBox plan tables: odc-stac's load, as a table

Design note, 2026-11-18. Working material for the collection -> GeoBox ->
cube layer in R (zaro, ndr, tabularberg, rwarp, with rgdal-dev/adgl as one
engine). Every number below comes from the scripts in this directory.

## What was done

`odc.stac.load(items, geobox = tile, groupby = "solar_day")` was run side by
side with an R plan table on the same GeoBoxes, and both were checked
against a single whole-tile `gdalwarp`.

The target is the tile from the December 2025 reticulate snippet:

```python
WGS84GRID30 = GeoboxTiles(GeoBox((600000, 1200000),
    Affine(0.0003, 0, -180, 0, 0.0003, -90), "epsg:4326"), (5000, 5000))
tile = WGS84GRID30[49, 239]   # lon 178.5..180, lat -16.5..-15, south-up
```

The workspace this was built in cannot reach Earth Search or the Landsat
buckets, so the sources are synthetic. There are three 30 m Landsat-like
scenes (uint16, nodata 0, DEFLATE, 512 x 512 tiles, COG layout), whose values
are a smooth function of lon/lat, so any engine that samples the right place
gets the right number:

| scene | date (UTC) | CRS | lon span | note |
|---|---|---|---|---|
| A | 2024-12-28 | EPSG:32760 | 178.3 .. 179.6 | |
| B | 2024-12-28 | EPSG:32701 | 179.3 .. 180.2 | crosses 180; values +500 so paint order shows |
| C | 2024-12-29 | EPSG:32760 | 178.4 .. 180.1 | crosses 180 |

There are two cases:

- **warp:** the tile above.
- **aligned:** a 1024 x 1024 GeoBox on scene C's own UTM grid, sitting on its
  512 blocks. This is the case where no resampling is needed.

## The plan table

The plan has one row per read: to fill this target chunk for this group, read
this window of this asset with this engine. Building it reads no pixels.

| column | meaning |
|---|---|
| `chunk`, `dst_xoff`, `dst_yoff`, `dst_ncol`, `dst_nrow` | target chunk, in pixels of the GeoBox |
| `group` | time step (here `solar_day`) |
| `fuse_order` | order within the group; the first valid value wins |
| `item`, `asset`, `href`, `src_crs`, `nodata` | the source |
| `engine` | `bytes`, `window` or `warp` |
| `src_xoff`, `src_yoff`, `src_xsize`, `src_ysize` | source window, in source pixels |

The warp case gives 21 rows (9 chunks of 2048 x 2048, against 3 scenes, in 2
groups) and takes 0.13 s to build. Scene A has no row for chunk `0_2`, because
it stops at 179.6 E and that chunk starts at 179.73 E. The aligned case gives
2 rows: C is `bytes` and A is `window` (same CRS and pixel size, integer
offset, not on A's blocks).

The plan is built in four steps:

1. **Chunk the GeoBox.** Chunk origins come from integer pixel offsets of the
   parent (vaster arithmetic).
2. **Group the population.** `solar_day` is the UTC time shifted by
   longitude / 15 hours, using an antimeridian-aware centre longitude.
3. **Find each source window by projecting the target into the source.**
   Each chunk's pixel corners (dense edges plus a 17 x 17 interior lattice)
   are transformed into each source CRS, then turned into source pixels and
   padded by the kernel's reach (2 px). The direction matters (finding 1). A
   UTM scene that straddles 180 is continuous in its own CRS, and its outline
   only breaks once it is projected into longitude.
4. **Choose the engine.** `bytes` if the source has the same CRS and pixel
   size, sits at an integer offset, and the window falls on whole source
   blocks. `window` if the grid matches but the blocks don't. `warp`
   otherwise.

Each plan is executed one chunk at a time:

- `warp` rows become `gdal_translate -srcwin` VRTs, so the plan, not the
  warper's own boundary walk, decides which source pixels are in play. They
  are then warped into an in-memory chunk with the GeoBox's own geotransform,
  which can be south-up.
- `bytes` rows read the tile bytes and inflate them in R with
  `memDecompress()`. Here the offsets come from GDAL's
  `BLOCK_OFFSET_x_y` / `BLOCK_SIZE_x_y` metadata; in the real thing they
  would be a byte-ref table such as tabularberg `refs_*`, a kerchunk Parquet
  or a get-refs harvest, and GDAL would not be involved at read time.
- Within a group the first valid value wins, which is odc-stac's default
  fuser.

## Results

Valid pixels per time step (warp case, 25,000,000 per step):

| group | R plan | whole gdalwarp | odc-stac (`chunks={}`) | odc-stac missing |
|---|---|---|---|---|
| day 1 (A + B) | 25,000,000 | 25,000,000 | 22,350,000 | 2,650,000, east of 179.841 E |
| day 2 (C) | 25,000,000 | 25,000,000 | 19,035,000 | 5,965,000, east of 179.642 E |

These counts are the same for nearest and bilinear.

Where both have data, the R plan and odc-stac agree on 99.90% of pixels
(nearest, day 1). Almost all of the rest differ by 1, a neighbouring source
pixel; a few on scene A's edge differ by about 500. The comparison against a
whole-tile gdalwarp has the same pattern, which finding 4 explains.

In the aligned case all three agree exactly, for both groups and both
resamplings. That includes the tile bytes decoded in R.

With the exact transformer (`-et 0`) on both sides, the chunked R plan is
identical to a single whole-tile gdalwarp for nearest (0 differing pixels on
both days). For bilinear it differs by integer rounding only, at most 1 count
on 0.06% (day 1) and 0.01% (day 2) of pixels.

Versions: odc-stac 0.5.3 with rasterio's GDAL 3.12.2. The R side used
gdalraster 2.7.0 with GDAL 3.8.4.

## Findings

### 1. odc-stac loses pixels near the antimeridian, depending on chunk size

The scenes that lose data are the ones whose footprint crosses 180 (B and C).
Scene A stops at 179.6 and comes through exactly. Moving the same GeoBox to
178.2..179.7 E gives a correct result. How much is lost depends on the dask
chunking:

| `chunks` | day 2 valid | east edge |
|---|---|---|
| none, or 2048 x 2048 | 23,978,248 | 179.940 E |
| `{}` or 5000 x 5000 (one chunk) | 19,035,000 | 179.642 E |

`python/odc_antimeridian.py` is a minimal repro on scene C alone, and it
needs nothing but the fixtures. The cause is not confirmed: odc-geo's source
code was not read. The pattern fits a source footprint or ROI computed by
projecting the source into longitude, where it wraps at 180, with coarser
sampling of a larger chunk giving a worse cut. A whole-tile gdalwarp on the
same inputs is complete.

The design rule for the R side is to project the target into the source, and
never the source into the target. The plan builder does this, and it covers
all 25 M pixels on both days.

### 2. One plan, three engines

The aligned case shows the opportunity that odc-stac doesn't take. When the
target grid lands on a source's own grid, there is nothing to warp: the chunk
is the source's bytes. odc-stac goes through the warper for every source. The
plan table makes the engine a column, so each (target chunk, source) pair gets
the cheapest correct read:

- `bytes`: a byte-ref table plus decode (tabularberg / zaro).
- `window`: a plain windowed read (GDAL7 / adgl `collect()`).
- `warp`: a GDAL warp (adgl `warp()`), and later rwarp.

### 3. The plan should carry the warp's error threshold

GDAL's approximate transformer (default `-et 0.125`) interpolates along each
output chunk's rows. A chunked warp therefore rounds differently from a
whole-tile warp, by one source pixel at about 0.14% of pixels here. With
`-et 0` the chunked plan is identical to the whole-tile warp. If plans are
meant to give the same answer whatever the chunking, the transformer settings
belong in the plan, either as a column or as plan-level metadata.

### 4. Group labels differ from odc-stac

Both group the same items together. odc-stac labels each group with the UTC
date of its first item (2024-12-28, 2024-12-29). The R plan labels it with the
solar date (2024-12-29, 2024-12-30), which at 179 E is the following day.
Labels are a choice to make deliberately, and to say in the output's
metadata.

### 5. Tile origins are floating point

`WGS84GRID30[49, 239]` has origin `178.49999999999994`, which is
`-180 + 1195000 * 0.0003`. Computing chunk origins from integer offsets of
the parent, rather than accumulating a step, keeps adjacent chunks exactly
abutting, and the R side should do the same so that its grids are equal to
odc-geo's.

## Where this would live

| odc-stac | R piece |
|---|---|
| GeoBox, GeoboxTiles | vaster (grid arithmetic), grout (tiles) |
| item population | starc rows; here `population()` from STAC JSON with `proj:*` |
| `groupby`, fuser | columns of the plan table, ordinary verbs |
| per-source warped read | adgl `warp()` on GDAL7; rwarp later |
| byte fetch (absent in odc-stac) | tabularberg `refs_*`, zaro |
| dask array | ndr / altarr; targets branching over plan rows |

adgl stays the single-source GDAL verb layer and supplies the `window` and
`warp` engines. adgl's `report()` carries straight over: a STAC item with no
`proj:*`, a `proj:shape` that disagrees with the asset, or a nodata missing
from `raster:bands` is a finding that changes the plan. The plan builder sits
above adgl, either in ndr or in a small package of its own, and that choice is
what settles the name.

## Next steps

1. **Run the real Landsat item.** Use LC08_L2SR_074071_20241228_02_T1 on the
   same tile, on Setonix or the Nectar VM, where Earth Search and the
   requester-pays bucket are reachable. Use `population()` on the fetched
   item JSON with hrefs rewritten to `/vsis3/` and `AWS_REQUEST_PAYER=requester`.
   This hasn't been run.
2. **Write plan rows to Parquet** and read them back as the targets branching
   input.
3. **Swap the byte engine's source** from GDAL metadata to a tabularberg
   `refs_*` table, so that reads need no GDAL at all.
4. **Add a `src_ovr` column**, so that a coarse GeoBox picks an overview level
   per source as adgl's `query(dim =)` does.
5. **Add rwarp as a fourth engine value**, and compare it against the `warp`
   rows chunk by chunk.

## Files

| path | what |
|---|---|
| `fixtures/make_fixtures.py` | builds `fixtures/data/{A,B,C}.tif` and `fixtures/items.json` |
| `python/odc_reference.py` | odc-stac loads for both cases -> `ref/*.tif`, `ref/meta.json` |
| `python/odc_antimeridian.py` | minimal repro for finding 1 |
| `R/plan.R` | `geobox()`, `population()`, `geobox_plan()`, `execute_plan()` |
| `R/run.R` | builds and runs both plans, compares -> `ref/comparison.csv`, `ref/plan_*.csv` |
| `R/check_exact.R` | chunked plan vs whole gdalwarp with `-et 0` |

Order: `python3 fixtures/make_fixtures.py`, `python3 python/odc_reference.py`,
`Rscript R/run.R`, `Rscript R/check_exact.R`. Needs odc-stac, rasterio and
pystac for Python, and gdalraster and jsonlite for R.
