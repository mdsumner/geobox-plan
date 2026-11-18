# The same comparison on real data: HLS L30

HLS L30 v2.0 (Landsat, harmonized onto Sentinel-2's MGRS tiles at 30 m)
over WGS84GRID30[49, 239]. Searching CMR-STAC needs no login; reading the
COGs needs an Earthdata login over HTTPS. There is no AWS account and no
egress charge.

## What is different from the fixtures

- Items carry no `proj:*` and no `raster:bands`. The R side reads each
  granule's grid and encoding from its header (`grid_from = "header"`);
  odc-stac is told the band's encoding through `stac_cfg`.
- B04 is int16 with nodata -9999.
- Southern-hemisphere granules are in UTM **north** zones with negative
  northings (EPSG:32660, origin y = -1,699,980 for T60LYH). A GeoBox in
  EPSG:32760 is the same grid moved by the 10,000,000 m false northing; the
  plan recognises that and still uses the byte engine (`aligned_s`).
- One Landsat pass over the tile is several granules (MGRS tiles in zones
  60 and 1), so each group has more members than in the fixtures.

## Run (on ecostats)

```sh
export GDAL_HTTP_COOKIEFILE=~/.urs_cookies
export GDAL_HTTP_COOKIEJAR=~/.urs_cookies
export GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR
export CPL_VSIL_CURL_ALLOWED_EXTENSIONS=.tif
# ~/.netrc: machine urs.earthdata.nasa.gov login ... password ...

python3 hls/search.py          # -> hls/items.json
python3 hls/odc_reference.py   # -> ref_hls/*.tif, ref_hls/meta.json
Rscript hls/run_hls.R          # -> ref_hls/plan_*.csv, ref_hls/comparison.csv
```

`HLS_DATES` changes the search window (default `2024-12-27/2024-12-30`).

## Cases

| case | GeoBox | expected engine |
|---|---|---|
| `warp` | WGS84GRID30[49, 239], EPSG:4326, south-up, 5000 x 5000 | warp |
| `aligned` | 1024 x 1024 on a zone 60 granule's own grid, 2 blocks in | bytes |
| `aligned_s` | the same window in EPSG:32760 | bytes |

`odc_reference.py` loads `warp` twice: without chunks, and with
`chunks={}` (one 5000 x 5000 chunk), which on the fixtures lost the most
pixels near 180.

## Dry run without Earthdata

`hls/make_mock.py` turns fixture scenes C and B into HLS-shaped granules
(int16, nodata -9999, DEFLATE with predictor 2, 256 tiles, UTM north
spelling, no `proj:*` on the items):

```sh
python3 fixtures/make_fixtures.py
python3 hls/make_mock.py
HLS_ITEMS=hls/mock/items.json HLS_OUT=ref_hls_mock python3 hls/odc_reference.py
HLS_ITEMS=hls/mock/items.json HLS_OUT=ref_hls_mock Rscript hls/run_hls.R
```

On the mock (GDAL 3.8.4 for R, 3.12.2 under odc-stac):

- `warp`: the R plan fills every pixel gdalwarp fills, plus 26 on one
  scene edge on the partial day (0 with `-et 0` on both sides). odc-stac
  loses about 1.02 M pixels per day without chunks and 2.65 M / 5.97 M
  with one chunk, all near 180.
- `aligned` and `aligned_s`: one `bytes` row each; identical to odc-stac
  and to gdalwarp. The day with no granule in the window comes back as an
  all-nodata layer on every side.

With `HLS_WARP_ARGS="-et 0"` (exact transformer on both sides) the `warp`
case is identical to gdalwarp for nearest, and the byte engine decodes the
predictor-2 tiles to the same values GDAL reads.

## First real run (2026-11-18)

Nine granules (five on 2024-12-28 22:06Z, four on 2024-12-29 22:00Z) in
UTM 60N and 1N. All are DEFLATE with predictor 2 and 256 x 256 tiles.

What it showed, before the fixes in the patch after it:

- **The byte engine was not used.** Predictor 2 was not decoded then, so the
  aligned cases fell back to `window` reads. These were identical to
  odc-stac and gdalwarp (1,030,298 valid pixels).
- **One pass was split across two days.** T01LAC on 2024-12-29 centres at
  -179.96 E while the rest of its pass centres at +179.7 to +179.98, so a
  per-item solar day put it a day early. Groups are now formed per pass
  (`pass_groups()`): items whose start times fall within 10 minutes are one
  pass, labelled by the solar day of the pass's circular-mean longitude.
  odc-stac kept that granule with its pass.
- **Coverage, R plan vs whole gdalwarp:** 15,191,348 vs 15,191,386 on day
  1, and 1,199,655 vs 1,199,671 on day 2, with the default transformer.
- **Values with the default transformer:** they agree on only 96.9%
  (nearest) and 28% (bilinear) of pixels. The smooth fixture field hid
  this: a neighbouring pixel there often has the same integer value, so the
  fixtures understated how often chunked and whole-tile warps pick
  differently. With `-et 0` the two should agree, as they do on the mock;
  this is to be confirmed on the real data.
- **odc-stac depends on chunking:** 14,453,383 valid pixels on day 1 with no
  chunks against 13,883,499 with one chunk, and 1,386,711 against 934,763 on
  day 2.
- **CMR-STAC bboxes for granules on 180 are degenerate**, for example
  `[180, -16.95, 179.9998, -16.25]`. `population()` takes its centre from the
  footprint vertices, so it is not affected.
