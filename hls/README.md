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
(int16, nodata -9999, 256 tiles, UTM north spelling, no `proj:*` on the
items):

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

The real run has not been done yet.
