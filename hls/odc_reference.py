"""odc-stac reference loads of HLS L30 B04 for the plan-table comparison.

Cases:
  warp       WGS84GRID30[49, 239] (south-up EPSG:4326, 5000 x 5000)
  aligned    1024 x 1024 on one granule's own grid (its native CRS, which
             for HLS is UTM north with negative northings in the south)
  aligned_s  the same window spelled in UTM south (EPSG:327xx), i.e. the
             same pixels with the false northing moved by 10,000,000 m

Reads hls/items.json (or $HLS_ITEMS). Writes ref_hls/*.tif, ref_hls/meta.json.
Earthdata access comes from the environment (GDAL_HTTP_* and ~/.netrc).
"""
import json, os, sys, time
import numpy as np
import pystac, rasterio
from affine import Affine
from odc.geo.geobox import GeoBox, GeoboxTiles
import odc.stac

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ITEMS = os.environ.get("HLS_ITEMS", os.path.join(ROOT, "hls", "items.json"))
OUT = os.environ.get("HLS_OUT", os.path.join(ROOT, "ref_hls"))
ASSET = "B04"
NODATA = -9999
os.makedirs(OUT, exist_ok=True)

fc = json.load(open(ITEMS))
for f in fc["features"]:   # local mock hrefs are relative to the repo root
    for a in f["assets"].values():
        if "://" not in a["href"] and not a["href"].startswith("/"):
            a["href"] = os.path.join(ROOT, a["href"])
items = [pystac.Item.from_dict(f) for f in fc["features"]]
coll = items[0].collection_id or "HLSL30.v2.0"
# HLS items carry no raster:bands: tell odc-stac the band's encoding
stac_cfg = {coll: {"assets": {ASSET: {"data_type": "int16", "nodata": NODATA}}},
            "*": {"assets": {ASSET: {"data_type": "int16", "nodata": NODATA}}}}

WGS84GRID30 = GeoboxTiles(
    GeoBox((600000, 1200000), Affine(0.0003, 0.0, -180.0, 0.0, 0.0003, -90.0), "epsg:4326"),
    (5000, 5000))
tile = WGS84GRID30[49, 239]

# aligned: the first granule whose id names an MGRS tile in zone 60
pick = next(i for i in items if ".T60" in i.id) if any(".T60" in i.id for i in items) else items[0]
with rasterio.open(pick.assets[ASSET].href) as d:
    tr, crs, (by, bx) = d.transform, d.crs, d.block_shapes[0]
aligned = GeoBox((1024, 1024), tr * Affine.translation(2 * bx, 2 * by), crs)
epsg = crs.to_epsg()
if epsg and 32601 <= epsg <= 32660:   # north zone: also spell it as the south zone
    s_crs = f"epsg:{epsg + 100}"
    a = aligned.affine
    aligned_s = GeoBox((1024, 1024), Affine(a.a, a.b, a.c, a.d, a.e, a.f + 10_000_000), s_crs)
else:
    aligned_s = None
print("aligned on", pick.id, crs, "block", (bx, by))

cases = [("warp", tile), ("aligned", aligned)] + ([("aligned_s", aligned_s)] if aligned_s else [])
meta = {"aligned_item": pick.id}
for name, gbox in cases:
    for rs in ["nearest", "bilinear"]:
        for ck_name, ck in [("nochunks", None), ("onechunk", {})]:
            if name != "warp" and ck_name == "onechunk":
                continue
            kw = {} if ck is None else {"chunks": ck}
            t0 = time.time()
            ds = odc.stac.load(items, bands=[ASSET], geobox=gbox, groupby="solar_day",
                               resampling=rs, stac_cfg=stac_cfg, **kw)
            arr = ds[ASSET].values
            dt = time.time() - t0
            key = f"{name}_{rs}" + ("" if ck_name == "nochunks" else "_onechunk")
            with rasterio.open(os.path.join(OUT, key + ".tif"), "w", driver="GTiff",
                               width=gbox.shape.x, height=gbox.shape.y, count=arr.shape[0],
                               dtype="int16", crs=str(gbox.crs), transform=gbox.affine,
                               nodata=NODATA, compress="DEFLATE", tiled=True) as d:
                d.write(arr.astype("int16"))
            meta[key] = {"times": [str(t)[:19] for t in ds.time.values], "shape": list(arr.shape),
                         "affine": list(gbox.affine)[:6], "crs": str(gbox.crs), "seconds": round(dt, 1),
                         "valid": [int((arr[i] != NODATA).sum()) for i in range(arr.shape[0])]}
            print(key, meta[key]["times"], meta[key]["valid"], f"{dt:.1f}s")

meta["versions"] = {"odc.stac": odc.stac.__version__, "gdal": rasterio.__gdal_version__}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
