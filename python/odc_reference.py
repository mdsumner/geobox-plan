"""odc-stac reference loads for the plan-table comparison.

Case 1 (warp): WGS84GRID30[49, 239], the GeoBox from the Dec 2025 snippet.
Case 2 (aligned): 1024 x 1024 on scene C's own UTM grid, on its 512 blocks.

Writes ref/*.tif (one band per time step) and ref/meta.json.
"""
import json, os, sys
import numpy as np
import pystac, rasterio
from affine import Affine
from odc.geo.geobox import GeoBox, GeoboxTiles
import odc.stac

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "ref")
os.makedirs(OUT, exist_ok=True)

fc = json.load(open(os.path.join(ROOT, "fixtures", "items.json")))
items = [pystac.Item.from_dict(f) for f in fc["features"]]

WGS84GRID30 = GeoboxTiles(
    GeoBox((600000, 1200000), Affine(0.0003, 0.0, -180.0, 0.0, 0.0003, -90.0), "epsg:4326"),
    (5000, 5000))
tile = WGS84GRID30[49, 239]

c = fc["features"][2]["properties"]
tc = Affine(*c["proj:transform"])
aligned = GeoBox((1024, 1024), tc * Affine.translation(2 * 512, 4 * 512), f"epsg:{c['proj:epsg']}")

meta = {}
for name, gbox in [("warp", tile), ("aligned", aligned)]:
    for rs in ["nearest", "bilinear"]:
        ds = odc.stac.load(items, bands=["red"], geobox=gbox, groupby="solar_day",
                           resampling=rs, chunks={})
        ds = ds.compute()
        arr = ds.red.values  # (time, y, x)
        path = os.path.join(OUT, f"{name}_{rs}.tif")
        a = gbox.affine
        with rasterio.open(path, "w", driver="GTiff", width=gbox.shape.x, height=gbox.shape.y,
                           count=arr.shape[0], dtype=arr.dtype, crs=str(gbox.crs),
                           transform=a, nodata=0, compress="DEFLATE", tiled=True) as d:
            d.write(arr)
        meta[f"{name}_{rs}"] = {
            "times": [str(t)[:10] for t in ds.time.values],
            "shape": list(arr.shape), "affine": list(a)[:6], "crs": str(gbox.crs),
            "valid": [int((arr[i] != 0).sum()) for i in range(arr.shape[0])],
        }
        print(name, rs, meta[f"{name}_{rs}"]["times"], arr.shape, meta[f"{name}_{rs}"]["valid"])

meta["versions"] = {"odc.stac": odc.stac.__version__, "gdal": rasterio.__gdal_version__}
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
