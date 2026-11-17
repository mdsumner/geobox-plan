"""Minimal repro: odc-stac loses pixels near 180 for a UTM scene that crosses it.

Scene C (EPSG:32760, lon 178.4..180.1) onto a 5000 x 5000 0.0003 deg GeoBox
at lon 178.5..180. gdalwarp fills every pixel; odc.stac.load does not, and
how much it loses depends on the dask chunk size. The same GeoBox moved to
lon 178.2..179.7 (away from 180) is filled correctly.
"""
import json, os, numpy as np, pystac, odc.stac, rasterio
from affine import Affine
from odc.geo.geobox import GeoBox
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fc = json.load(open(os.path.join(ROOT, "fixtures", "items.json")))
C = [pystac.Item.from_dict(f) for f in fc["features"] if f["id"] == "C"]
print("odc.stac", odc.stac.__version__, "GDAL", rasterio.__gdal_version__)
for x0 in (178.5, 178.2):
    g = GeoBox((5000, 5000), Affine(0.0003, 0, x0, 0, -0.0003, -15.0), "epsg:4326")
    for chunks in (None, {"x": 2048, "y": 2048}, {"x": 5000, "y": 5000}):
        kw = {} if chunks is None else {"chunks": chunks}
        a = odc.stac.load(C, bands=["red"], geobox=g, **kw).red.values[0]
        v = a != 0; c = np.where(v.any(0))[0]
        print(f"lon {x0}..{x0 + 1.5:.1f} chunks={chunks}: valid {v.sum()} of 25000000, "
              f"east edge {x0 + (c.max() + 1) * 0.0003:.4f}")
