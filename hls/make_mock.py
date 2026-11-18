"""HLS-shaped mock items from the fixtures, to dry-run the HLS scripts offline.

Scenes C and B become int16 / nodata -9999 COGs (DEFLATE, 256 tiles) in UTM
*north* zones with negative northings (as HLS writes the southern
hemisphere), and the items carry no proj:* or raster:bands, like CMR-STAC
HLS items. Needs fixtures/make_fixtures.py to have run.

  python3 hls/make_mock.py
  HLS_ITEMS=hls/mock/items.json HLS_OUT=ref_hls_mock python3 hls/odc_reference.py
  HLS_ITEMS=hls/mock/items.json HLS_OUT=ref_hls_mock Rscript hls/run_hls.R
"""
import copy, json, os
import numpy as np
import rasterio
from affine import Affine

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "hls", "mock")
os.makedirs(OUT, exist_ok=True)
fx = json.load(open(os.path.join(ROOT, "fixtures", "items.json")))
src = {f["id"]: f for f in fx["features"]}

feats = []
for sid, tile, epsg_n, date in [("C", "T60KYG", 32660, "2024363T220615"),
                                ("B", "T01KAB", 32601, "2024363T220615")]:
    f = src[sid]
    with rasterio.open(os.path.join(ROOT, f["assets"]["red"]["href"])) as d:
        v = d.read(1).astype("int32"); t = d.transform
    v = np.where(v == 0, -9999, v).astype("int16")
    tn = Affine(t.a, t.b, t.c, t.d, t.e, t.f - 10_000_000)   # UTM S -> UTM N spelling
    iid = f"HLS.L30.{tile}.{date}.v2.0"
    path = os.path.join(OUT, f"{iid}.B04.tif")
    with rasterio.open(path, "w", driver="COG", width=v.shape[1], height=v.shape[0], count=1,
                       dtype="int16", crs=f"EPSG:{epsg_n}", transform=tn, nodata=-9999,
                       compress="DEFLATE", predictor=2, blocksize=256) as d:
        d.write(v, 1)
    g = copy.deepcopy(f)
    g["id"] = iid
    g["collection"] = "HLSL30.v2.0"
    g["properties"] = {"datetime": f["properties"]["datetime"].replace("Z", ".123Z"),
                       "eo:cloud_cover": 88}
    g["stac_extensions"] = []
    g["assets"] = {"B04": {"href": os.path.relpath(path, ROOT), "type": "image/tiff; application=geotiff",
                           "roles": ["data"]}}
    feats.append(g)
    print(iid, epsg_n, v.shape)
json.dump({"type": "FeatureCollection", "features": feats}, open(os.path.join(OUT, "items.json"), "w"), indent=1)
