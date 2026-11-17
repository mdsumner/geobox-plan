"""HLS L30 granules over WGS84GRID30[49, 239], from CMR-STAC.

Writes hls/items.json (a FeatureCollection of the full items). Searching
needs no login; reading the files does (Earthdata, see hls/README.md).
"""
import json, os
import pystac_client

HERE = os.path.dirname(os.path.abspath(__file__))
BBOX = [178.5, -16.5, 180.0, -15.0]
DATES = os.environ.get("HLS_DATES", "2024-12-27/2024-12-30")

cat = pystac_client.Client.open("https://cmr.earthdata.nasa.gov/stac/LPCLOUD")
items = cat.search(collections=["HLSL30.v2.0"], bbox=BBOX, datetime=DATES).item_collection()
fc = {"type": "FeatureCollection", "features": [i.to_dict() for i in items]}
with open(os.path.join(HERE, "items.json"), "w") as f:
    json.dump(fc, f, indent=1)
for i in items:
    p = i.properties
    print(i.id, p.get("datetime") or p.get("start_datetime"), i.bbox,
          "proj:" + ",".join(k for k in p if k.startswith("proj:")) or "none")
print(len(items), "items ->", os.path.join(HERE, "items.json"))
