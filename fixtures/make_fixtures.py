"""Synthetic Landsat-like scenes around WGS84GRID30[49, 239].

Three 30 m UTM scenes (uint16, nodata 0, DEFLATE, 512 tiles, COG layout):
  A  2024-12-28  EPSG:32760  lon 178.3 .. 179.6
  B  2024-12-28  EPSG:32701  lon 179.3 .. 180.2 (crosses the antimeridian)
  C  2024-12-29  EPSG:32760  lon 178.4 .. 180.1
Values are a smooth function of lon/lat, so any engine that samples the right
place gets the right number. Scene B is offset by +500 so the paint order of
the A/B overlap is visible in the result.

Writes fixtures/data/*.tif and fixtures/items.json (STAC items with proj:*).
"""
import json, os, sys
import numpy as np
import rasterio
from rasterio.transform import from_origin
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

LAT = (-16.7, -14.8)
RES = 30.0

def field(lon, lat):
    lon = np.where(lon < 0, lon + 360, lon)
    return (2000 + 400 * np.sin(np.radians(lon - 178) * 40)
            + 300 * np.cos(np.radians(lat + 15) * 55))

def scene(name, epsg, lon0, lon1, date, offset):
    fwd = Transformer.from_crs(4326, epsg, always_xy=True)
    inv = Transformer.from_crs(epsg, 4326, always_xy=True)
    lo, la = np.meshgrid(np.linspace(lon0, lon1, 200), np.linspace(*LAT, 200))
    x, y = fwd.transform(lo.ravel(), la.ravel())
    xmin = np.floor(x.min() / RES) * RES
    xmax = np.ceil(x.max() / RES) * RES
    ymin = np.floor(y.min() / RES) * RES
    ymax = np.ceil(y.max() / RES) * RES
    ncol, nrow = int((xmax - xmin) / RES), int((ymax - ymin) / RES)
    tr = from_origin(xmin, ymax, RES, RES)
    v = np.zeros((nrow, ncol), dtype="uint16")
    step = 512
    for r0 in range(0, nrow, step):
        rr = np.arange(r0, min(r0 + step, nrow))
        cx = xmin + (np.arange(ncol) + 0.5) * RES
        cy = ymax - (rr + 0.5) * RES
        X, Y = np.meshgrid(cx, cy)
        lon, lat = inv.transform(X, Y)
        f = field(lon, lat) + offset
        # scene footprint is the lon/lat box, so the corners are nodata
        lonw = np.where(lon < 0, lon + 360, lon)
        inside = (lonw >= lon0) & (lonw <= lon1) & (lat >= LAT[0]) & (lat <= LAT[1])
        v[rr, :] = np.where(inside, np.round(f), 0).astype("uint16")
    path = os.path.join(DATA, f"{name}.tif")
    with rasterio.open(path, "w", driver="COG", width=ncol, height=nrow, count=1,
                       dtype="uint16", crs=f"EPSG:{epsg}", transform=tr, nodata=0,
                       compress="DEFLATE", blocksize=512, overview_resampling="average") as d:
        d.write(v, 1)
    # footprint in lon/lat for the item geometry
    b = [lon0, LAT[0], lon1, LAT[1]]
    if lon1 > 180:
        geom = {"type": "MultiPolygon", "coordinates": [
            [[[lon0, LAT[0]], [180, LAT[0]], [180, LAT[1]], [lon0, LAT[1]], [lon0, LAT[0]]]],
            [[[-180, LAT[0]], [lon1 - 360, LAT[0]], [lon1 - 360, LAT[1]], [-180, LAT[1]], [-180, LAT[0]]]]]}
        bbox = [lon0, LAT[0], lon1 - 360, LAT[1]]
    else:
        geom = {"type": "Polygon", "coordinates": [[[lon0, LAT[0]], [lon1, LAT[0]],
                [lon1, LAT[1]], [lon0, LAT[1]], [lon0, LAT[0]]]]}
        bbox = b
    return {
        "type": "Feature", "stac_version": "1.0.0",
        "stac_extensions": ["https://stac-extensions.github.io/projection/v1.1.0/schema.json"],
        "id": name, "bbox": bbox, "geometry": geom,
        "properties": {"datetime": f"{date}T22:10:00Z", "proj:epsg": epsg,
                       "proj:shape": [nrow, ncol], "proj:transform": list(tr)[:6]},
        "assets": {"red": {"href": os.path.relpath(path, os.path.dirname(HERE)), "type": "image/tiff; application=geotiff; profile=cloud-optimized",
                           "roles": ["data"], "raster:bands": [{"nodata": 0, "data_type": "uint16"}]}},
        "links": [],
    }

items = [scene("A", 32760, 178.3, 179.6, "2024-12-28", 0),
         scene("B", 32701, 179.3, 180.2, "2024-12-28", 500),
         scene("C", 32760, 178.4, 180.1, "2024-12-29", 0)]
with open(os.path.join(HERE, "items.json"), "w") as f:
    json.dump({"type": "FeatureCollection", "features": items}, f, indent=1)
for it in items:
    print(it["id"], it["properties"]["proj:epsg"], it["properties"]["proj:shape"])
