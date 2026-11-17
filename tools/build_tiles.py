"""Cut NASA/USGS high-resolution landing-site imagery into globe tiles for news.bsvkey.com.

Tiles use Cesium's GeographicTilingScheme: level z is 2^(z+1) x 2^z tiles of 256 px over the whole
body, x from 180 W eastward, y from the north pole down. Output (one folder per site):
    <out>/<body>/<site>/<z>/<x>/<y>.webp      greyscale + alpha (transparent where there is no data)
    <out>/index.json                          sites, bounds, levels, tile lists, credits

Each site: "outer" levels cover the whole source area; "detail" levels only where it matters
(near Perseverance's route at Jezero, around the landing site on the Moon), to stay within
GitHub Pages' ~1 GB site limit. Usage (from the Enceladus pixi env):
    python -I build_tiles.py <work_dir> [site ...]
"""
import json, math, os, subprocess, sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from PIL import Image

WORK = Path(sys.argv[1])
SRC, TMP, OUT = WORK / "sources", WORK / "tmp", WORK / "out"
ENV = Path(sys.executable).parent / "Library"
os.environ.setdefault("GDAL_DATA", str(ENV / "share" / "gdal"))
os.environ.setdefault("PROJ_DATA", str(ENV / "share" / "proj"))
GDALWARP = str(ENV / "bin" / "gdalwarp.exe")

MARS_R, MOON_R = 3396190.0, 1737400.0
SITES = {
    "jezero": dict(body="mars", name="Jezero crater (Perseverance)", radius=MARS_R,
        src=str(SRC / "JEZ_hirise_soc_006_orthoMosaic_25cm.tif"),
        outer=(77.2227, 18.3068, 77.5838, 18.6693), outer_levels=(9, 14), detail_levels=(15, 17),
        detail_route=str(WORK / "m20_traverse.json"), route_buffer_m=350, res="25 cm",
        credit="HiRISE 25 cm orthomosaic for Mars 2020 Terrain Relative Navigation (NASA/JPL/University of Arizona/USGS)",
        date="2020"),
    "apollo11": dict(body="moon", name="Apollo 11, Tranquility Base", radius=MOON_R,
        src="/vsicurl/https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/BROWSE/NAC_ROI/APOLLO11HIB/NAC_ROI_APOLLO11HIB_E006N0235.TIF",
        outer=(23.452, 0.40, 23.538, 0.95), outer_levels=(9, 14), detail_levels=(15, 16),
        detail_box=(23.452, 0.62, 23.538, 0.73), res="40 cm",
        credit="LROC NAC controlled mosaic NAC_ROI_APOLLO11HIB, 0.40 m/px (NASA/GSFC/Arizona State University)", date=""),
    "apollo15": dict(body="moon", name="Apollo 15, Hadley-Apennine", radius=MOON_R,
        src="/vsicurl/https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/BROWSE/NAC_ROI/APOLLO15HIA/NAC_ROI_APOLLO15HIA_E259N0036.TIF",
        outer=(3.40, 25.95, 3.90, 26.32), outer_levels=(9, 13), detail_levels=(14, 15),
        detail_box=(3.52, 26.06, 3.74, 26.20), res="55 cm",
        credit="LROC NAC mosaic NAC_ROI_APOLLO15HIA, 0.55 m/px (NASA/GSFC/Arizona State University)", date=""),
    "apollo17": dict(body="moon", name="Apollo 17, Taurus-Littrow", radius=MOON_R,
        src="/vsicurl/https://asc-astropedia.s3.us-west-2.amazonaws.com/Moon/Apollo/Traverse/Apollo17/ancillary/APOLLO17_ORTHOMOSAIC_50CM.TIFF",
        outer=None, outer_levels=(9, 14), detail_levels=(15, 16),
        detail_box=(30.69, 20.15, 30.85, 20.23), res="50 cm",
        credit="LROC NAC orthomosaic, 50 cm, Apollo 17 traverse map (USGS Astrogeology; NASA/GSFC/ASU)", date=""),
}

def deg_per_px(z): return 360.0 / (2 ** (z + 1) * 256)
def tile_range(lon0, lat0, lon1, lat1, z):
    s = 360.0 / 2 ** (z + 1)
    return (int(math.floor((lon0 + 180) / s)), int(math.floor((90 - lat1) / s)),
            int(math.floor((lon1 + 180) / s - 1e-9)), int(math.floor((90 - lat0) / s - 1e-9)))
def tile_bounds(x, y, z):
    s = 360.0 / 2 ** (z + 1)
    return (-180 + x * s, 90 - (y + 1) * s, -180 + (x + 1) * s, 90 - y * s)

def warp(site, key, bounds, z, out):
    """Reproject the needed part of the source into lat/lon at level z's resolution."""
    if out.exists(): return out
    srs = f"+proj=longlat +R={site['radius']} +no_defs"
    d = deg_per_px(z)
    cmd = [GDALWARP, "-q", "-overwrite", "-t_srs", srs, "-te", *map(str, bounds), "-tr", str(d), str(d),
           "-r", "bilinear", "-dstalpha", "-multi", "-wo", "NUM_THREADS=ALL_CPUS",
           "--config", "GDAL_HTTP_USERAGENT", "Mozilla/5.0 (BSVKey planet imagery)", "--config", "GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR",
           "--config", "GDAL_HTTP_MAX_RETRY", "5", "--config", "GDAL_HTTP_RETRY_DELAY", "3", "--config", "CPL_VSIL_CURL_CHUNK_SIZE", "4194304",
           "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "BIGTIFF=YES", site["src"], str(out)]
    print("  warp", key, bounds, "level", z, flush=True)
    subprocess.run(cmd, check=True)
    return out

def stretch_of(path):
    with rasterio.open(path) as d:
        f = max(1, max(d.width, d.height) // 2000)
        a = d.read(1, out_shape=(d.height // f, d.width // f), resampling=Resampling.average)
        m = d.read(2, out_shape=(d.height // f, d.width // f), resampling=Resampling.nearest) > 0
    v = a[m]
    return (float(np.percentile(v, 0.5)), float(np.percentile(v, 99.5))) if v.size else (0.0, 255.0)

def route_mask_fn(site):
    """True for tiles within route_buffer_m of the rover's route."""
    if "detail_route" not in site: return None
    j = json.load(open(site["detail_route"], encoding="utf-8"))
    pts = []
    for f in j["features"]:
        g = f["geometry"]; lines = [g["coordinates"]] if g["type"] == "LineString" else g["coordinates"]
        for l in lines: pts.extend((p[0], p[1]) for p in l)
    pts = np.array(pts[::3])
    m_per_deg = site["radius"] * math.pi / 180
    buf = site["route_buffer_m"] / m_per_deg
    def near(b):
        lon0, lat0, lon1, lat1 = b
        cx = np.clip(pts[:, 0], lon0, lon1); cy = np.clip(pts[:, 1], lat0, lat1)
        d = np.hypot((pts[:, 0] - cx) * math.cos(math.radians((lat0 + lat1) / 2)), pts[:, 1] - cy)
        return bool((d <= buf).any())
    return near

def cut(site_key, site, src_path, levels, region, keep, stretch, tiles):
    lo, hi = stretch
    with rasterio.open(src_path) as d:
        for z in range(levels[0], levels[1] + 1):
            x0, y0, x1, y1 = tile_range(*region, z); n = 0
            for x in range(x0, x1 + 1):
                for y in range(y0, y1 + 1):
                    b = tile_bounds(x, y, z)
                    if keep and not keep(b): continue
                    win = from_bounds(*b, transform=d.transform)
                    a = d.read([1, 2], window=win, out_shape=(2, 256, 256), boundless=True, fill_value=0,
                               resampling=Resampling.average)
                    if a[1].max() == 0: continue
                    g = np.clip((a[0].astype(np.float32) - lo) * (255.0 / max(hi - lo, 1)), 0, 255).astype(np.uint8)
                    p = OUT / site["body"] / site_key / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                    Image.fromarray(np.dstack([g, g, g, a[1]]), "RGBA").save(p / f"{y}.webp", "WEBP", quality=72, method=5)
                    tiles.append(f"{z}/{x}/{y}"); n += 1
            print(f"  level {z}: {n} tiles", flush=True)

def build(key):
    site = SITES[key]; TMP.mkdir(parents=True, exist_ok=True)
    print(key, flush=True)
    if site["outer"] is None:   # whole source area, in lon/lat
        with rasterio.open(site["src"]) as d:
            from rasterio.warp import transform_bounds
            site["outer"] = tuple(round(v, 4) for v in transform_bounds(d.crs, f"+proj=longlat +R={site['radius']} +no_defs", *d.bounds))
    tiles = []
    outer = warp(site, key, site["outer"], site["outer_levels"][1], TMP / f"{key}_outer.tif")
    stretch = stretch_of(outer)
    cut(key, site, outer, site["outer_levels"], site["outer"], None, stretch, tiles)
    if "detail_route" in site:
        route = route_mask_fn(site)
        o = site["outer"]; rb = json.load(open(site["detail_route"], encoding="utf-8"))
        xs = [p[0] for f in rb["features"] for l in ([f["geometry"]["coordinates"]] if f["geometry"]["type"] == "LineString" else f["geometry"]["coordinates"]) for p in l]
        ys = [p[1] for f in rb["features"] for l in ([f["geometry"]["coordinates"]] if f["geometry"]["type"] == "LineString" else f["geometry"]["coordinates"]) for p in l]
        pad = 0.01
        box = (max(o[0], min(xs) - pad), max(o[1], min(ys) - pad), min(o[2], max(xs) + pad), min(o[3], max(ys) + pad))
        keep = route
    else:
        box, keep = site["detail_box"], None
    detail = warp(site, key, box, site["detail_levels"][1], TMP / f"{key}_detail.tif")
    cut(key, site, detail, site["detail_levels"], box, keep, stretch, tiles)
    size = sum(f.stat().st_size for f in (OUT / site["body"] / key).rglob("*.webp"))
    print(f"  {key}: {len(tiles)} tiles, {size/1e6:.1f} MB", flush=True)
    return dict(key=key, body=site["body"], name=site["name"], res=site.get("res", ""), outer=site["outer"], detail=list(box),
                levels=[site["outer_levels"][0], site["detail_levels"][1]], credit=site["credit"], tiles=tiles, mb=round(size / 1e6, 1))

if __name__ == "__main__":
    keys = sys.argv[2:] or list(SITES)
    idx_path = OUT / "index.json"
    idx = json.load(open(idx_path)) if idx_path.exists() else {"sites": {}}
    for k in keys:
        idx["sites"][k] = build(k)
        OUT.mkdir(parents=True, exist_ok=True)
        json.dump(idx, open(idx_path, "w"), separators=(",", ":"))
    print("done", {k: v["mb"] for k, v in idx["sites"].items()})
