"""Terrain height tiles (and, for the Moon, colour topography tiles) for the globes on news.bsvkey.com.

Reads the global NASA/USGS elevation maps through HTTP range requests, every Nth row only (they are
stored one row per block, so this fetches a few hundred MB instead of several GB):
    Moon: LRO LOLA LDEM 118 m (Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif), values x 0.5 m
    Mars: MGS MOLA DEM 463 m (Mars_MGS_MOLA_DEM_mosaic_global_463m.tif), metres
Writes, under <out>/terrain/<body>/:
    <z>/<x>/<y>.bin   65 x 65 little-endian Int16 heights in metres, levels 0..5, Cesium geographic
                      tiling (2^(z+1) x 2^z tiles; samples include the tile edges), like the Enceladus globe
and for the Moon, <out>/relief/moon/<z>/<x>/<y>.jpg colour shaded relief, levels 0..4.
Usage (Enceladus pixi env on PATH): python -I build_terrain.py <work_dir> moon mars
"""
import json, math, sys
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from scipy.ndimage import map_coordinates
from PIL import Image

WORK = Path(sys.argv[1]); OUT = WORK / "out"
BODIES = {
    "moon": dict(url="https://planetarymaps.usgs.gov/mosaic/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif", scale=0.5, radius=1737400,
                 credit="LRO LOLA global DEM 118 m (NASA/GSFC), via USGS Astrogeology", relief=True),
    "mars": dict(url="https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif", scale=1.0, radius=3396190,
                 credit="MGS MOLA global DEM 463 m (NASA/GSFC), via USGS Astrogeology", relief=False),
}
GW, GH = 8192, 4096          # global grid read from the DEM (about 1.3 km/px Moon, 2.6 km/px Mars)
TMAX, N = 5, 65              # terrain levels 0..5, 65 x 65 samples per tile
RMAX = 4                     # relief imagery levels 0..4 (256 px tiles)

def read_grid(b):
    env = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_USERAGENT="Mozilla/5.0 (BSVKey planet imagery)",
               GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="3", CPL_VSIL_CURL_CHUNK_SIZE="262144")
    with rasterio.Env(**env), rasterio.open("/vsicurl/" + b["url"]) as d:
        print("  reading", d.width, "x", d.height, "->", GW, "x", GH, flush=True)
        a = d.read(1, out_shape=(GH, GW), resampling=Resampling.nearest).astype(np.float32)
        a[a == d.nodata] = np.nan
    a *= b["scale"]
    if np.isnan(a).any():   # fill the few holes from neighbours
        m = np.isnan(a); a[m] = np.nanmean(a)
    return a   # row 0 = 90 N, col 0 = 180 W; pixel centres

def sample(grid, lons, lats):
    """Bilinear heights at lon/lat (degrees) from the global grid, wrapping in longitude."""
    col = (lons + 180.0) / 360.0 * GW - 0.5
    row = (90.0 - lats) / 180.0 * GH - 0.5
    return map_coordinates(grid, [np.clip(row, 0, GH - 1), np.mod(col, GW)], order=1, mode="wrap")

def build(key):
    b = BODIES[key]; print(key, flush=True)
    grid = read_grid(b)
    print(f"  heights {np.nanmin(grid):.0f} .. {np.nanmax(grid):.0f} m", flush=True)
    n = 0
    for z in range(TMAX + 1):
        nx, ny = 2 ** (z + 1), 2 ** z; s = 180.0 / ny
        for x in range(nx):
            for y in range(ny):
                lon0, lat1 = -180 + x * s, 90 - y * s
                lons = lon0 + np.linspace(0, s, N)[None, :].repeat(N, 0)
                lats = lat1 - np.linspace(0, s, N)[:, None].repeat(N, 1)
                h = np.clip(np.round(sample(grid, lons, lats)), -32767, 32767).astype("<i2")
                p = OUT / "terrain" / key / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                (p / f"{y}.bin").write_bytes(h.tobytes()); n += 1
    print(f"  terrain: {n} tiles", flush=True)
    if b["relief"]:
        # hypsometric colours (low = blue-violet, high = warm) x hillshade from the north-west
        lo, hi = np.percentile(grid, 1), np.percentile(grid, 99)
        stops = np.array([[0, 46, 52, 110], [0.25, 60, 110, 170], [0.45, 120, 170, 140], [0.65, 215, 200, 120], [0.85, 200, 120, 70], [1, 245, 240, 235]], float)
        mpp = 2 * math.pi * b["radius"] / GW
        gy, gx = np.gradient(grid, mpp)
        az, alt = math.radians(315), math.radians(35)
        shade = np.clip(math.sin(alt) - np.cos(alt) * (np.sin(az) * gx * 6 - np.cos(az) * gy * 6) / np.sqrt(1 + (gx * 6) ** 2 + (gy * 6) ** 2), 0, 1)
        t = np.clip((grid - lo) / (hi - lo), 0, 1)
        rgb = np.stack([np.interp(t, stops[:, 0], stops[:, k]) for k in (1, 2, 3)], -1) * (0.35 + 0.65 * shade[..., None])
        rgb = np.clip(rgb, 0, 255).astype(np.uint8)
        img = Image.fromarray(rgb)
        m = 0
        for z in range(RMAX + 1):
            nx, ny = 2 ** (z + 1), 2 ** z
            lvl = img.resize((nx * 256, ny * 256), Image.LANCZOS) if nx * 256 != GW else img
            for x in range(nx):
                for y in range(ny):
                    p = OUT / "relief" / key / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                    lvl.crop((x * 256, y * 256, x * 256 + 256, y * 256 + 256)).save(p / f"{y}.jpg", quality=82)
                    m += 1
        print(f"  relief: {m} tiles", flush=True)
    return dict(maxLevel=TMAX, samples=N, credit=b["credit"], relief=RMAX if b["relief"] else None)

if __name__ == "__main__":
    idx_path = OUT / "index.json"; idx = json.load(open(idx_path))
    idx.setdefault("terrain", {})
    for k in sys.argv[2:]:
        idx["terrain"][k] = build(k)
        json.dump(idx, open(idx_path, "w"), separators=(",", ":"))
    print("done", list(idx["terrain"]))
