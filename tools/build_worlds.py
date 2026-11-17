"""Globe tiles for more bodies (planets, moons) for news.bsvkey.com, from USGS Astrogeology global mosaics.

For each body: reads the global mosaic through HTTP range requests (every Nth row only), turns it into a
lon/lat grid (-180..180, 90..-90), and writes under <out>/worlds/<key>/:
    <z>/<x>/<y>.jpg   256 px tiles, Cesium geographic tiling (2^(z+1) x 2^z tiles), levels 0..maxLevel
    full.jpg          2048 x 1024 whole-body image (for the Orbit page's spheres)
and records it in <out>/index.json under "worlds".
Usage (Enceladus pixi env on PATH): python -I build_worlds.py <work_dir> [key ...]
"""
import io, json, math, sys, zipfile, urllib.request
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from PIL import Image

WORK = Path(sys.argv[1]); OUT = WORK / "out"
S3 = "https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/"
SSS = "https://www.solarsystemscope.com/textures/download/"
W = {   # key: name, mean radius km, source, credit
    "mercury": ("Mercury", 2439.7, S3 + "Mercury_MESSENGER_MDIS_Basemap_EnhancedColor_Mosaic_Global_665m.tif", "MESSENGER MDIS enhanced colour mosaic, 665 m (NASA/JHUAPL/Carnegie), via USGS"),
    "venus": ("Venus", 6051.8, S3 + "Venus_Magellan_C3-MDIR_Colorized_Global_Mosaic_4641m.tif", "Magellan radar mosaic, colourised, 4.6 km (NASA/JPL), via USGS"),
    "phobos": ("Phobos", 11.08, S3 + "Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif", "Viking Orbiter mosaic, DLR control (NASA/JPL/DLR), via USGS"),
    "deimos": ("Deimos", 6.2, None, "No global map of Deimos exists; shown plain"),
    "io": ("Io", 1821.6, S3 + "Io_GalileoSSI-Voyager_Global_Mosaic_ClrMerge_1km.tif", "Galileo SSI and Voyager colour mosaic, 1 km (NASA/JPL), via USGS"),
    "europa": ("Europa", 1560.8, S3 + "Europa_Voyager_GalileoSSI_global_mosaic_500m.tif", "Voyager and Galileo SSI mosaic, 500 m (NASA/JPL), via USGS"),
    "ganymede": ("Ganymede", 2634.1, S3 + "Ganymede_Voyager_GalileoSSI_Global_ClrMosaic_1435m.tif", "Voyager and Galileo SSI colour mosaic, 1.4 km (NASA/JPL), via USGS"),
    "callisto": ("Callisto", 2410.3, S3 + "Callisto_Voyager_GalileoSSI_global_mosaic_1km.tif", "Voyager and Galileo SSI mosaic, 1 km (NASA/JPL), via USGS"),
    "titan": ("Titan", 2574.7, S3 + "Titan_ISS_P19658_Mosaic_Global_4km.tif", "Cassini ISS near-infrared mosaic, 4 km (NASA/JPL/SSI), via USGS"),
    "enceladus": ("Enceladus", 252.1, S3 + "Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif", "Cassini ISS mosaic, 100 m (NASA/JPL/SSI), via USGS"),
    "mimas": ("Mimas", 198.2, "zip:" + S3 + "Mimas/Cassini_DLR_Mimas.zip", "Cassini ISS mosaic (NASA/JPL/SSI/DLR), via USGS"),
    "tethys": ("Tethys", 531.1, S3 + "Tethys_Cassini_mosaic_global_293m.tif", "Cassini ISS mosaic, 293 m (NASA/JPL/SSI), via USGS"),
    "dione": ("Dione", 561.4, S3 + "Dione_Cassini_Voyager_mosaic_global_154m.tif", "Cassini and Voyager mosaic, 154 m (NASA/JPL/SSI), via USGS"),
    "rhea": ("Rhea", 763.8, S3 + "Rhea_Cassini_Voyager_mosaic_global_417m.tif", "Cassini and Voyager mosaic, 417 m (NASA/JPL/SSI), via USGS"),
    "iapetus": ("Iapetus", 734.5, S3 + "Iapetus_Cassini_Voyager_mosaic_global_783m.tif", "Cassini and Voyager mosaic, 783 m (NASA/JPL/SSI), via USGS"),
    "jupiter": ("Jupiter", 69911, "img:" + SSS + "2k_jupiter.jpg", "Solar System Scope texture (CC BY 4.0), from NASA Juno/Cassini imagery"),
    "saturn": ("Saturn", 58232, "img:" + SSS + "2k_saturn.jpg", "Solar System Scope texture (CC BY 4.0), from NASA Cassini imagery"),
}
ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_USERAGENT="Mozilla/5.0 (BSVKey planet imagery)",
           GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="3", CPL_VSIL_CURL_CHUNK_SIZE="262144")

def to8(a):
    """Any numeric band(s) -> uint8 with a 0.5-99.5 percentile stretch (8-bit data is kept as is)."""
    if a.dtype == np.uint8: return a
    a = a.astype(np.float32); v = a[np.isfinite(a)]
    lo, hi = np.percentile(v, 0.5), np.percentile(v, 99.5)
    return np.clip((a - lo) * 255 / max(hi - lo, 1e-6), 0, 255).astype(np.uint8)

def read_global(src, width):
    """Read a global equirectangular raster decimated to (width x width/2), columns from 180 W, rows from 90 N."""
    with rasterio.Env(**ENV), rasterio.open(src) as d:
        bands = [1, 2, 3] if d.count >= 3 else [1]
        H = width // 2
        b = d.bounds; crs = d.crs.to_dict() if d.crs else {}
        print(f"  {d.width} x {d.height} x{d.count} {d.dtypes[0]} bounds {tuple(round(x) for x in b)} crs {crs}", flush=True)
        # longitude of the left edge and the span covered, from the projection
        if d.crs and d.crs.is_geographic: lon0, span_x, lat_top, span_y = b.left, b.right - b.left, b.top, b.top - b.bottom
        else:
            R = crs.get("R") or crs.get("a") or (b.right - b.left) / (2 * math.pi)
            k = R * math.pi / 180; c = crs.get("lon_0", 0)
            lon0, span_x, lat_top, span_y = c + b.left / k, (b.right - b.left) / k, b.top / k, (b.top - b.bottom) / k
        # read at a resolution matching the target grid
        rw = max(16, round(width * span_x / 360))
        rh = max(8, round(H * span_y / 180))
        a = d.read(bands, out_shape=(len(bands), rh, rw), resampling=Resampling.nearest, masked=True)
        nod = a.mask if np.ma.isMaskedArray(a) else None
        a = to8(np.ma.filled(a, 0) if nod is not None else a)
    if a.shape[0] == 1: a = np.repeat(a, 3, 0)
    img = np.zeros((3, H, width), np.uint8)
    # place into the global grid (rows from lat_top, columns from lon0, wrapping in longitude)
    r0 = round((90 - lat_top) / 180 * H); c0 = round((lon0 + 180) / 360 * width)
    rows = slice(max(r0, 0), min(r0 + a.shape[1], H))
    for i in range(a.shape[2]):
        img[:, rows, (c0 + i) % width] = a[:, rows.start - r0: rows.stop - r0, i]
    return np.moveaxis(img, 0, -1)

def get_source(key, src):
    if src is None:   # plain grey with a little texture
        rng = np.random.default_rng(7); g = (120 + rng.normal(0, 6, (1024, 2048))).clip(0, 255).astype(np.uint8)
        return np.stack([g, g, g], -1), 2
    if src.startswith("img:"):
        im = Image.open(io.BytesIO(urllib.request.urlopen(urllib.request.Request(src[4:], headers={"User-Agent": "Mozilla/5.0"})).read())).convert("RGB")
        return np.asarray(im.resize((2048, 1024), Image.LANCZOS)), 2
    if src.startswith("zip:"):
        z = WORK / "tmp" / (key + ".zip"); z.parent.mkdir(parents=True, exist_ok=True)
        if not z.exists(): urllib.request.urlretrieve(src[4:], z)
        with zipfile.ZipFile(z) as f:
            name = max((n for n in f.namelist() if n.lower().endswith((".tif", ".tiff"))), key=lambda n: f.getinfo(n).file_size)
            print("  zip member", name, flush=True)
            dst = WORK / "tmp" / Path(name).name
            if not dst.exists(): dst.write_bytes(f.read(name))
        src = str(dst)
    else:
        local = WORK / "tmp" / Path(src).name   # downloaded in full beforehand (the server drops long range reads)
        src = str(local) if local.exists() else "/vsicurl/" + src
    with rasterio.Env(**ENV), rasterio.open(src) as d: full_w = d.width
    zmax = max(2, min(4, int(math.floor(math.log2(full_w / 512)))))
    return read_global(src, 256 * 2 ** (zmax + 1)), zmax

def build(key):
    name, r_km, src, credit = W[key]; print(key, flush=True)
    a, zmax = get_source(key, src)
    img = Image.fromarray(a); base = OUT / "worlds" / key; n = 0
    for z in range(zmax + 1):
        nx, ny = 2 ** (z + 1), 2 ** z
        lvl = img if img.size == (nx * 256, ny * 256) else img.resize((nx * 256, ny * 256), Image.LANCZOS)
        for x in range(nx):
            for y in range(ny):
                p = base / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                lvl.crop((x * 256, y * 256, x * 256 + 256, y * 256 + 256)).save(p / f"{y}.jpg", quality=84); n += 1
    img.resize((2048, 1024), Image.LANCZOS).save(base / "full.jpg", quality=86)
    if key == "saturn":
        ring = urllib.request.urlopen(urllib.request.Request(SSS + "2k_saturn_ring_alpha.png", headers={"User-Agent": "Mozilla/5.0"})).read()
        (base / "ring.png").write_bytes(ring)
    size = sum(f.stat().st_size for f in base.rglob("*.*")) / 1e6
    print(f"  {key}: levels 0..{zmax}, {n} tiles, {size:.1f} MB", flush=True)
    return dict(name=name, radius=r_km, maxLevel=zmax, credit=credit, plain=src is None)

if __name__ == "__main__":
    idx_path = OUT / "index.json"; idx = json.load(open(idx_path))
    idx.setdefault("worlds", {})
    for k in (sys.argv[2:] or list(W)):
        try: idx["worlds"][k] = build(k)
        except Exception as e: print("  FAILED", k, e, flush=True); continue
        json.dump(idx, open(idx_path, "w"), separators=(",", ":"))
    print("done", list(idx["worlds"]))
