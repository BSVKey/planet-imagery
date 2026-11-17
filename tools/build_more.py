"""Sharper imagery, extra map styles and 3D terrain for the Worlds globes (after build_worlds.py).

Sources are downloaded in full first (the USGS server drops long range reads) into <work>/tmp, then:
  styles   <out>/worlds/<body>/<style>/<z>/<x>/<y>.jpg  (style "" = the main imagery, replacing the old tiles)
  terrain  <out>/terrain/<body>/<z>/<x>/<y>.bin       65 x 65 Int16 heights (m), levels 0..5, like the Moon and Mars
and records them in index.json: worlds[body].styles = [{name, path, maxLevel, credit}], terrain[body] = {...}.
Usage (Enceladus pixi env on PATH): python -I build_more.py <work_dir> [job ...]
"""
import json, math, sys, subprocess, urllib.request
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from scipy.ndimage import map_coordinates
from PIL import Image
sys.path.insert(0, str(Path(__file__).parent))
import build_worlds as bw

WORK = Path(sys.argv[1]); OUT = WORK / "out"; TMP = WORK / "tmp"; TMP.mkdir(exist_ok=True)
S3 = "https://asc-pds-services.s3.us-west-2.amazonaws.com/mosaic/"
STOOKE = "https://sbnarchive.psi.edu/pds3/non_mission/EAR_A_3_RDR_STOOKEMAPS_V1_0/maps/"
Image.MAX_IMAGE_PIXELS = None

# imagery jobs: body, style path ("" = main), shown name, source, max level, credit
IMG = [
    ("deimos", "", "Imagery", "img:" + STOOKE + "m2deimos/newmap.jpg", 3, "Viking Orbiter mosaic by P. Stooke, C. Jongkind and M. Arntz (1999), NASA PDS Small Bodies Node"),
    ("mercury", "mono", "Natural (greyscale)", S3 + "Mercury_MESSENGER_mosaic_global_250m_2013.tif", 5, "MESSENGER MDIS global mosaic, 250 m (NASA/JHUAPL/Carnegie), via USGS"),
    ("venus", "radar", "Radar", S3 + "Venus_Magellan_C3-MDIR_Global_Mosaic_2025m.tif", 5, "Magellan synthetic-aperture radar mosaic, 2 km (NASA/JPL), via USGS"),
    ("europa", "", "Imagery", S3 + "Europa_Voyager_GalileoSSI_global_mosaic_500m.tif", 5, "Voyager and Galileo SSI mosaic, 500 m (NASA/JPL), via USGS"),
    ("ganymede", "mono", "Detail (greyscale)", S3 + "Ganymede_Voyager_GalileoSSI_global_mosaic_1km.tif", 5, "Voyager and Galileo SSI mosaic, 1 km (NASA/JPL), via USGS"),
    ("dione", "", "Imagery", S3 + "Dione_Cassini_Voyager_mosaic_global_154m.tif", 5, "Cassini and Voyager mosaic, 154 m (NASA/JPL/SSI), via USGS"),
]
# terrain jobs: body, source DEM, height scale to metres, credit
DEM = [
    ("mercury", S3 + "Mercury_Messenger_USGS_DEM_Global_665m_v2.tif", 1.0, "MESSENGER global DEM, 665 m (NASA/JHUAPL/Carnegie/USGS)"),
    ("venus", S3 + "Venus_Magellan_Topography_Global_4641m_v02.tif", 1.0, "Magellan global topography, 4.6 km (NASA/JPL/USGS)"),
    ("phobos", S3 + "Phobos_ME_HRSC_DEM_Global_2ppd.tif", 1.0, "Mars Express HRSC global DEM of Phobos (ESA/DLR/FU Berlin), via USGS"),
]

def fetch(url):
    """Whole file into tmp (curl resumes and retries); returns the local path."""
    dst = TMP / url.rsplit("/", 1)[1]
    for _ in range(8):
        r = subprocess.run(["curl", "-s", "-L", "-C", "-", "--retry", "5", "--retry-delay", "3", "-o", str(dst), url])
        if r.returncode == 0: break
    print(f"  have {dst.name} {dst.stat().st_size / 1e6:.0f} MB", flush=True)
    return dst

def write_tiles(a, base, zmax):
    img = Image.fromarray(a); n = 0
    for z in range(zmax + 1):
        nx, ny = 2 ** (z + 1), 2 ** z
        lvl = img if img.size == (nx * 256, ny * 256) else img.resize((nx * 256, ny * 256), Image.LANCZOS)
        for x in range(nx):
            for y in range(ny):
                p = base / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                lvl.crop((x * 256, y * 256, x * 256 + 256, y * 256 + 256)).save(p / f"{y}.jpg", quality=84); n += 1
    return n

def imagery(job):
    body, style, name, src, zmax, credit = job; print(body, style or "main", flush=True)
    if src.startswith("img:"):
        local = TMP / src.rsplit("/", 1)[1]
        if not local.exists(): local.write_bytes(urllib.request.urlopen(urllib.request.Request(src[4:], headers={"User-Agent": "Mozilla/5.0"})).read())
        im = Image.open(local).convert("RGB"); a = np.asarray(im.resize((256 * 2 ** (zmax + 1), 256 * 2 ** zmax), Image.LANCZOS))
    else:
        a = bw.read_global(str(fetch(src)), 256 * 2 ** (zmax + 1))
    base = OUT / "worlds" / body / style if style else OUT / "worlds" / body
    if not style:   # replacing the main imagery: clear old levels first
        for z in range(8):
            d = base / str(z)
            if d.exists():
                for f in d.rglob("*.jpg"): f.unlink()
    n = write_tiles(a, base, zmax)
    if not style: Image.fromarray(a).resize((2048, 1024), Image.LANCZOS).save(OUT / "worlds" / body / "full.jpg", quality=86)
    print(f"  {n} tiles", flush=True)
    return dict(name=name, path=style, maxLevel=zmax, credit=credit)

def terrain(job):
    body, src, scale, credit = job; print(body, "terrain", flush=True)
    GW, GH, TMAX, N = 4096, 2048, 5, 65
    with rasterio.open(str(fetch(src))) as d:
        a = d.read(1, out_shape=(GH, GW), resampling=Resampling.average if d.width <= 2 * GW else Resampling.nearest, masked=True).astype(np.float32)
        a = np.ma.filled(a, np.nan) * scale
        b = d.bounds; crs = d.crs.to_dict() if d.crs else {}
        print(f"  {d.width} x {d.height}, bounds {tuple(round(x) for x in b)}, crs {crs}", flush=True)
    if np.isnan(a).any(): a[np.isnan(a)] = np.nanmean(a)
    print(f"  heights {a.min():.0f} .. {a.max():.0f} m", flush=True)
    def sample(lons, lats):
        col = (lons + 180) / 360 * GW - 0.5; row = (90 - lats) / 180 * GH - 0.5
        return map_coordinates(a, [np.clip(row, 0, GH - 1), np.mod(col, GW)], order=1, mode="wrap")
    n = 0
    for z in range(TMAX + 1):
        nx, ny = 2 ** (z + 1), 2 ** z; s = 180 / ny
        for x in range(nx):
            for y in range(ny):
                lons = -180 + x * s + np.linspace(0, s, N)[None, :].repeat(N, 0); lats = 90 - y * s - np.linspace(0, s, N)[:, None].repeat(N, 1)
                h = np.clip(np.round(sample(lons, lats)), -32767, 32767).astype("<i2")
                p = OUT / "terrain" / body / str(z) / str(x); p.mkdir(parents=True, exist_ok=True)
                (p / f"{y}.bin").write_bytes(h.tobytes()); n += 1
    print(f"  {n} terrain tiles", flush=True)
    return dict(maxLevel=TMAX, samples=N, credit=credit, relief=None)

if __name__ == "__main__":
    want = set(sys.argv[2:])
    idx_path = OUT / "index.json"
    for job in IMG:
        key = job[0] + ("/" + job[1] if job[1] else "")
        if want and key not in want and job[0] not in want: continue
        try: r = imagery(job)
        except Exception as e: print("  FAILED", key, e, flush=True); continue
        idx = json.load(open(idx_path)); w = idx["worlds"][job[0]]
        if job[1]: w["styles"] = [s for s in w.get("styles", []) if s["path"] != job[1]] + [r]
        else: w.update(maxLevel=r["maxLevel"], credit=r["credit"], plain=False)
        json.dump(idx, open(idx_path, "w"), separators=(",", ":"))
    for job in DEM:
        if want and job[0] + "/terrain" not in want and job[0] not in want: continue
        try: r = terrain(job)
        except Exception as e: print("  FAILED", job[0], "terrain", e, flush=True); continue
        idx = json.load(open(idx_path)); idx.setdefault("terrain", {})[job[0]] = r
        json.dump(idx, open(idx_path, "w"), separators=(",", ":"))
    print("done", flush=True)
