"""Earth feature names for the globes, from Natural Earth (public domain), in the same compact format as the IAU
names: <out>/names/earth.json {"credit", "types": {code: name}, "f": [[name, code, lat, lon, size_km, note, wikidata], ...]}
largest first. Oceans, seas, bays, mountain ranges, deserts, plains, islands, peaks and major cities.
Usage: python -I -X utf8 build_earth_names.py <work_dir>"""
import json, math, sys
from pathlib import Path
WORK = Path(sys.argv[1]); SRC = WORK / "names" / "earth"; OUT = WORK / "out" / "names"
TYPES = {}
def code(cla):
    c = "".join(w[0] for w in cla.split()).upper()[:2] + str(len(TYPES))
    for k, v in TYPES.items():
        if v == cla: return k
    TYPES[c] = cla; return c
def ring_centroid(r):
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(r, r[1:] + r[:1]):
        f = x0 * y1 - x1 * y0; a += f; cx += (x0 + x1) * f; cy += (y0 + y1) * f
    if abs(a) < 1e-12: xs, ys = zip(*r); return sum(xs) / len(xs), sum(ys) / len(ys), 0
    return cx / (3 * a), cy / (3 * a), abs(a) / 2
def label_point(g):
    polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
    cs = [ring_centroid(p[0]) for p in polys]; best = max(cs, key=lambda c: c[2])
    area = sum(c[2] for c in cs) * 111.32 ** 2 * math.cos(math.radians(best[1]))   # km^2, roughly
    return best[0], best[1], 2 * math.sqrt(area / math.pi)
def size_km(bbox):
    w, s, e, n = bbox; lat = (s + n) / 2
    return math.hypot((e - w) * 111.32 * math.cos(math.radians(lat)), (n - s) * 111.32)
rows = []
def add(name, cla, lat, lon, d, note, wiki): rows.append([name, code(cla), round(lat, 3), round(lon, 3), round(d, 1), note, wiki or ""])
for f in ("ne_10m_geography_marine_polys", "ne_10m_geography_regions_polys"):
    for ft in json.load(open(SRC / f"{f}.geojson", encoding="utf-8"))["features"]:
        p = ft["properties"]; name = p.get("name") or p.get("NAME") or p.get("LABEL")
        if not name or not ft.get("geometry"): continue
        cla = (p.get("featurecla") or p.get("FEATURECLA") or "Feature").replace("_", " ").title()
        lon, lat, d = label_point(ft["geometry"])   # size from area: bounding boxes break across the 180 degree line
        d *= 1.6 if cla in ("Ocean", "Sea") else 1
        region = p.get("REGION") or p.get("region") or ""
        add(name.title() if name.isupper() else name, cla, lat, lon, d, f"{cla} · {region}".strip(" ·"), p.get("WIKIDATAID") or p.get("wikidataid"))
for ft in json.load(open(SRC / "ne_10m_geography_regions_elevation_points.geojson", encoding="utf-8"))["features"]:
    p = ft["properties"]; name = p.get("NAME") or p.get("name")
    if not name: continue
    lon, lat = ft["geometry"]["coordinates"][:2]; el = p.get("ELEVATION") or p.get("elevation")
    add(name, "Peak", lat, lon, max(5, 60 / (1 + (p.get("SCALERANK") or p.get("scalerank") or 5))), f"Peak · {int(el):,} m" if el else "Peak", p.get("WIKIDATAID") or p.get("wikidataid"))
for ft in json.load(open(SRC / "ne_10m_populated_places_simple.geojson", encoding="utf-8"))["features"]:
    p = ft["properties"]
    if (p.get("scalerank") or 10) > 3: continue
    lon, lat = ft["geometry"]["coordinates"][:2]; pop = p.get("pop_max") or 0
    add(p["name"], "City", lat, lon, max(20, math.sqrt(pop) / 12), f"City · {p.get('adm0name', '')}" + (f" · {pop:,} people" if pop else ""), p.get("wikidataid"))
rows.sort(key=lambda r: -r[4])
json.dump({"credit": "Natural Earth (public domain)", "types": TYPES, "f": rows}, open(OUT / "earth.json", "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
print("earth", len(rows), round((OUT / "earth.json").stat().st_size / 1e6, 2), "MB", list(TYPES.values())[:30])
