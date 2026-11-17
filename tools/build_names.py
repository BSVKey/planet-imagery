"""IAU feature names (USGS Gazetteer of Planetary Nomenclature) as compact JSON for the globes.
Writes <out>/names/<body>.json: {"credit", "types": {code: name}, "f": [[name, code, lat, lon, diam_km, origin, id], ...]}
sorted largest first. Usage: python -I -X utf8 build_names.py <work_dir> moon mars"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from read_dbf import read_dbf
WORK = Path(sys.argv[1]); OUT = WORK / "out" / "names"; OUT.mkdir(parents=True, exist_ok=True)
for body in sys.argv[2:]:
    rows, types = [], {}
    for r in read_dbf(WORK / "names" / body / f"{body.upper()}_nomenclature_center_pts.dbf"):
        if not r["approval"].startswith("Adopted"): continue
        lon = float(r["center_lon"]); lon = lon - 360 if lon > 180 else lon
        types[r["code"]] = r["type"].split(",")[0]
        origin = r["origin"].strip()
        if len(origin) > 180: origin = origin[:177].rsplit(" ", 1)[0] + "…"
        link = r["link"].rsplit("/", 1)[-1]
        rows.append([r["clean_name"] or r["name"], r["code"], round(float(r["center_lat"]), 3), round(lon, 3),
                     round(float(r["diameter"] or 0), 1), origin, int(link) if link.isdigit() else 0])
    rows.sort(key=lambda x: -x[4])
    json.dump({"credit": "IAU / USGS Gazetteer of Planetary Nomenclature", "types": types, "f": rows},
              open(OUT / f"{body}.json", "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
    print(body, len(rows), round((OUT / f"{body}.json").stat().st_size / 1e6, 2), "MB", sorted(types.items())[:40])
