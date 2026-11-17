# planet-imagery

High-resolution imagery of Mars and Moon landing sites, cut into map tiles for the globes on
[news.bsvkey.com](https://news.bsvkey.com/) (Moon and Mars tabs), in the same way the Enceladus
globe shows its Cassini close-ups. Served with GitHub Pages at `https://bsvkey.github.io/planet-imagery/`.

| Site | Body | Source | Detail |
|---|---|---|---|
| Jezero crater (Perseverance) | Mars | HiRISE orthomosaic for Mars 2020 Terrain Relative Navigation, USGS Astrogeology | 25 cm/px along Perseverance's route, about 2.5 m elsewhere |
| Apollo 11, Tranquility Base | Moon | LROC NAC controlled mosaic NAC_ROI_APOLLO11HIB | 40 cm/px |
| Apollo 15, Hadley-Apennine | Moon | LROC NAC mosaic NAC_ROI_APOLLO15HIA | 55 cm/px |
| Apollo 17, Taurus-Littrow | Moon | LROC NAC orthomosaic, Apollo 17 traverse map, USGS Astrogeology | 50 cm/px |

## Layout

- `<body>/<site>/<z>/<x>/<y>.webp`: 256 px tiles in Cesium's geographic tiling scheme (level `z` has
  2^(z+1) × 2^z tiles over the whole body; x from 180° W eastward, y from the north pole down).
  Greyscale with alpha: transparent where the source has no data.
- `index.json`: each site's bounds, zoom levels, tile list and credit.

## Credits and licence

Imagery: NASA/JPL/University of Arizona (HiRISE); NASA/GSFC/Arizona State University (LROC);
processed products from the U.S. Geological Survey Astrogeology Science Center. NASA and USGS
imagery is generally not subject to copyright in the United States; credit the sources above when
you reuse it. The tiles here are reprojected, contrast-stretched and compressed for web display and
are not science data; use the original products for measurement.

## Terrain heights and Moon topography

- `terrain/<body>/<z>/<x>/<y>.bin`: 65 x 65 little-endian Int16 heights in metres, levels 0 to 5, using the same geographic tiling as the imagery. Moon: LRO LOLA global DEM 118 m. Mars: MGS MOLA global DEM 463 m (both NASA/GSFC, via USGS Astrogeology).
- `relief/moon/<z>/<x>/<y>.jpg`: colour shaded relief made from LOLA, levels 0 to 4.
- Built with `tools/build_terrain.py`.

## Feature names

- `names/<body>.json`: IAU-approved feature names from the USGS Gazetteer of Planetary Nomenclature (name, type code, lat, lon, diameter km, origin, gazetteer id), largest first. Built with `tools/build_names.py`.
