# Kalimantan Wildfire Monitoring

Satellite-based wildfire monitoring across all five Kalimantan provinces (Indonesian Borneo) using Google Earth Engine. Integrates VIIRS active-fire detections, Landsat 8/9 burn severity mapping, and Sentinel-1 SAR cloud-penetrating change detection for the August 2026 fire season.

---

## Live Apps & Demo

| App Version | Demo |
|---|---|
| **Desktop** ([Open App](https://mwahyur46.users.earthengine.app/view/kalimantan-fire-monitoring-desktop)) | [![Watch on YouTube](https://img.youtube.com/vi/mK6-vJiLtpg/hqdefault.jpg)](https://youtu.be/mK6-vJiLtpg) |
| **Mobile** ([Open App](https://mwahyur46.users.earthengine.app/view/kalimantan-fire-monitoring)) | [![Watch on YouTube](https://img.youtube.com/vi/HjclsARBXyo/hqdefault.jpg)](https://youtube.com/shorts/HjclsARBXyo) |

---

## Scripts

| File | Purpose |
|---|---|
| `kalimantan_wildfire_monitoring_aug2026.js` | Desktop GEE App -- three-panel layout with layer toggles, legends, opacity sliders, pixel inspector, and async statistics panel |
| `kalimantan_wildfire_monitoring_aug2026_mobile.js` | Mobile GEE App -- full-screen map with floating bottom tab bar; identical data pipeline, portrait-optimized UI |

Both scripts share the same data sources, burn index methodology, and configuration constants. Only the UI layer differs.

---

## Data Sources

| Source | GEE Collection | Use |
|---|---|---|
| VIIRS NRT 375m (SNPP) | `NASA/LANCE/SNPP_VIIRS/C2` | Active fire hotspot detections (2023-09-03 to present) |
| Landsat 8 C2 L2 SR | `LANDSAT/LC08/C02/T1_L2` | Pre/post composites, ATBI, dATBI, burn severity |
| Landsat 9 C2 L2 SR | `LANDSAT/LC09/C02/T1_L2` | Merged with L8 for ~8-day revisit |
| Sentinel-1 GRD IW | `COPERNICUS/S1_GRD` | SAR backscatter change (cloud-penetrating) |
| FAO GAUL 2024 Level 1 | `projects/sat-io/open-datasets/FAO/GAUL/GAUL_2024_L1` | Province boundaries (includes Kalimantan Utara) |
| Global Peatland Map 2.0 | Greifswald Mire Centre | Peatland overlay for contextual fire risk |

---

## Methodology

### Burn Index: ATBI / dATBI

The Automated Temporal Burn Index (Waleed & Bilal 2026) is used instead of the more common NBR/dNBR:

```
ATBI  = ((NIR - SWIR2) / (NIR + SWIR2)) * (SWIR2 / NIR)
dATBI = ATBI_pre - ATBI_post       (positive = burned)
```

Landsat bands: SR_B5 (NIR 865 nm) and SR_B7 (SWIR2 2200 nm). The multiplicative SWIR2/NIR term simultaneously amplifies NIR decrease and SWIR2 increase from combustion, yielding stronger spectral separation than dNBR in 13 of 15 wildfire test events across six continents (Waleed & Bilal 2026).

### Adaptive Otsu Thresholding

The burn threshold is computed server-side via `ee.Reducer.autoHistogram` and cumulative between-class variance maximization. A conservative bias (+0.02) and minimum floor (0.03) are applied to reduce commission errors. Fixed FIREMON dNBR thresholds (Key & Benson 2006) were designed for temperate conifer forests and require local calibration in tropical peatland contexts (Afira 2022).

### Burn Severity Classes

| Class | Threshold |
|---|---|
| Low | `dATBI >= T` |
| Moderate | `dATBI >= T + 0.30` |
| High | `dATBI >= T + 0.60` |

Where `T` is the scene-specific Otsu threshold.

### VIIRS Hotspots

Active fire detections from `NASA/LANCE/SNPP_VIIRS/C2`. The `confidence` band encodes `0 = low / 1 = nominal / 2 = high`. Low-confidence detections are excluded -- Urbanski (2018) demonstrates they substantially increase false positives. Nominal and high detections are max-composited for display and counted per province.

### SAR Backscatter Change

Sentinel-1 GRD IW descending orbit, VV + VH polarizations. Mean composites computed for pre-fire and post-fire windows. `dVV = post - pre` (dB): negative dVV (decreased backscatter) indicates canopy loss and surface exposure consistent with burning. Siegert (2000) established this approach for Kalimantan fires using ERS-2 C-band SAR.

### Landsat Processing

Cloud masking uses QA_PIXEL bits 3 (cloud) and 4 (cloud shadow). C2 L2 scale factor applied: `DN * 0.0000275 - 0.2`. Merged L8 + L9 median composites provide approximately 8-day revisit at 30 m. A water mask (NIR > 0.05) is applied inside ATBI computation to prevent SWIR2/NIR division blowup near open water.

---

## Configuration

All tuneable parameters are defined at the top of each script:

| Variable | Default | Purpose |
|---|---|---|
| `END_DATE` | `2026-08-21` | Post-fire window end |
| `START_DATE` | `2026-08-01` | Post-fire window start (21-day lookback) |
| `PRE_FIRE_START` | `2026-07-01` | Pre-fire reference composite start |
| `PRE_FIRE_END` | `2026-07-31` | Pre-fire reference composite end |
| `L_CLOUD_MAX` | `80` | Maximum scene-level cloud % for Landsat filter |
| `OTSU_BIAS` | `0.02` | Conservative upward bias added to raw Otsu threshold |
| `OTSU_MIN` | `0.03` | Minimum dATBI floor for Otsu threshold |
| `SEV_MOD_OFFSET` | `0.30` | Low/Moderate boundary = T + 0.30 |
| `SEV_HIGH_OFFSET` | `0.60` | Moderate/High boundary = T + 0.60 |
| `BURN_MIN_HA` | `10` | Minimum patch area for vectorized burn polygons |
| `VECTORIZE_SCALE` | `60` | Vectorization scale in metres (30 m risks compute timeout at regional scale) |

---

## Running the Scripts

Scripts run in the [GEE Code Editor](https://code.earthengine.google.com/). Paste the script contents or use the `earthengine` CLI:

```bash
earthengine run kalimantan_wildfire_monitoring_aug2026.js
```

**Desktop script**: requires an `aoi` geometry imported via the Code Editor "Imports" panel (polygon covering Kalimantan). For standalone GEE App deployment, replace the import with:

```javascript
var aoi = ee.Geometry.BBox(108.0, -4.1, 119.0, 4.2);
```

**Mobile script**: `aoi` is hardcoded as the BBox above -- ready for direct App deployment without imports.

---

## Layer Stack

The following layers are available in both scripts (desktop via checkboxes, mobile via tab bar):

1. Province Boundaries
2. Landsat Pre-fire True Color (July 2026)
3. Landsat Post-fire True Color (August 2026)
4. Landsat Post-fire False Color SWIR (SR_B7, SR_B5, SR_B4)
5. dATBI continuous (diagnostic)
6. Burn Severity classified (Low / Moderate / High)
7. VIIRS Hotspots -- Nominal confidence
8. VIIRS Hotspots -- High confidence
9. SAR Backscatter Change dVV
10. Peatland overlay (Global Peatland Map 2.0)

---

## UI Layout

**Desktop** (`kalimantan_wildfire_monitoring_aug2026.js`)
- Left panel: layer toggles, legends, Otsu threshold display, opacity sliders, pixel inspector
- Center: full-screen map
- Right panel: async statistics (total VIIRS count, per-province breakdown, burn area by severity class, methodology notes)

**Mobile** (`kalimantan_wildfire_monitoring_aug2026_mobile.js`)
- Full-screen map with no side panels
- Floating tab bar at bottom: Layers / Stats / Info
- Native GEE layer list retained for opacity control
- Map controls reduced (no map type selector, no fullscreen button) to maximize screen real estate on portrait smartphones

---

## References

- Afira, N., & Wijayanto, A. W. (2022). Mono-temporal and multi-temporal approaches for burnt area detection using Sentinel-2 satellite imagery (a case study of Rokan Hilir Regency, Indonesia). *Ecological Informatics, 69*, 101677. https://doi.org/10.1016/j.ecoinf.2022.101677

- Franceschini, G., Khan, A., Moretti, L., Nyabuti, K., Asif, M., Bezuidenhoudt, E., & Morteo, K. (2025). The Global Administrative Unit Layers (GAUL) 2024. Technical guidelines. Rome, FAO. https://doi.org/10.4060/cd4262en

- Giglio, L., Boschetti, L., Roy, D. P., Hall, J. V., Zubkova, M., Humber, M., Huang, H., & Oles, V. (2025). The NASA VIIRS burned area product, global validation, and intercomparison with the NASA MODIS burned area product. *Remote Sensing of Environment, 331*, 115006. https://doi.org/10.1016/j.rse.2025.115006

- Greifswald Mire Centre (2022). Global Peatland Map 2.0. Underlying dataset of the UNEP Global Peatland Assessment -- The State of the World's Peatlands: Evidence for action toward the conservation, restoration, and sustainable management of peatlands. Global Peatlands Initiative, United Nations Environment Programme, Nairobi. https://www.greifswald-moor-centrum.de/en/services/gis-data/global-peatland-map-2-0/

- Kadir, E. A., Rosa, S. L., Syukur, A., Othman, M., & Daud, H. (2021). Forest fire spreading and carbon concentration identification in tropical region Indonesia. *Alexandria Engineering Journal, 61*(2), 1551-1561. https://doi.org/10.1016/j.aej.2021.06.064

- Kurbanov, E., Vorobev, O., Lezhnin, S., Sha, J., Wang, J., Li, X., Cole, J., Dergunov, D., & Wang, Y. (2022). Remote sensing of forest burnt area, burn severity, and post-fire recovery: A review. *Remote Sensing, 14*(19), 4714. https://doi.org/10.3390/rs14194714

- Pinto, M. M., Trigo, R. M., Trigo, I. F., & DaCamara, C. C. (2021). A practical method for high-resolution burned area monitoring using Sentinel-2 and VIIRS. *Remote Sensing, 13*(9), 1608. https://doi.org/10.3390/rs13091608

- Siegert, F., & Hoffmann, A. A. (2000). The 1998 forest fires in East Kalimantan (Indonesia). *Remote Sensing of Environment, 72*(1), 64-77. https://doi.org/10.1016/s0034-4257(99)00092-9

- Urbanski, S., Nordgren, B., Albury, C., Schwert, B., Peterson, D., Quayle, B., & Hao, W. M. (2018). A VIIRS direct broadcast algorithm for rapid response mapping of wildfire burned area in the western United States. *Remote Sensing of Environment, 219*, 271-283. https://doi.org/10.1016/j.rse.2018.10.007

- Waleed, M., & Bilal, M. (2026). BAM: A physics-informed self-supervised framework for near-real-time wildfire burned area mapping from multi-source earth observation. *International Journal of Applied Earth Observation and Geoinformation, 153*, 105517. https://doi.org/10.1016/j.jag.2026.105517

---

## Author

**Muhammad Wahyu Ramadhan**
[GitHub](https://github.com/mwahyur46) · [LinkedIn](https://linkedin.com/in/mwahyur)
