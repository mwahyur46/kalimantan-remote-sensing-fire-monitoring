# Kalimantan Wildfire Monitoring

Satellite-based wildfire monitoring across all five Kalimantan provinces (Indonesian Borneo) using Google Earth Engine. Integrates VIIRS active-fire detections, Landsat 8/9 burn severity mapping, and Sentinel-1 SAR cloud-penetrating change detection for the August 2026 fire season.

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

| Citation | Contribution |
|---|---|
| Waleed & Bilal (2026) | ATBI index, dATBI, adaptive Otsu thresholding, BAM framework |
| Giglio et al. (2025) | NASA VIIRS burned area product VNP64A1 validation |
| Pinto et al. (2021) | Sentinel-2 + VIIRS practical burned area, GEE cloud-free compositing |
| Afira (2022) | Multi-index burned area mapping, Indonesian peatland context |
| Kurbanov et al. (2022) | Review of 329 RS burn severity studies (2000-2020) |
| Urbanski (2018) | VIIRS rapid-response burned area, confidence filtering |
| Siegert (2000) | ERS-2 SAR burn mapping, East Kalimantan 1998 El Nino fires |

---

## Author

**Muhammad Wahyu Ramadhan**
[GitHub](https://github.com/mwahyur46) · [LinkedIn](https://linkedin.com/in/mwahyur)
