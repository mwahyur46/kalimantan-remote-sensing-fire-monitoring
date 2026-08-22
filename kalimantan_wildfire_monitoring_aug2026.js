/* ============================================================================
 *  Kalimantan Wildfire Monitoring App -- August 2026
 *
 *  Monitors widespread wildfire events across Kalimantan (Indonesian Borneo)
 *  using multi-source satellite data:
 *    - VIIRS NRT 375m (SNPP) : active fire hotspot detection (NASA/LANCE/SNPP_VIIRS/C2)
 *    - Landsat 8/9 C2 L2 SR  : pre/post-fire composites, ATBI, dATBI burn severity
 *    - Sentinel-1 GRD        : SAR backscatter change (cloud-penetrating layer)
 *
 *  Burn index methodology:
 *    ATBI  = ((NIR-SWIR2)/(NIR+SWIR2)) * (SWIR2/NIR)   -- Waleed & Bilal (2026)
 *    dATBI = ATBIpre - ATBIpost  (positive = burned)
 *    Burn threshold from adaptive Otsu method (Waleed & Bilal 2026, server-side)
 *
 *  Methodology references:
 *    Waleed & Bilal (2026)    -- ATBI, dATBI, adaptive Otsu thresholding
 *    Pinto et al. (2021)      -- S2 + VIIRS practical burned area method
 *    Afira (2022)             -- Multi-temporal burn detection, Indonesia context
 *    Kurbanov et al. (2022)   -- Review: forest burn severity RS methods
 *    Giglio et al. (2025)     -- NASA VIIRS burned area product validation
 *    Urbanski (2018)          -- VIIRS rapid response burned area mapping
 *    Siegert (2000)           -- SAR-based fire monitoring, East Kalimantan
 *
 *  Author  : Muhammad Wahyu Ramadhan
 *  GitHub  : github.com/mwahyur46
 *  LinkedIn: linkedin.com/in/mwahyur
 * ============================================================================
 *  Required asset (imported via the Code Editor "Imports" panel):
 *    aoi -- ee.Geometry polygon covering the Kalimantan region of interest
 *
 *  For standalone GEE App deployment, replace the aoi import with:
 *    var aoi = ee.Geometry.BBox(108.0, -4.1, 119.0, 4.2);
 * ============================================================================ */

// ============================================================================
// 1. CONFIGURATION -- all tuneable parameters in one place
// ============================================================================
var END_DATE       = '2026-08-21';   // Post-fire window end (analysis date)
var START_DATE     = '2026-08-01';   // Post-fire window start (21-day lookback)
var PRE_FIRE_START = '2026-07-01';   // Pre-fire reference composite start
var PRE_FIRE_END   = '2026-07-31';   // Pre-fire reference composite end

// Landsat scene-level cloud cover filter (Landsat C2 metadata field)
var L_CLOUD_MAX    = 80;             // Drop scenes with cloud cover > 80%

// Adaptive Otsu threshold tuning (Waleed & Bilal 2026)
var OTSU_BIAS      = 0.02;           // Conservative upward bias added to raw Otsu T
var OTSU_MIN       = 0.03;           // Minimum dATBI floor (avoids near-zero thresholds)

// Severity class offsets relative to Otsu threshold T (Waleed & Bilal 2026)
var SEV_MOD_OFFSET = 0.30;           // Low/Moderate boundary = T + SEV_MOD_OFFSET
var SEV_HIGH_OFFSET = 0.60;          // Moderate/High boundary = T + SEV_HIGH_OFFSET

var BURN_MIN_HA    = 10;             // Minimum patch size for vectorized burn polygons

// NOTE: VECTORIZE_SCALE at 60 m balances spatial detail against GEE compute
// timeout risk for a region as large as Kalimantan. For smaller sub-AOIs,
// reduce to 30 m.
var VECTORIZE_SCALE = 60;

// Province names matching FAO GAUL 2024 Level 1 gaul1_name field.
var KALIMANTAN_ADM1 = [
  'Kalimantan Barat',
  'Kalimantan Tengah',
  'Kalimantan Selatan',
  'Kalimantan Timur',
  'Kalimantan Utara'
];

// Visual palette constants -- referenced by both layers and legends
var PALETTE_SEVERITY = ['#ffffb2', '#fd8d3c', '#e31a1c'];  // low / moderate / high
var PALETTE_HOTSPOT  = {nominal: '#f7dc6f', high: '#e74c3c'};
var PALETTE_SAR      = ['#d73027', '#f7f7f7', '#1a9641'];   // decrease / stable / increase
var PALETTE_DATBI    = ['#4575b4', '#f7f7f7', '#d73027'];   // negative / zero / positive
var PALETTE_PEAT     = ['#c7e9c0', '#238b45'];              // mosaic / peat-dominated

// ============================================================================
// 2. MAP INITIALIZATION
// ============================================================================
Map.setCenter(114, 0, 6);
Map.setOptions('HYBRID');
Map.style().set('cursor', 'crosshair');

// ============================================================================
// 3. PROVINCE BOUNDARIES
// ============================================================================
/**
 * Loads Kalimantan provincial boundaries from FAO GAUL 2024 Level 1.
 * GAUL 2024 includes Kalimantan Utara unlike the 2015 edition.
 * Citation: Franceschini et al. (2025) FAO. https://doi.org/10.4060/cd4262en
 *
 * @returns {ee.FeatureCollection} Five Kalimantan province polygons
 */
function loadProvinces() {
  return ee.FeatureCollection('projects/sat-io/open-datasets/FAO/GAUL/GAUL_2024_L1')
    .filter(ee.Filter.eq('iso3_code', 'IDN'))
    .filter(ee.Filter.inList('gaul1_name', KALIMANTAN_ADM1));
}

// ============================================================================
// 4. VIIRS NRT -- Active fire detections (375m, Suomi-NPP VIIRS C2)
// ============================================================================
/**
 * Loads VIIRS 375m NRT active fire data from NASA/LANCE/SNPP_VIIRS/C2.
 * This collection supersedes MODIS FIRMS for NRT monitoring: 375m resolution
 * detects smaller fires missed at MODIS 1km, and Suomi-NPP VIIRS is actively
 * operating in 2026 (Giglio et al. 2025).
 *
 * Key band: 'Mask' (fire detection categorisation per pixel)
 *   0  = non-fire land
 *   3  = non-fire water
 *   5  = fire, low confidence
 *   7  = fire, nominal confidence   <- included
 *   8  = fire, high confidence      <- included
 *   9  = fire, high confidence      <- included
 *
 * Additional bands: MaxFRP (fire radiative power, MW), QA, sample.
 *
 * Returns two masked images (nominal and high) plus ee.Number pixel counts.
 *
 * @param {ee.Geometry} geometry  - Spatial filter and clip geometry
 * @param {string}      startDate - Inclusive start date (YYYY-MM-DD)
 * @param {string}      endDate   - Inclusive end date (YYYY-MM-DD)
 * @returns {{nominal: ee.Image, high: ee.Image,
 *            nominalCount: ee.Number, highCount: ee.Number,
 *            totalCount: ee.Number, provinceCount: function}}
 */
function loadVIIRS(geometry, startDate, endDate) {
  var col = ee.ImageCollection('NASA/LANCE/SNPP_VIIRS/C2')
    .filterDate(startDate, endDate)
    .filterBounds(geometry)
    .select('confidence');

  print('VIIRS NRT daily image count (' + startDate + ' to ' + endDate + '):', col.size());

  // Max-composite: each 375m pixel gets the highest confidence seen across days.
  var maxConf = col.max().clip(geometry);

  // Diagnostic: print confidence value range to confirm encoding (0-100 integer expected)
  print('VIIRS confidence value range:', maxConf.reduceRegion({
    reducer: ee.Reducer.minMax(), geometry: geometry, scale: 10000, maxPixels: 1e9
  }));

  // VIIRS NRT confidence encoding (NASA/LANCE/SNPP_VIIRS/C2):
  //   0 = low confidence  (excluded per Urbanski 2018)
  //   1 = nominal confidence
  //   2 = high confidence
  var nominal = maxConf.updateMask(maxConf.eq(1));
  var high    = maxConf.updateMask(maxConf.eq(2));
  var allFire = maxConf.updateMask(maxConf.gte(1));

  // Pixel counts at 375m native resolution
  function pixelCount(img) {
    return img.reduceRegion({
      reducer  : ee.Reducer.count(),
      geometry : geometry,
      scale    : 375,
      maxPixels: 1e10,
      tileScale: 4
    }).getNumber('confidence');
  }

  var nomCount   = pixelCount(nominal);
  var highCount  = pixelCount(high);
  var totalCount = pixelCount(allFire);

  print('VIIRS fire pixels -- nominal:', nomCount, '| high:', highCount, '| total:', totalCount);

  // Province-level count function (called per province in panel builder)
  function provinceCount(provGeom) {
    return allFire.reduceRegion({
      reducer  : ee.Reducer.count(),
      geometry : provGeom,
      scale    : 375,
      maxPixels: 1e10,
      tileScale: 4
    }).getNumber('confidence');
  }

  return {
    nominal      : nominal,
    high         : high,
    nominalCount : nomCount,
    highCount    : highCount,
    totalCount   : totalCount,
    provinceCount: provinceCount
  };
}

// ============================================================================
// 5. LANDSAT 8/9 -- Cloud masking + median composite
// ============================================================================
/**
 * Masks cloud and cloud shadow pixels in a Landsat Collection 2 Level-2 image
 * using QA_PIXEL bitmask (bit 3 = cloud, bit 4 = cloud shadow). Applies
 * the C2 L2 scale factor: SR = DN * 0.0000275 - 0.2.
 * Retains bands SR_B2 (Blue), SR_B3 (Green), SR_B4 (Red),
 * SR_B5 (NIR 865nm), SR_B7 (SWIR2 2200nm).
 *
 * @param {ee.Image} image - Raw Landsat C2 L2 image
 * @returns {ee.Image} Cloud-masked, reflectance-scaled image
 */
function maskLandsatClouds(image) {
  var qa         = image.select('QA_PIXEL');
  var cloudMask  = qa.bitwiseAnd(1 << 3).eq(0);
  var shadowMask = qa.bitwiseAnd(1 << 4).eq(0);
  return image
    .updateMask(cloudMask.and(shadowMask))
    .select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
    .multiply(0.0000275).add(-0.2)
    .copyProperties(image, ['system:time_start']);
}

/**
 * Builds a cloud-masked Landsat 8+9 SR median composite for a given period.
 * Landsat 8 and 9 share identical band designations in C2 L2; merging them
 * gives an ~8-day effective revisit frequency over Kalimantan.
 *
 * @param {ee.Geometry} geometry  - Clip and filter geometry
 * @param {string}      startDate - Composite start date (YYYY-MM-DD)
 * @param {string}      endDate   - Composite end date (YYYY-MM-DD)
 * @returns {ee.Image} Median composite clipped to geometry
 */
function getLandsatComposite(geometry, startDate, endDate) {
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(geometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', L_CLOUD_MAX))
    .map(maskLandsatClouds);
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(geometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', L_CLOUD_MAX))
    .map(maskLandsatClouds);
  var merged = l8.merge(l9);
  print('Landsat 8+9 scene count (' + startDate + ' to ' + endDate + '):', merged.size());
  return merged.median().clip(geometry);
}

// ============================================================================
// 6. ATBI AND dATBI COMPUTATION
// ============================================================================
/**
 * Computes the Automated Temporal Burn Index (ATBI) for a Landsat image.
 *
 * ATBI = ((NIR - SWIR2) / (NIR + SWIR2)) * (SWIR2 / NIR)
 *
 * The first term is analogous to NBR; multiplying by SWIR2/NIR applies
 * multiplicative scaling that simultaneously amplifies NIR decrease and
 * SWIR2 increase from combustion, improving separation of burned from
 * unburned pixels relative to dNBR in 13 of 15 test events
 * (Waleed & Bilal 2026).
 *
 * Bands used: SR_B5 (NIR 865nm), SR_B7 (SWIR2 2200nm).
 *
 * @param {ee.Image} lsImage - Scaled Landsat SR composite
 * @returns {ee.Image} Single-band ATBI image
 */
function computeATBI(lsImage) {
  var nir   = lsImage.select('SR_B5');
  var swir2 = lsImage.select('SR_B7');
  // Valid land pixels: both bands must be positive and below 1.0 (physical SR range).
  // NIR < 0.1 covers open water (very low NIR), where SWIR2/NIR diverges.
  // Requiring NIR > 0.05 removes ocean, deep water, and cloud-shadow artefacts.
  var validMask = nir.gt(0.05)
    .and(swir2.gt(0))
    .and(nir.lt(1.0))
    .and(swir2.lt(1.0));
  var nbr = nir.subtract(swir2).divide(nir.add(swir2));
  return nbr.multiply(swir2.divide(nir))
    .updateMask(validMask)
    .rename('ATBI');
}

/**
 * Computes delta ATBI (dATBI) as the difference between pre-fire and post-fire
 * ATBI images. Positive dATBI indicates fire impact.
 *
 * dATBI = ATBIpre - ATBIpost
 *
 * @param {ee.Image} atbiPre  - Pre-fire ATBI image
 * @param {ee.Image} atbiPost - Post-fire ATBI image
 * @returns {ee.Image} Single-band dATBI image
 */
function computeDATBI(atbiPre, atbiPost) {
  return atbiPre.subtract(atbiPost).rename('dATBI');
}

// ============================================================================
// 7. ADAPTIVE OTSU THRESHOLDING
// ============================================================================
/**
 * Estimates an adaptive burn threshold from the dATBI histogram using
 * Otsu's method (maximise between-class variance) computed server-side
 * via ee.Array operations.
 *
 * The raw Otsu threshold is adjusted upward by OTSU_BIAS (conservative bias
 * to reduce commission errors) and clamped to a minimum of OTSU_MIN.
 * Histogrammed at 90 m scale (3x Landsat native) to reduce computation cost.
 *
 * Reference: Waleed & Bilal (2026) -- adaptive Otsu for dATBI.
 *
 * @param {ee.Image}    datbi    - dATBI image (single band 'dATBI')
 * @param {ee.Geometry} geometry - Region over which to compute histogram
 * @returns {ee.Number} Burn onset threshold T
 */
function otsuThreshold(datbi, geometry) {
  var hist = datbi.reduceRegion({
    reducer  : ee.Reducer.autoHistogram({maxBuckets: 256, cumulative: false}),
    geometry : geometry,
    scale    : 90,
    maxPixels: 1e10,
    tileScale: 4
  });

  var arr    = ee.Array(hist.get('dATBI'));
  var vals   = arr.slice(1, 0, 1).project([0]);
  var counts = arr.slice(1, 1, 2).project([0]);

  var n      = vals.length().get([0]);
  var total  = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sumAll = vals.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);

  var cumW   = counts.accum(0);
  var cumSum = vals.multiply(counts).accum(0);

  var totalArr  = ee.Array(ee.List.repeat(total,  n));
  var sumAllArr = ee.Array(ee.List.repeat(sumAll, n));

  var wB = cumW;
  var wF = totalArr.subtract(cumW);
  var mB = cumSum.divide(wB.add(1e-10));
  var mF = sumAllArr.subtract(cumSum).divide(wF.add(1e-10));

  var between = wB.multiply(wF).multiply(mB.subtract(mF).pow(2));

  // Find the maximum between-class variance, then locate the bin value where
  // between equals that max. Avoids argmax indexing entirely -- .get() on an
  // ee.Array returns a plain JS primitive when GEE can resolve it, which breaks
  // ee.Number() when the value is 0. Using max+eq keeps everything server-side.
  var maxB    = between.reduce(ee.Reducer.max(), [0]);        // shape-[1] ee.Array
  var maxBRep = maxB.repeat(0, ee.Number(n));                 // broadcast to length n
  var isMax   = between.eq(maxBRep);                          // 1 at argmax, 0 elsewhere
  var rawT    = vals.multiply(isMax).reduce(ee.Reducer.sum(), [0]).get([0]);

  // Apply bias and floor
  return rawT.add(OTSU_BIAS).max(OTSU_MIN);
}

// ============================================================================
// 8. BURN SEVERITY CLASSIFICATION
// ============================================================================
/**
 * Classifies a dATBI image into three burn severity categories relative to
 * the adaptive Otsu threshold T (Waleed & Bilal 2026):
 *   1 = Low      severity  [T,        T + SEV_MOD_OFFSET)
 *   2 = Moderate severity  [T + 0.30, T + SEV_HIGH_OFFSET)
 *   3 = High     severity  [T + 0.60, inf)
 * Pixels below T are masked as unburned.
 *
 * @param {ee.Image}  datbi          - dATBI image
 * @param {ee.Number} burnThreshold  - Otsu-derived threshold T
 * @returns {ee.Image} Integer image (1-3), masked where dATBI < T
 */
function classifyBurnSeverity(datbi, burnThreshold) {
  var T   = burnThreshold;
  var mod = T.add(SEV_MOD_OFFSET);
  var hi  = T.add(SEV_HIGH_OFFSET);

  // Start from datbi's own mask so cloud/water gaps are not rendered as gray.
  return datbi.multiply(0)
    .where(datbi.gte(T).and(datbi.lt(ee.Image(mod))),  1)
    .where(datbi.gte(ee.Image(mod)).and(datbi.lt(ee.Image(hi))), 2)
    .where(datbi.gte(ee.Image(hi)),                    3)
    .updateMask(datbi.gte(ee.Image(T)))
    .rename('burn_severity')
    .toInt();
}

// ============================================================================
// 9. SENTINEL-1 SAR -- Cloud-penetrating backscatter change
// ============================================================================
/**
 * Loads Sentinel-1 GRD IW descending-orbit VV+VH composites and computes
 * backscatter change as (post - pre) in dB. Fire-affected surfaces typically
 * show negative dVV due to loss of volume scatterers (canopy/vegetation).
 *
 * SAR's cloud-penetration capability makes it a critical complement to
 * optical data in persistently cloudy tropical regions, as demonstrated
 * by Siegert (2000) for the 1998 East Kalimantan fires.
 *
 * @param {ee.Geometry} geometry  - Area of interest
 * @param {string}      preStart  - Pre-fire start date
 * @param {string}      preEnd    - Pre-fire end date
 * @param {string}      postStart - Post-fire start date
 * @param {string}      postEnd   - Post-fire end date
 * @returns {{pre: ee.Image, post: ee.Image, change: ee.Image}}
 */
function loadSAR(geometry, preStart, preEnd, postStart, postEnd) {
  var s1Col = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(geometry)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .select(['VV', 'VH']);

  var s1Pre    = s1Col.filterDate(preStart, preEnd).mean().clip(geometry);
  var s1Post   = s1Col.filterDate(postStart, postEnd).mean().clip(geometry);
  var s1Change = s1Post.subtract(s1Pre).rename(['dVV', 'dVH']);

  print('S1 pre collection size  (' + preStart + ' to ' + preEnd + '):', s1Col.filterDate(preStart, preEnd).size());
  print('S1 post collection size (' + postStart + ' to ' + postEnd + '):', s1Col.filterDate(postStart, postEnd).size());

  return {pre: s1Pre, post: s1Post, change: s1Change};
}

// ============================================================================
// 10. BURN AREA CALCULATION PER SEVERITY CLASS
// ============================================================================
/**
 * Computes burned area in hectares per severity class using pixel-area
 * reduction. Pixel area (m2) is divided by 10 000 to yield hectares.
 * Runs at 30 m scale to match Landsat native resolution.
 *
 * @param {ee.Image}    burnSeverity - Classified severity image (int 1-3)
 * @param {ee.Geometry} geometry    - Reduce region geometry
 * @returns {{lowHa: ee.Number, modHa: ee.Number, highHa: ee.Number}}
 */
function computeBurnAreas(burnSeverity, geometry) {
  var pixelArea = ee.Image.pixelArea().divide(10000);

  function areaForClass(cls) {
    return pixelArea
      .updateMask(burnSeverity.eq(cls))
      .rename('area_ha')
      .reduceRegion({
        reducer  : ee.Reducer.sum(),
        geometry : geometry,
        scale    : 30,
        maxPixels: 1e10,
        tileScale: 4
      })
      .getNumber('area_ha');
  }

  return {
    lowHa : areaForClass(1),
    modHa : areaForClass(2),
    highHa: areaForClass(3)
  };
}

// ============================================================================
// 11. VECTORIZE BURN SCAR (optional -- may timeout on very large AOIs)
// ============================================================================
/**
 * Converts the burn severity raster into a vector FeatureCollection.
 * Uses VECTORIZE_SCALE (default 60 m) to avoid GEE computation timeout.
 * Filters output polygons to those >= BURN_MIN_HA in area to remove noise.
 *
 * @param {ee.Image}    burnSeverity - Classified severity image
 * @param {ee.Geometry} geometry    - AOI for reduceToVectors
 * @returns {ee.FeatureCollection} Burn polygon features with 'burn_severity' property
 */
function vectorizeBurnScars(burnSeverity, geometry) {
  var vectors = burnSeverity.reduceToVectors({
    geometry      : geometry,
    scale         : VECTORIZE_SCALE,
    geometryType  : 'polygon',
    eightConnected: true,
    labelProperty : 'burn_severity',
    maxPixels     : 1e10,
    tileScale     : 4
  });

  return vectors.filter(ee.Filter.gte('area', BURN_MIN_HA * 10000));
}

// ============================================================================
// 12. HELPER -- divider line widget (reused throughout panels)
// ============================================================================
function divider() {
  return ui.Label('', {
    height         : '1px',
    backgroundColor: '#ccc',
    margin         : '8px 0',
    stretch        : 'horizontal'
  });
}

// ============================================================================
// 13. BUILD LEFT SIDEBAR
// ============================================================================
/**
 * Constructs the left sidebar panel containing:
 *   - App title and metadata
 *   - Layer visibility checkboxes (auto-generated from Map.layers())
 *   - Burn severity legend (gradient thumbnail)
 *   - VIIRS confidence legend (colored swatches)
 *   - SAR backscatter change legend (gradient thumbnail)
 *   - Opacity sliders for the severity and SAR layers
 *   - Pixel inspector with click-to-sample handler
 *
 * Must be called AFTER all Map.addLayer() calls so that the checkbox loop
 * reflects the correct layer list.
 *
 * @param {ee.Image}     datbi           - dATBI image (for inspector)
 * @param {ee.Image}     burnSeverity    - Classified severity image (for inspector)
 * @param {ee.Image}     sarChange       - SAR dVV/dVH change image (for inspector)
 * @param {ui.Map.Layer} severityLayer   - Reference to severity map layer object
 * @param {ui.Map.Layer} sarLayer        - Reference to SAR change map layer object
 * @param {ee.Number}    otsuVal         - Computed Otsu threshold (for display)
 * @param {ee.Image}     peatMask        - Binary peatland mask (1 = peat present)
 * @param {ee.Image}     burnSeverityRaw - Unmasked severity image (before peat filter)
 * @param {ee.Image}     datbiRaw        - Unmasked dATBI image (before peat filter)
 * @param {Object}       viirsObj        - {nominal, high, ...} from loadVIIRS()
 * @returns {ui.Panel} Constructed left panel
 */
function buildLeftPanel(datbi, burnSeverity, sarChange, severityLayer, sarLayer, otsuVal,
                        peatMask, burnSeverityRaw, datbiRaw, viirsObj) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style : {width: '300px', padding: '12px', backgroundColor: 'white'}
  });

  // --- Title block ---
  panel.add(ui.Label('Kalimantan Wildfire Monitor',
    {fontWeight: 'bold', fontSize: '16px', margin: '0 0 2px 0'}));
  panel.add(ui.Label('August 1-21, 2026 -- Active Fires & Burn Scars',
    {fontSize: '12px', color: '#444', margin: '0 0 2px 0'}));
  panel.add(ui.Label('VIIRS NRT 375m  |  Landsat 8/9 dATBI  |  Sentinel-1 SAR',
    {fontSize: '11px', color: '#666', margin: '0 0 6px 0'}));
  panel.add(ui.Label('Muhammad Wahyu Ramadhan',
    {fontSize: '11px', margin: '0'}));
  panel.add(ui.Label('github.com/mwahyur46',
    {fontSize: '11px', color: '#1a73e8', margin: '0 0 4px 0'}));
  panel.add(divider());

  // --- Layer toggles with per-layer descriptions ---
  panel.add(ui.Label('Map Layers',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 2px 0'}));

  // Helper: checkbox + description pair for a layer by index
  function layerRow(idx, description) {
    var layer = Map.layers().get(idx);
    var cb = ui.Checkbox({
      label: layer.getName(),
      value: layer.getShown(),
      style: {fontSize: '11px', fontWeight: 'bold', margin: '3px 0 0 0'}
    });
    cb.onChange(function(checked) { layer.setShown(checked); });
    panel.add(cb);
    panel.add(ui.Label(description,
      {fontSize: '10px', color: '#666', margin: '0 0 2px 8px'}));
  }

  layerRow(0, 'Outlines the five Kalimantan provinces as a geographic reference.');

  panel.add(ui.Label('Landsat Optical Imagery',
    {fontWeight: 'bold', fontSize: '11px', color: '#555', margin: '4px 0 0 0'}));
  layerRow(1,
    'True-color image from July 2026 (pre-fire baseline). Appears green because ' +
    'Kalimantan is mostly intact tropical forest. Compare with the post-fire image.'
  );
  layerRow(2,
    'True-color image from August 2026. Green = intact forest; darker brownish ' +
    'patches = burn scars. White patches = cloud cover blocking the satellite view.'
  );
  layerRow(3,
    'False-color using shortwave infrared (SWIR). Active burn scars appear ' +
    'bright red-orange, making them easier to spot than in the true-color view.'
  );

  panel.add(ui.Label('Peatland Data',
    {fontWeight: 'bold', fontSize: '11px', color: '#555', margin: '4px 0 0 0'}));
  layerRow(4,
    'Global Peatland Map 2.0 (1 km). Light green = peat in soil mosaic; ' +
    'dark green = peat-dominated. Source: Global Peatlands Initiative / COP26.'
  );

  panel.add(ui.Label('Burn Scar Analysis',
    {fontWeight: 'bold', fontSize: '11px', color: '#555', margin: '4px 0 0 0'}));
  layerRow(5,
    'Diagnostic: shows only where the satellite burn signal is positive, before ' +
    'the severity threshold is applied. Yellow = weak signal; red = strong burn signal. ' +
    'Transparent areas have no burn signal (intact forest, water, or cloud shadow).'
  );
  layerRow(6,
    'Fire damage intensity in three levels. Low = partial scorching; ' +
    'High = intense burn with major vegetation loss. Clouds and unburned areas ' +
    'are transparent.'
  );

  panel.add(ui.Label('Active Fire Detections',
    {fontWeight: 'bold', fontSize: '11px', color: '#555', margin: '4px 0 0 0'}));
  layerRow(7,
    'NASA VIIRS 375m pixel with a nominal (likely) active fire signal ' +
    'detected during August 2026.'
  );
  layerRow(8,
    'Active fire detection with a strong, high-confidence thermal signal -- ' +
    'the most certain fire locations.'
  );

  panel.add(ui.Label('Radar (Cloud-Penetrating)',
    {fontWeight: 'bold', fontSize: '11px', color: '#555', margin: '4px 0 0 0'}));
  layerRow(9,
    'Sentinel-1 radar signal change between July and August. Unlike optical ' +
    'cameras, radar passes through clouds. Red = signal loss, suggesting forest ' +
    'canopy was lost. Speckle-filtered to reduce noise.'
  );

  panel.add(divider());

  // --- Peatland mask toggle ---
  panel.add(ui.Label('Peatland Mask',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'When enabled, burn severity, dATBI, and VIIRS hotspot layers are restricted ' +
    'to peatland areas only (Global Peatland Map 2.0, 1 km). Peatland fires carry ' +
    'the greatest carbon emission and longest recovery implications (Afira 2022).',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
  ));
  panel.add(ui.Label(
    'Note: the peatland layer is at 1 km resolution. Burn and hotspot layer ' +
    'edges within peat zones may show blocky 1 km boundaries as a result. ' +
    'SAR and Landsat layers are unaffected.',
    {fontSize: '10px', color: '#888', margin: '0 0 6px 0'}
  ));

  var peatToggle = ui.Checkbox({
    label: 'Restrict fire layers to peatland areas',
    value: false,
    style: {fontSize: '12px', fontWeight: 'bold', margin: '2px 0 4px 0'}
  });
  panel.add(peatToggle);

  // Layer indices after reordering:
  //   5 = dATBI diagnostic, 6 = burn severity (severityLayer ref),
  //   7 = VIIRS nominal,    8 = VIIRS high confidence
  peatToggle.onChange(function(checked) {
    if (checked) {
      severityLayer.setEeObject(burnSeverityRaw.updateMask(peatMask));
      Map.layers().get(5).setEeObject(datbiRaw.updateMask(datbiRaw.gt(0)).updateMask(peatMask));
      Map.layers().get(7).setEeObject(viirsObj.nominal.updateMask(peatMask));
      Map.layers().get(8).setEeObject(viirsObj.high.updateMask(peatMask));
    } else {
      severityLayer.setEeObject(burnSeverityRaw);
      Map.layers().get(5).setEeObject(datbiRaw.updateMask(datbiRaw.gt(0)));
      Map.layers().get(7).setEeObject(viirsObj.nominal);
      Map.layers().get(8).setEeObject(viirsObj.high);
    }
  });

  panel.add(divider());

  // --- Burn severity legend ---
  panel.add(ui.Label('Burn Severity (Landsat dATBI)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'Areas where post-fire satellite reflectance changed significantly ' +
    'relative to the pre-fire baseline, indicating vegetation loss from burning.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}));
  panel.add(ui.Thumbnail({
    image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
               .visualize({min: 0, max: 1, palette: PALETTE_SEVERITY}),
    params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
    style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
  }));
  panel.add(ui.Panel([
    ui.Label('Low',      {fontSize: '10px', margin: '2px 0'}),
    ui.Label('Moderate', {fontSize: '10px', margin: '2px 0',
      stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('High',     {fontSize: '10px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

  // Async display of computed Otsu threshold
  var otsuLabel = ui.Label('Burn threshold T: computing...',
    {fontSize: '10px', color: '#888', margin: '2px 0'});
  panel.add(otsuLabel);
  otsuVal.evaluate(function(v) {
    otsuLabel.setValue('Burn threshold T: ' + (v != null ? v.toFixed(4) : 'N/A') +
      '  (scene-adaptive Otsu)');
  });

  panel.add(divider());

  // --- VIIRS confidence legend (colored swatches) ---
  panel.add(ui.Label('VIIRS Active Fire Hotspots (375m)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'Each colored cell is a 375 m x 375 m area where a satellite ' +
    'thermal anomaly consistent with an active fire was detected.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}));

  function swatch(color, text) {
    return ui.Panel([
      ui.Label('', {
        backgroundColor: color, width: '14px', height: '14px',
        margin: '2px 6px 2px 0', padding: '0', border: '1px solid #ccc'
      }),
      ui.Label(text, {fontSize: '11px', margin: '2px 0'})
    ], ui.Panel.Layout.flow('horizontal'), {margin: '1px 0'});
  }
  panel.add(swatch(PALETTE_HOTSPOT.nominal, 'Nominal confidence -- likely fire'));
  panel.add(swatch(PALETTE_HOTSPOT.high,    'High confidence -- strong fire signal'));

  panel.add(divider());

  // --- SAR backscatter change legend ---
  panel.add(ui.Label('SAR Radar Change (Sentinel-1)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'Radar signal change between July and August. Unlike optical sensors, ' +
    'radar penetrates cloud cover. Red areas lost vegetation structure, ' +
    'which can indicate fire damage.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}));
  panel.add(ui.Thumbnail({
    image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
               .visualize({min: 0, max: 1, palette: PALETTE_SAR}),
    params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
    style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
  }));
  panel.add(ui.Panel([
    ui.Label('Loss (red)',    {fontSize: '10px', margin: '2px 0'}),
    ui.Label('No change',    {fontSize: '10px', margin: '2px 0',
      stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('Gain (green)', {fontSize: '10px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

  panel.add(divider());

  // --- Opacity sliders ---
  panel.add(ui.Label('Burn Severity Opacity',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
  panel.add(ui.Slider({
    min: 0, max: 1, value: 0.85, step: 0.05,
    onChange: function(v) { severityLayer.setOpacity(v); },
    style: {stretch: 'horizontal'}
  }));

  panel.add(ui.Label('SAR Change Opacity',
    {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));
  panel.add(ui.Slider({
    min: 0, max: 1, value: 0.7, step: 0.05,
    onChange: function(v) { sarLayer.setOpacity(v); },
    style: {stretch: 'horizontal'}
  }));

  panel.add(divider());

  // --- Pixel inspector ---
  panel.add(ui.Label('Pixel Inspector',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));

  var inspContent = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'), style: {margin: '0'}
  });
  inspContent.add(ui.Label('Click anywhere on the map.',
    {fontSize: '11px', color: '#888'}));
  panel.add(inspContent);

  var SEVERITY_LABELS = ['', 'Low', 'Moderate', 'High'];

  Map.onClick(function(coords) {
    inspContent.clear();
    inspContent.add(ui.Label(
      'Lon: ' + coords.lon.toFixed(5) + '  Lat: ' + coords.lat.toFixed(5),
      {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
    ));
    inspContent.add(ui.Label('Sampling...', {fontSize: '11px', color: '#888'}));

    var point     = ee.Geometry.Point([coords.lon, coords.lat]);
    var sampleImg = datbi
      .addBands(burnSeverity)
      .addBands(sarChange.select('dVV'));

    sampleImg.reduceRegion({
      reducer: ee.Reducer.first(), geometry: point, scale: 30
    }).evaluate(function(vals) {
      inspContent.clear();
      inspContent.add(ui.Label(
        'Lon: ' + coords.lon.toFixed(5) + '  Lat: ' + coords.lat.toFixed(5),
        {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
      ));
      if (!vals || vals['dATBI'] === null || vals['dATBI'] === undefined) {
        inspContent.add(ui.Label('No optical data at this location.',
          {fontSize: '11px', color: '#c00'}));
        return;
      }
      var sevCode  = vals['burn_severity'];
      var sevLabel = (sevCode === null || sevCode === undefined)
        ? 'Unburned' : SEVERITY_LABELS[sevCode];
      var sevColor = (sevCode === null || sevCode === undefined)
        ? '#27ae60' : ['', '#fd8d3c', '#e67e22', '#e31a1c'][sevCode];

      inspContent.add(ui.Label(
        'dATBI: ' + vals['dATBI'].toFixed(4),
        {fontSize: '12px', fontWeight: 'bold', color: '#c0392b', margin: '0 0 2px 0'}
      ));
      inspContent.add(ui.Label(
        'Severity: ' + sevLabel,
        {fontSize: '12px', fontWeight: 'bold', color: sevColor, margin: '0 0 2px 0'}
      ));
      if (vals['dVV'] !== null && vals['dVV'] !== undefined) {
        inspContent.add(ui.Label(
          'SAR dVV: ' + vals['dVV'].toFixed(2) + ' dB',
          {fontSize: '12px', fontWeight: 'bold', color: '#1a73e8', margin: '0 0 2px 0'}
        ));
      }
    });
  });

  return panel;
}

// ============================================================================
// 14. BUILD RIGHT SIDEBAR (statistics panel)
// ============================================================================
/**
 * Constructs the right statistics panel showing:
 *   - Analysis date range summary
 *   - Total FIRMS fire pixel count
 *   - Fire pixel breakdown by Kalimantan province
 *   - Estimated burned area (ha) per severity class
 *   - Total estimated burned area
 *   - SAR supplementary note
 *
 * All statistics are computed server-side and populated asynchronously via
 * .evaluate() callbacks, identical to the pattern used in the AGB/canopy
 * height script for regression metrics.
 *
 * @param {Object}               viirs      - {totalCount, provinceCount, ...} from loadVIIRS()
 * @param {ee.FeatureCollection} provinces  - Kalimantan province polygons
 * @param {Object}               burnAreas  - {lowHa, modHa, highHa} ee.Numbers
 * @param {ee.Number}            otsuVal    - Computed Otsu threshold
 * @returns {ui.Panel} Constructed right panel
 */
function buildRightPanel(viirs, provinces, burnAreas, otsuVal) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style : {width: '300px', padding: '12px', backgroundColor: 'white'}
  });

  // --- Header ---
  panel.add(ui.Label('Fire Statistics',
    {fontWeight: 'bold', fontSize: '17px', margin: '0 0 2px 0'}));
  panel.add(ui.Label('Kalimantan -- August 2026 rapid fire assessment',
    {fontSize: '11px', color: '#666', margin: '0 0 4px 0'}));
  panel.add(ui.Label(
    'In August 2026, widespread fires were detected across Kalimantan ' +
    '(Indonesian Borneo) during the annual dry season. Fires affected peatland ' +
    'and forest areas across multiple provinces. The causes are under ' +
    'investigation by relevant authorities. This map combines three satellite ' +
    'systems to track where fires occurred, how severely vegetation was burned, ' +
    'and where radar data reveals fire signals hidden under cloud cover.',
    {fontSize: '10px', color: '#333', margin: '4px 0 6px 0'}
  ));
  panel.add(divider());

  // --- Date range ---
  panel.add(ui.Label('Analysis Period',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));
  panel.add(ui.Label('Post-fire : ' + START_DATE + ' to ' + END_DATE,
    {fontSize: '11px', margin: '1px 0'}));
  panel.add(ui.Label('Pre-fire reference: ' + PRE_FIRE_START + ' to ' + PRE_FIRE_END,
    {fontSize: '11px', margin: '1px 0'}));
  panel.add(divider());

  // --- VIIRS total fire pixel count ---
  panel.add(ui.Label('Active Fire Detections (VIIRS 375m)',
    {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 2px 0', color: '#222'}));
  panel.add(ui.Label(
    'Counted from NASA VIIRS satellite. Each fire pixel represents a ' +
    '375 m x 375 m ground area (about 14 hectares) where an active fire ' +
    'thermal signal was detected. Nominal and high confidence only.',
    {fontSize: '10px', color: '#555', margin: '0 0 6px 0'}));

  var totalLabel = ui.Label('Total fire pixels: computing...',
    {fontSize: '13px', margin: '2px 0', fontWeight: 'bold', color: '#e74c3c'});
  panel.add(totalLabel);
  viirs.totalCount.evaluate(function(n) {
    totalLabel.setValue('Total fire pixels: ' + (n != null ? n.toLocaleString() : '0'));
  });

  // --- Fire pixel count by province ---
  panel.add(ui.Label('By Province:',
    {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));

  KALIMANTAN_ADM1.forEach(function(name) {
    var shortName = name.replace('Kalimantan ', 'Kal. ');
    var lbl = ui.Label(shortName + ': computing...',
      {fontSize: '11px', margin: '1px 0'});
    panel.add(lbl);

    var filtered = provinces.filter(ee.Filter.eq('gaul1_name', name));
    var provGeom = ee.Geometry(ee.Algorithms.If(
      filtered.size().gt(0),
      filtered.first().geometry(),
      ee.Geometry.Point([0, 0])
    ));
    viirs.provinceCount(provGeom).evaluate(function(n) {
      lbl.setValue(shortName + ': ' + (n != null ? n.toLocaleString() : '0') + ' fire pixels');
    });
  });

  panel.add(divider());

  // --- Burn scar area by severity ---
  panel.add(ui.Label('Estimated Burned Area (Landsat)',
    {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 2px 0', color: '#222'}));
  panel.add(ui.Label(
    'Mapped from Landsat 8/9 satellite imagery by comparing vegetation ' +
    'reflectance before and after the fire using the ATBI burn index ' +
    '(Waleed & Bilal, 2026), which outperforms standard indices in 13 of ' +
    '15 wildfire test events globally. Areas covered by clouds are not ' +
    'mapped and will appear as gaps.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}));

  // Show the computed Otsu threshold value asynchronously
  var otsuStatLabel = ui.Label('Detection threshold: computing...',
    {fontSize: '10px', color: '#888', margin: '0 0 6px 0'});
  panel.add(otsuStatLabel);
  otsuVal.evaluate(function(v) {
    otsuStatLabel.setValue('Scene-adaptive threshold T = ' + (v != null ? v.toFixed(4) : 'N/A'));
  });

  var lowLabel  = ui.Label('Low severity: computing...',
    {fontSize: '12px', margin: '2px 0', color: '#e67e22'});
  var modLabel  = ui.Label('Moderate severity: computing...',
    {fontSize: '12px', margin: '2px 0', color: '#e74c3c'});
  var highLabel = ui.Label('High severity: computing...',
    {fontSize: '12px', margin: '2px 0', color: '#922b21'});
  var totLabel  = ui.Label('Total burned area: computing...',
    {fontSize: '13px', margin: '6px 0 2px 0', fontWeight: 'bold', color: '#c0392b'});

  panel.add(lowLabel);
  panel.add(modLabel);
  panel.add(highLabel);
  panel.add(totLabel);

  burnAreas.lowHa.evaluate(function(v) {
    lowLabel.setValue('Low severity (partial burn): ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.modHa.evaluate(function(v) {
    modLabel.setValue('Moderate severity: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.highHa.evaluate(function(v) {
    highLabel.setValue('High severity (intense burn): ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });

  burnAreas.lowHa.add(burnAreas.modHa).add(burnAreas.highHa).evaluate(function(v) {
    totLabel.setValue('Total burned area: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });

  panel.add(ui.Label(
    'Note: cloud-covered areas are excluded. Actual burned extent may be larger.',
    {fontSize: '10px', color: '#c0392b', margin: '4px 0'}
  ));

  panel.add(divider());

  // --- SAR note ---
  panel.add(ui.Label('Radar Imagery (Sentinel-1 SAR)',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));
  panel.add(ui.Label(
    'Radar signals pass through clouds, providing fire-related information ' +
    'where optical satellites cannot. A significant drop in radar signal ' +
    '(shown in red) after the fire period may indicate loss of forest canopy ' +
    'structure -- a technique demonstrated for Kalimantan fires by Siegert (2000). ' +
    'Use this layer to look for fire signals in cloud-covered areas. ' +
    'Cross-reference with the burn severity layer for confirmation.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
  ));

  panel.add(divider());

  // --- Limitations ---
  panel.add(ui.Label('Rapid Mapping Limitations',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0', color: '#922b21'}));
  panel.add(ui.Label(
    'This is a near-real-time rapid mapping product, not a validated fire ' +
    'damage assessment. Key limitations to be aware of:',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
  ));
  panel.add(ui.Label(
    '1. Cloud gaps: Kalimantan in August has persistent cloud cover. ' +
    'Burned areas under clouds are missed entirely.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));
  panel.add(ui.Label(
    '2. Timing lag: VIIRS fire data and Landsat imagery each have a 1 to 2 day ' +
    'ingestion delay. Very recent fires may not yet appear.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));
  panel.add(ui.Label(
    '3. Resolution: each mapped burn pixel covers 30 m x 30 m (0.09 ha). ' +
    'Small isolated fires below this size are not resolved.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));
  panel.add(ui.Label(
    '4. False positives: bare agricultural land and recently cleared areas ' +
    'can produce similar spectral signals to burn scars.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));

  panel.add(divider());

  // --- Future improvements ---
  panel.add(ui.Label('Potential Improvements',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0', color: '#1a5276'}));
  panel.add(ui.Label(
    '1. SAR-optical fusion: combining radar and optical burn signals ' +
    'would fill cloud gaps and reduce false positives (partially addressed ' +
    'by the SAR layer).',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));
  panel.add(ui.Label(
    '2. Carbon emission estimate: combining mapped burn area with the ' +
    'companion mangrove biomass analysis would allow a first-order ' +
    'estimate of CO₂ released.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));
  panel.add(ui.Label(
    '3. Time-series monitoring: tracking VIIRS fire counts day by day ' +
    'through the dry season would show how the event evolved over time.',
    {fontSize: '10px', color: '#555', margin: '0 0 3px 8px'}
  ));

  panel.add(divider());

  // --- Data & Methods ---
  panel.add(ui.Label('Data Sources',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));
  [
    'VIIRS NRT: NASA/LANCE/SNPP_VIIRS/C2 (375m)',
    'Landsat 8/9: USGS Collection 2 Level-2 SR (30m)',
    'Sentinel-1: ESA Copernicus GRD IW (10m)',
    'Boundaries: FAO GAUL 2024 Level 1 (Franceschini et al. 2025)',
    'Peatlands: Global Peatland Map 2.0 -- Global Peatlands Initiative / COP26 (1 km)'
  ].forEach(function(s) {
    panel.add(ui.Label(s, {fontSize: '10px', color: '#444', margin: '1px 0 1px 4px'}));
  });

  panel.add(ui.Label('Key References',
    {fontWeight: 'bold', fontSize: '13px', margin: '8px 0 4px 0'}));
  [
    {text: 'Giglio, L., Boschetti, L., Roy, D. P., Hall, J. V., Zubkova, M., Humber, M., Huang, H., & Oles, V. (2025). The NASA VIIRS burned area product, global validation, and intercomparison with the NASA MODIS burned area product. Remote Sensing of Environment, 331, 115006.',
     doi:  'https://doi.org/10.1016/j.rse.2025.115006'},
    {text: 'Kadir, E. A., Rosa, S. L., Syukur, A., Othman, M., & Daud, H. (2021). Forest fire spreading and carbon concentration identification in tropical region Indonesia. Alexandria Engineering Journal, 61(2), 1551-1561.',
     doi:  'https://doi.org/10.1016/j.aej.2021.06.064'},
    {text: 'Kurbanov, E., Vorobev, O., Lezhnin, S., Sha, J., Wang, J., Li, X., Cole, J., Dergunov, D., & Wang, Y. (2022). Remote Sensing of forest burnt area, burn Severity, and Post-Fire Recovery: A review. Remote Sensing, 14(19), 4714.',
     doi:  'https://doi.org/10.3390/rs14194714'},
    {text: 'Pinto, M. M., Trigo, R. M., Trigo, I. F., & DaCamara, C. C. (2021). A practical method for High-Resolution burned Area monitoring using Sentinel-2 and VIIRS. Remote Sensing, 13(9), 1608.',
     doi:  'https://doi.org/10.3390/rs13091608'},
    {text: 'Siegert, F., & Hoffmann, A. A. (2000). The 1998 forest fires in East Kalimantan (Indonesia). Remote Sensing of Environment, 72(1), 64-77.',
     doi:  'https://doi.org/10.1016/s0034-4257(99)00092-9'},
    {text: 'Urbanski, S., Nordgren, B., Albury, C., Schwert, B., Peterson, D., Quayle, B., & Hao, W. M. (2018). A VIIRS direct broadcast algorithm for rapid response mapping of wildfire burned area in the western United States. Remote Sensing of Environment, 219, 271-283.',
     doi:  'https://doi.org/10.1016/j.rse.2018.10.007'},
    {text: 'Waleed, M., & Bilal, M. (2026). BAM: A physics-informed self-supervised framework for near-real-time wildfire burned area mapping from multi-source earth observation. International Journal of Applied Earth Observation and Geoinformation, 153, 105517.',
     doi:  'https://doi.org/10.1016/j.jag.2026.105517'},
    {text: 'Greifswald Mire Centre (2022). Global Peatland Map 2.0. Underlying dataset of the UNEP Global Peatland Assessment -- The State of the World\'s Peatlands: Evidence for action toward the conservation, restoration, and sustainable management of peatlands, Global Peatlands Initiative, United Nations Environment Programme, Nairobi.',
     doi:  'https://www.greifswald-moor-centrum.de/en/services/gis-data/global-peatland-map-2-0/'},
    {text: 'Franceschini, G., Khan, A., Moretti, L., Nyabuti, K., Asif, M., Bezuidenhoudt, E., & Morteo, K. (2025). The Global Administrative Unit Layers (GAUL) 2024. Technical guidelines. Rome, FAO.',
     doi:  'https://doi.org/10.4060/cd4262en'}
  ].forEach(function(ref) {
    panel.add(ui.Label(ref.text,
      {fontSize: '10px', color: '#333', margin: '4px 0 0 0', fontWeight: 'bold'}));
    panel.add(ui.Label(ref.doi,
      {fontSize: '9px', color: '#888', margin: '0 0 0 4px', fontFamily: 'monospace'}));
  });

  return panel;
}

// ============================================================================
// 15. MAIN -- load data, compute products, add layers, build UI
// ============================================================================

// --- Province boundaries ---
var provinces = loadProvinces();

// --- Global Peatland Map 2.0 (1 km resolution) ---
// Pixel values: 1 = peat-dominated, 2 = peat in soil mosaic (Global Peatlands Initiative / COP26)
// Clipped to Kalimantan bounding box; unmasked so gte(1) picks up both classes.
var peatRaw  = ee.Image('projects/sat-io/open-datasets/ML-GLOBAL-PEATLAND-EXTENT')
                 .clip(aoi)
                 .unmask(0);
var peatMask = peatRaw.gte(1);   // binary: 1 where any peat class present

// --- VIIRS NRT active fire (375m raster) ---
var viirs = loadVIIRS(aoi, START_DATE, END_DATE);

// --- Landsat 8+9 composites ---
var lsPre  = getLandsatComposite(aoi, PRE_FIRE_START, PRE_FIRE_END);
var lsPost = getLandsatComposite(aoi, START_DATE,     END_DATE);

// --- ATBI and dATBI ---
var atbiPre  = computeATBI(lsPre);
var atbiPost = computeATBI(lsPost);
var datbi    = computeDATBI(atbiPre, atbiPost);

// --- Adaptive Otsu threshold ---
var otsuVal = otsuThreshold(datbi, aoi);

// --- Burn severity ---
var burnSeverity = classifyBurnSeverity(datbi, otsuVal);

// --- Burn area statistics (pixel-area reduction) ---
var burnAreas = computeBurnAreas(burnSeverity, aoi);

// --- Sentinel-1 SAR ---
var sar         = loadSAR(aoi, PRE_FIRE_START, PRE_FIRE_END, START_DATE, END_DATE);
// 60m focal median reduces SAR speckle while preserving real fire-patch signals.
var sarChangeSm = sar.change.focal_median(60, 'circle', 'meters');

// ============================================================================
// ADD MAP LAYERS
// Layer order here drives the checkbox order in the left panel.
// ============================================================================

// Province boundaries (outline only)
Map.addLayer(
  provinces.style({color: '#ffffff', fillColor: '00000000', width: 1.5}),
  {}, 'Province Boundaries', true
);

// Landsat pre-fire true color (reference baseline)
Map.addLayer(
  lsPre,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0, max: 0.25, gamma: 0.9},
  'Landsat Pre-fire True Color (Jul 2026)', false
);

// Landsat post-fire true color (primary optical layer)
Map.addLayer(
  lsPost,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0, max: 0.25, gamma: 0.9},
  'Landsat Post-fire True Color (Aug 2026)', true
);

// Landsat post-fire false color SWIR -- highlights active burn scars and bare soil
// RGB: SR_B7 (SWIR2), SR_B5 (NIR), SR_B4 (Red)
Map.addLayer(
  lsPost,
  {bands: ['SR_B7','SR_B5','SR_B4'], min: 0, max: 0.4, gamma: 0.9},
  'Landsat Post-fire False Color SWIR', false
);

// Peatland extent (Global Peatland Map 2.0, 1 km) -- index 4
// Placed above Landsat base imagery but below analysis layers so it serves as
// a contextual reference without obscuring burn products or hotspot detections.
// Value 1 = peat-dominated (dark green), value 2 = peat in soil mosaic (light green).
var peatLayer = ui.Map.Layer(
  peatRaw.updateMask(peatMask)
         .visualize({min: 1, max: 2, palette: PALETTE_PEAT, opacity: 0.35}),
  {}, 'Peatland Extent (Global Peatland Map 2.0)', false
);
Map.layers().add(peatLayer);

// dATBI continuous (diagnostic layer) -- index 5
// dATBI diagnostic: mask out negative values (cloud shadow / water artefacts)
// so only the positive burn signal is shown. Yellow-to-red palette avoids
// the misleading blue-dominant display from cloud shadow in the pre-fire composite.
Map.addLayer(
  datbi.updateMask(datbi.gt(0)),
  {min: 0, max: 0.4, palette: PALETTE_SEVERITY},
  'dATBI (burn signal only)', false
);

// Burn severity -- classified (primary output) -- index 6
// Added via ui.Map.Layer to retain a reference for the opacity slider
var severityLayer = ui.Map.Layer(
  burnSeverity,
  {min: 1, max: 3, palette: PALETTE_SEVERITY},
  'Burn Severity (classified)', true, 0.85
);
Map.layers().add(severityLayer);

// VIIRS 375m hotspots -- nominal confidence (confidence == 1) -- index 7
Map.addLayer(
  viirs.nominal,
  {min: 1, max: 1, palette: [PALETTE_HOTSPOT.nominal]},
  'VIIRS Hotspots -- Nominal (conf=1)', true
);

// VIIRS 375m hotspots -- high confidence (confidence == 2) -- index 8
Map.addLayer(
  viirs.high,
  {min: 2, max: 2, palette: [PALETTE_HOTSPOT.high]},
  'VIIRS Hotspots -- High Confidence (conf=2)', true
);

// SAR backscatter change dVV: speckle-filtered + ±3 dB range so fire-related
// canopy loss (-1 to -3 dB) renders as visible red rather than near-white.
// index 9
var sarLayer = ui.Map.Layer(
  sarChangeSm.select('dVV'),
  {min: -3, max: 3, palette: PALETTE_SAR},
  'SAR Backscatter Change dVV (optional)', false, 0.7
);
Map.layers().add(sarLayer);

// ============================================================================
// BUILD AND MOUNT UI PANELS
// ============================================================================

// Left panel must be built AFTER all Map.addLayer() calls
var leftPanel  = buildLeftPanel(datbi, burnSeverity, sarChangeSm, severityLayer, sarLayer, otsuVal,
                                peatMask, burnSeverity, datbi, viirs);
var rightPanel = buildRightPanel(viirs, provinces, burnAreas, otsuVal);

ui.root.insert(0, leftPanel);
ui.root.add(rightPanel);

// ============================================================================
// END OF SCRIPT
// ============================================================================

/* ============================================================================
 * KNOWN LIMITATIONS
 * ============================================================================
 *
 * 1. Cloud contamination in optical composites
 *    Kalimantan in August typically has persistent cloud cover, especially
 *    over peatland areas. The 14-day post-fire window and 80% cloud threshold
 *    may still result in cloud-gaps in the dATBI, causing burn area
 *    underestimation. Extend the post-fire window or lower L_CLOUD_MAX if
 *    the composite is sparse.
 *
 * 2. Otsu threshold sensitivity
 *    The server-side Otsu method works best when burned pixels constitute
 *    a meaningful fraction of the histogram. In scenes with very little
 *    burning, the bimodal assumption may not hold and T will be unreliable.
 *    Print 'dATBI min/max' and the Otsu T value to verify plausibility.
 *    The OTSU_BIAS and OTSU_MIN constants provide conservative fallbacks.
 *
 * 3. VIIRS NRT collection availability
 *    NASA/LANCE/SNPP_VIIRS/C2 is NRT (near-real-time) data with a ~6-24h
 *    ingestion lag. Standard science-quality VIIRS fire data (VNP14IMG)
 *    is available with a longer delay but with improved geolocation and
 *    calibration. For retrospective analysis, prefer VNP14IMG over NRT.
 *
 * 4. SAR orbit gaps
 *    Sentinel-1 descending orbit coverage over Kalimantan is not daily.
 *    If the post-fire window is short (< 12 days), there may be few or no
 *    S1 acquisitions. Check the print() output for collection sizes.
 *    Switch to ascending orbit or combine both passes if needed.
 *
 * 5. dATBI sensitivity in tropical peatlands
 *    ATBI was originally validated across 15 wildfire events globally
 *    (Waleed & Bilal 2026). Performance for sub-surface peatland fires
 *    (where optical signal is partially obscured) has not been independently
 *    validated. SAR layer provides a cloud-penetrating complement.
 *
 * 6. Vectorization timeout
 *    vectorizeBurnScars() is defined but NOT called by default due to
 *    GEE compute timeout risk at Kalimantan scale. Enable it for sub-AOIs
 *    or after narrowing the analysis to specific hotspot clusters.
 *
 * 7. Area statistics exclude patch-size filter
 *    computeBurnAreas() uses pixel-area reduction without patch-size
 *    filtering. Small salt-and-pepper noise pixels are included. The
 *    BURN_MIN_HA threshold is only applied in vectorizeBurnScars().
 *
 * ============================================================================
 * SUGGESTED IMPROVEMENTS
 * ============================================================================
 *
 * 1. VIIRS-seeded burned area (Pinto et al., 2021)
 *    Use VIIRS hotspot locations as seeds for region growing on the dATBI
 *    image. This links optically mapped burn scars to confirmed active fires,
 *    reducing commission errors from agricultural burning or bare soil.
 *
 * 2. Full BAM framework (Waleed & Bilal, 2026)
 *    Add the GTB (Gradient Tree Boost) refinement step: generate pseudo-labels
 *    from dATBI+Otsu, train GTB with GLCM texture features and FABDEM terrain
 *    variables, and refine the burn mask. Improves F1 from ~0.85 to ~0.99.
 *
 * 3. Peatland mask overlay (implemented)
 *    Global Peatland Map 2.0 (projects/sat-io/open-datasets/ML-GLOBAL-PEATLAND-EXTENT)
 *    is loaded and exposed as a left-panel checkbox toggle. When enabled, burn
 *    severity and dATBI layers are restricted to peatland extent. The peatland
 *    fill layer (layer index 9) is also available as a separate map layer toggle.
 *
 * 4. NBR from Sentinel-1 (SAR-NBR)
 *    Compute a SAR-based burn index (e.g., VH/VV ratio change) following
 *    Siegert (2000) as a cloud-independent index proxy. Fuse with optical
 *    dATBI for cloud-gap filling.
 *
 * 5. Time-series animation
 *    Build a VIIRS hotspot density time series chart (ee.ImageCollection
 *    daily composites) to show the temporal evolution of the fire event.
 *
 * 6. Emissions estimation
 *    Combine burn area with fuel load data (AGB map from the companion
 *    agb_canopy_height_west_kalimantan.js script) to estimate CO2 emissions
 *    using IPCC Tier 1 combustion factors.
 *
 * ============================================================================
 * GEE-SPECIFIC CAVEATS
 * ============================================================================
 *
 * - NASA/LANCE/SNPP_VIIRS/C2 is the Suomi-NPP VIIRS NRT active fire product
 *   at 375m resolution, updated with a ~6-24h lag. It is distinct from the
 *   MODIS 'FIRMS' ImageCollection (1km). The Mask band encodes fire confidence:
 *   7 = nominal, 8-9 = high. Non-fire pixels (Mask 0, 3, 5) are masked out.
 *
 * - The Otsu histogram is computed at 90 m scale to avoid memory limits.
 *   If the histogram is empty (all NaN), lower L_CLOUD_MAX or widen the
 *   date range so the composite has coverage.
 *
 * - reduceRegion() at 30 m over all of Kalimantan may take 1-3 minutes to
 *   evaluate. If the statistics panel remains on "computing..." for more
 *   than 5 minutes, the computation has likely timed out. Increase tileScale
 *   to 8 or reduce scale to 60 m in computeBurnAreas().
 *
 * - GEE Apps deployed via the Apps menu run in a user-isolated environment.
 *   The 'aoi' import must be embedded as a hardcoded geometry (e.g.,
 *   ee.Geometry.BBox(108.0, -4.1, 119.0, 4.2)) rather than an import
 *   variable for a fully standalone deployable App.
 * ============================================================================ */
