/* ============================================================================
 *  Kalimantan Wildfire Monitoring App -- August 2026
 *
 *  Monitors widespread wildfire events across Kalimantan (Indonesian Borneo)
 *  using multi-source satellite data:
 *    - VIIRS FIRMS 375m      : active fire hotspot detection (last 14 days)
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
var START_DATE     = '2026-08-07';   // Post-fire window start (14-day lookback)
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

// Kalimantan province names matching FAO/GAUL/2015/level1 ADM1_NAME field
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
 * Loads Kalimantan provincial boundaries from FAO GAUL Level 1 dataset.
 * Filters Indonesia features to the five Kalimantan administrative units.
 *
 * @returns {ee.FeatureCollection} Five Kalimantan province polygons
 */
function loadProvinces() {
  return ee.FeatureCollection('FAO/GAUL/2015/level1')
    .filter(ee.Filter.eq('ADM0_NAME', 'Indonesia'))
    .filter(ee.Filter.inList('ADM1_NAME', KALIMANTAN_ADM1));
}

// ============================================================================
// 4. VIIRS FIRMS -- Active fire hotspot detections
// ============================================================================
/**
 * Loads VIIRS 375m active fire detections from the GEE FIRMS FeatureCollection.
 * Filters by date range, spatial bounds, and confidence level.
 *
 * VIIRS confidence encoding in the FIRMS FeatureCollection:
 *   'l' = low  |  'n' = nominal  |  'h' = high
 * Only nominal and high detections are loaded to suppress false positives,
 * consistent with Urbanski (2018) and Giglio et al. (2025).
 *
 * @param {ee.Geometry} geometry  - Spatial filter geometry
 * @param {string}      startDate - Inclusive start date (YYYY-MM-DD)
 * @param {string}      endDate   - Inclusive end date (YYYY-MM-DD)
 * @returns {ee.FeatureCollection} Filtered VIIRS fire detection points
 */
function loadVIIRS(geometry, startDate, endDate) {
  var firms = ee.FeatureCollection('FIRMS')
    .filterDate(startDate, endDate)
    .filterBounds(geometry)
    .filter(ee.Filter.inList('confidence', ['n', 'h']));
  print('VIIRS detections (nominal + high, ' + startDate + ' to ' + endDate + '):', firms.size());
  return firms;
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
  var nbr   = nir.subtract(swir2).divide(nir.add(swir2));
  return nbr.multiply(swir2.divide(nir)).rename('ATBI');
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

  var idx     = between.argmax().get([0]);
  var rawT    = ee.Number(vals.get([idx]));

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

  return ee.Image(0)
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
 * @param {ee.Image}     datbi         - dATBI image (for inspector)
 * @param {ee.Image}     burnSeverity  - Classified severity image (for inspector)
 * @param {ee.Image}     sarChange     - SAR dVV/dVH change image (for inspector)
 * @param {ui.Map.Layer} severityLayer - Reference to severity map layer object
 * @param {ui.Map.Layer} sarLayer      - Reference to SAR change map layer object
 * @param {ee.Number}    otsuVal       - Computed Otsu threshold (for display)
 * @returns {ui.Panel} Constructed left panel
 */
function buildLeftPanel(datbi, burnSeverity, sarChange, severityLayer, sarLayer, otsuVal) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style : {width: '300px', padding: '12px', backgroundColor: 'white'}
  });

  // --- Title block ---
  panel.add(ui.Label('Kalimantan Wildfire Monitor',
    {fontWeight: 'bold', fontSize: '16px', margin: '0 0 2px 0'}));
  panel.add(ui.Label('August 7-21, 2026 -- Active Fires & Burn Scars',
    {fontSize: '12px', color: '#444', margin: '0 0 2px 0'}));
  panel.add(ui.Label('VIIRS 375m  |  Landsat 8/9 dATBI  |  Sentinel-1 SAR',
    {fontSize: '11px', color: '#666', margin: '0 0 6px 0'}));
  panel.add(ui.Label('Muhammad Wahyu Ramadhan',
    {fontSize: '11px', margin: '0'}));
  panel.add(ui.Label('github.com/mwahyur46',
    {fontSize: '11px', color: '#1a73e8', margin: '0 0 4px 0'}));
  panel.add(divider());

  // --- Layer toggles (driven by Map.layers()) ---
  panel.add(ui.Label('Layers',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));

  for (var i = 0; i < Map.layers().length(); i++) {
    (function(idx) {
      var layer = Map.layers().get(idx);
      var cb = ui.Checkbox({
        label: layer.getName(),
        value: layer.getShown(),
        style: {fontSize: '11px', margin: '2px 0'}
      });
      cb.onChange(function(checked) { layer.setShown(checked); });
      panel.add(cb);
    })(i);
  }

  panel.add(divider());

  // --- Burn severity legend ---
  panel.add(ui.Label('Burn Severity (dATBI, Otsu adaptive)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));
  panel.add(ui.Thumbnail({
    image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
               .visualize({min: 0, max: 1, palette: PALETTE_SEVERITY}),
    params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
    style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
  }));
  panel.add(ui.Panel([
    ui.Label('Low (>T)',          {fontSize: '10px', margin: '2px 0'}),
    ui.Label('Moderate (>T+0.30)',{fontSize: '10px', margin: '2px 0',
      stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('High (>T+0.60)',   {fontSize: '10px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

  // Async display of computed Otsu threshold
  var otsuLabel = ui.Label('Burn threshold T: computing...',
    {fontSize: '10px', color: '#888', margin: '2px 0'});
  panel.add(otsuLabel);
  otsuVal.evaluate(function(v) {
    otsuLabel.setValue('Burn threshold T: ' + (v != null ? v.toFixed(4) : 'N/A') +
      '  (Otsu + ' + OTSU_BIAS + ' bias)');
  });

  panel.add(divider());

  // --- VIIRS confidence legend (colored swatches) ---
  panel.add(ui.Label('VIIRS Hotspot Confidence',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));

  function swatch(color, text) {
    return ui.Panel([
      ui.Label('', {
        backgroundColor: color, width: '14px', height: '14px',
        margin: '2px 6px 2px 0', padding: '0', border: '1px solid #ccc'
      }),
      ui.Label(text, {fontSize: '11px', margin: '2px 0'})
    ], ui.Panel.Layout.flow('horizontal'), {margin: '1px 0'});
  }
  panel.add(swatch(PALETTE_HOTSPOT.nominal, 'Nominal confidence'));
  panel.add(swatch(PALETTE_HOTSPOT.high,    'High confidence'));

  panel.add(divider());

  // --- SAR backscatter change legend ---
  panel.add(ui.Label('SAR Backscatter Change (dVV)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));
  panel.add(ui.Thumbnail({
    image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
               .visualize({min: 0, max: 1, palette: PALETTE_SAR}),
    params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
    style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
  }));
  panel.add(ui.Panel([
    ui.Label('< -3 dB',  {fontSize: '10px', margin: '2px 0'}),
    ui.Label('0 dB',     {fontSize: '10px', margin: '2px 0',
      stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('> +3 dB',  {fontSize: '10px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));
  panel.add(ui.Label('Negative dVV = possible surface change from fire.',
    {fontSize: '10px', color: '#888', margin: '2px 0'}));

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
 *   - Total VIIRS hotspot count
 *   - Hotspot breakdown by Kalimantan province
 *   - Estimated burned area (ha) per severity class
 *   - Total estimated burned area
 *   - SAR supplementary note
 *
 * All statistics are computed server-side and populated asynchronously via
 * .evaluate() callbacks, identical to the pattern used in the AGB/canopy
 * height script for regression metrics.
 *
 * @param {ee.FeatureCollection} hotspots   - VIIRS fire detections
 * @param {ee.FeatureCollection} provinces  - Kalimantan province polygons
 * @param {Object}               burnAreas  - {lowHa, modHa, highHa} ee.Numbers
 * @param {ee.Number}            otsuVal    - Computed Otsu threshold
 * @returns {ui.Panel} Constructed right panel
 */
function buildRightPanel(hotspots, provinces, burnAreas, otsuVal) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style : {width: '300px', padding: '12px', backgroundColor: 'white'}
  });

  // --- Header ---
  panel.add(ui.Label('Fire Statistics',
    {fontWeight: 'bold', fontSize: '17px', margin: '0 0 2px 0'}));
  panel.add(ui.Label('Kalimantan -- 14-day active fire analysis',
    {fontSize: '11px', color: '#666', margin: '0 0 4px 0'}));
  panel.add(divider());

  // --- Date range ---
  panel.add(ui.Label('Analysis Period',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));
  panel.add(ui.Label('Post-fire : ' + START_DATE + ' to ' + END_DATE,
    {fontSize: '11px', margin: '1px 0'}));
  panel.add(ui.Label('Pre-fire  : ' + PRE_FIRE_START + ' to ' + PRE_FIRE_END,
    {fontSize: '11px', margin: '1px 0'}));
  panel.add(ui.Label('VIIRS confidence: nominal + high only',
    {fontSize: '11px', color: '#666', margin: '1px 0'}));
  panel.add(divider());

  // --- VIIRS total hotspot count ---
  panel.add(ui.Label('VIIRS Active Fire Detections',
    {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 4px 0', color: '#222'}));

  var totalLabel = ui.Label('Total hotspots: computing...',
    {fontSize: '13px', margin: '2px 0', fontWeight: 'bold', color: '#e74c3c'});
  panel.add(totalLabel);
  hotspots.size().evaluate(function(n) {
    totalLabel.setValue('Total hotspots: ' + (n != null ? n.toLocaleString() : '0'));
  });

  // --- Hotspot count by province ---
  panel.add(ui.Label('Breakdown by Province:',
    {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));

  KALIMANTAN_ADM1.forEach(function(name) {
    var shortName = name.replace('Kalimantan ', 'Kal. ');
    var lbl = ui.Label(shortName + ': computing...',
      {fontSize: '11px', margin: '1px 0'});
    panel.add(lbl);

    var provGeom = provinces
      .filter(ee.Filter.eq('ADM1_NAME', name))
      .first()
      .geometry();
    hotspots.filterBounds(provGeom).size().evaluate(function(n) {
      lbl.setValue(shortName + ': ' + (n != null ? n.toLocaleString() : '0') + ' hotspots');
    });
  });

  panel.add(divider());

  // --- Burn scar area by severity ---
  panel.add(ui.Label('Estimated Burn Scar Area',
    {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 4px 0', color: '#222'}));
  panel.add(ui.Label('Landsat 8/9 dATBI  |  Waleed & Bilal (2026) Otsu threshold',
    {fontSize: '10px', color: '#666', margin: '0 0 2px 0'}));

  // Show the computed Otsu threshold value asynchronously
  var otsuStatLabel = ui.Label('Burn onset threshold T: computing...',
    {fontSize: '10px', color: '#888', margin: '0 0 6px 0', fontFamily: 'monospace'});
  panel.add(otsuStatLabel);
  otsuVal.evaluate(function(v) {
    otsuStatLabel.setValue('Burn onset threshold T = ' + (v != null ? v.toFixed(4) : 'N/A'));
  });

  var lowLabel  = ui.Label('Low  severity  (T to T+0.30): computing...',
    {fontSize: '12px', margin: '2px 0', color: '#e67e22'});
  var modLabel  = ui.Label('Moderate (T+0.30 to T+0.60): computing...',
    {fontSize: '12px', margin: '2px 0', color: '#e74c3c'});
  var highLabel = ui.Label('High     (> T+0.60): computing...',
    {fontSize: '12px', margin: '2px 0', color: '#922b21'});
  var totLabel  = ui.Label('Total burned area: computing...',
    {fontSize: '13px', margin: '6px 0 2px 0', fontWeight: 'bold', color: '#c0392b'});

  panel.add(lowLabel);
  panel.add(modLabel);
  panel.add(highLabel);
  panel.add(totLabel);

  burnAreas.lowHa.evaluate(function(v) {
    lowLabel.setValue('Low  severity  (T to T+0.30): ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.modHa.evaluate(function(v) {
    modLabel.setValue('Moderate (T+0.30 to T+0.60): ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.highHa.evaluate(function(v) {
    highLabel.setValue('High     (> T+0.60): ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });

  burnAreas.lowHa.add(burnAreas.modHa).add(burnAreas.highHa).evaluate(function(v) {
    totLabel.setValue('Total burned area: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });

  panel.add(ui.Label(
    'Scale: 30 m pixels  |  All dATBI >= T pixels included.',
    {fontSize: '10px', color: '#888', margin: '4px 0'}
  ));

  panel.add(divider());

  // --- SAR note ---
  panel.add(ui.Label('Sentinel-1 SAR (supplementary)',
    {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));
  panel.add(ui.Label(
    'VV backscatter change (post - pre, dB) is shown as a cloud-penetrating ' +
    'indicator layer. Negative dVV (red) may indicate fire-related surface ' +
    'alteration. Cross-reference with dATBI for confirmation. ' +
    'Ref: Siegert (2000), IW descending orbit, mean composite.',
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
  ));

  panel.add(divider());

  // --- Methodology note ---
  panel.add(ui.Label('Methodology',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'ATBI = ((NIR-SWIR2)/(NIR+SWIR2)) * (SWIR2/NIR)',
    {fontSize: '10px', color: '#444', margin: '1px 0', fontFamily: 'monospace'}
  ));
  panel.add(ui.Label(
    'dATBI = ATBIpre - ATBIpost  (positive = burned)',
    {fontSize: '10px', color: '#444', margin: '1px 0', fontFamily: 'monospace'}
  ));
  panel.add(ui.Label(
    'Landsat 8/9 C2 L2: SR_B5 (NIR 865nm), SR_B7 (SWIR2 2200nm). ' +
    'Cloud mask: QA_PIXEL bits 3+4. Combined ~8-day revisit. ' +
    'Otsu bias: +' + OTSU_BIAS + ', floor: ' + OTSU_MIN + '. ' +
    'Refs: Waleed & Bilal (2026); Kurbanov et al. (2022).',
    {fontSize: '10px', color: '#888', margin: '2px 0'}
  ));

  return panel;
}

// ============================================================================
// 15. MAIN -- load data, compute products, add layers, build UI
// ============================================================================

// --- Province boundaries ---
var provinces = loadProvinces();

// --- VIIRS hotspots ---
var hotspots = loadVIIRS(aoi, START_DATE, END_DATE);

// Separate by confidence for styled visualization
var hotspotsNominal = hotspots.filter(ee.Filter.eq('confidence', 'n'));
var hotspotsHigh    = hotspots.filter(ee.Filter.eq('confidence', 'h'));

// --- Landsat 8+9 composites ---
var lsPre  = getLandsatComposite(aoi, PRE_FIRE_START, PRE_FIRE_END);
var lsPost = getLandsatComposite(aoi, START_DATE,     END_DATE);

// --- ATBI and dATBI ---
var atbiPre  = computeATBI(lsPre);
var atbiPost = computeATBI(lsPost);
var datbi    = computeDATBI(atbiPre, atbiPost);

print('dATBI min/max:', datbi.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoi, scale: 500, maxPixels: 1e9
}));

// --- Adaptive Otsu threshold ---
var otsuVal = otsuThreshold(datbi, aoi);
print('Adaptive Otsu threshold T:', otsuVal);

// --- Burn severity ---
var burnSeverity = classifyBurnSeverity(datbi, otsuVal);

// --- Burn area statistics (pixel-area reduction) ---
var burnAreas = computeBurnAreas(burnSeverity, aoi);

// --- Sentinel-1 SAR ---
var sar = loadSAR(aoi, PRE_FIRE_START, PRE_FIRE_END, START_DATE, END_DATE);

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

// dATBI continuous (diagnostic layer)
Map.addLayer(
  datbi,
  {min: -0.5, max: 0.8, palette: PALETTE_DATBI},
  'dATBI (continuous)', false
);

// Burn severity -- classified (primary output)
// Added via ui.Map.Layer to retain a reference for the opacity slider
var severityLayer = ui.Map.Layer(
  burnSeverity,
  {min: 1, max: 3, palette: PALETTE_SEVERITY},
  'Burn Severity (classified)', true, 0.85
);
Map.layers().add(severityLayer);

// VIIRS hotspots -- nominal confidence
Map.addLayer(
  hotspotsNominal.style({color: PALETTE_HOTSPOT.nominal, pointSize: 3}),
  {}, 'VIIRS Hotspots -- Nominal', true
);

// VIIRS hotspots -- high confidence (larger dot, more saturated)
Map.addLayer(
  hotspotsHigh.style({color: PALETTE_HOTSPOT.high, pointSize: 4}),
  {}, 'VIIRS Hotspots -- High Confidence', true
);

// SAR backscatter change dVV (cloud-penetrating supplementary layer)
var sarLayer = ui.Map.Layer(
  sar.change.select('dVV'),
  {min: -6, max: 6, palette: PALETTE_SAR},
  'SAR Backscatter Change dVV (optional)', false, 0.7
);
Map.layers().add(sarLayer);

// ============================================================================
// BUILD AND MOUNT UI PANELS
// ============================================================================

// Left panel must be built AFTER all Map.addLayer() calls
var leftPanel  = buildLeftPanel(datbi, burnSeverity, sar.change, severityLayer, sarLayer, otsuVal);
var rightPanel = buildRightPanel(hotspots, provinces, burnAreas, otsuVal);

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
 * 3. VIIRS confidence encoding
 *    The FIRMS FeatureCollection confidence field encoding ('l', 'n', 'h')
 *    is specific to VIIRS 375m data. If your FIRMS dataset version stores
 *    confidence as integers (0-100), update the loadVIIRS() filter to use
 *    ee.Filter.gte('confidence', 30) for nominal-equivalent filtering.
 *    Check the first few features with: print(hotspots.first()).
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
 * 3. Peatland mask overlay
 *    Overlay the CIFOR/WRI global peatland layer (or Indonesian MoEF peat
 *    data) to flag hotspots and burn scars occurring on peat, which have
 *    disproportionate carbon emission implications (Afira 2022).
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
 * - The FIRMS FeatureCollection is updated with a ~24h lag relative to
 *   real-time satellite overpasses. For truly near-real-time monitoring,
 *   use the NASA FIRMS Web Map Service or the LANCE NRT FIRMS API directly.
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
