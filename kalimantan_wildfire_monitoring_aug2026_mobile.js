/* ============================================================================
 *  Kalimantan Wildfire Monitor -- Mobile App (August 2026)
 *
 *  Mobile-optimised layout: map fills full screen; layer controls, fire
 *  statistics, and references are accessible via a floating tab bar at the
 *  bottom of the map. No side panels. Designed for portrait smartphone use.
 *
 *  Data sources and methodology are identical to the desktop version
 *  (kalimantan_wildfire_monitoring_aug2026.js). Only the UI layer differs.
 *
 *  Burn index methodology:
 *    ATBI  = ((NIR-SWIR2)/(NIR+SWIR2)) * (SWIR2/NIR)   -- Waleed & Bilal (2026)
 *    dATBI = ATBIpre - ATBIpost  (positive = burned)
 *    Burn threshold from adaptive Otsu method (server-side)
 *
 *  Author  : Muhammad Wahyu Ramadhan
 *  GitHub  : github.com/mwahyur46
 *  LinkedIn: linkedin.com/in/mwahyur
 * ============================================================================ */

// ============================================================================
// 1. CONFIGURATION
// ============================================================================
var aoi = ee.Geometry.BBox(108.0, -4.1, 119.0, 4.2);

var _cutoffMs  = new Date('2026-09-30').getTime();
var _endMs     = Math.min(Date.now(), _cutoffMs);
function _fmt(ms) { return new Date(ms).toISOString().slice(0, 10); }
var END_DATE        = _fmt(_endMs);
var START_DATE      = _fmt(_endMs - 20 * 86400000);
var PRE_FIRE_START  = _fmt(_endMs - 51 * 86400000);
var PRE_FIRE_END    = _fmt(_endMs - 21 * 86400000);

var L_CLOUD_MAX     = 80;
var OTSU_BIAS       = 0.02;
var OTSU_MIN        = 0.03;
var SEV_MOD_OFFSET  = 0.30;
var SEV_HIGH_OFFSET = 0.60;
var BURN_MIN_HA     = 10;
var VECTORIZE_SCALE = 60;

// Province names matching FAO GAUL 2024 Level 1 gaul1_name field.
var KALIMANTAN_ADM1 = [
  'Kalimantan Barat',
  'Kalimantan Tengah',
  'Kalimantan Selatan',
  'Kalimantan Timur',
  'Kalimantan Utara'
];

var PALETTE_SEVERITY = ['#ffffb2', '#fd8d3c', '#e31a1c'];
var PALETTE_HOTSPOT  = {nominal: '#f7dc6f', high: '#e74c3c'};
var PALETTE_SAR      = ['#d73027', '#f7f7f7', '#1a9641'];
var PALETTE_DATBI    = ['#4575b4', '#f7f7f7', '#d73027'];
var PALETTE_PEAT     = ['#c7e9c0', '#238b45'];              // mosaic / peat-dominated

// ============================================================================
// 2. MAP INITIALIZATION
// ============================================================================
Map.setCenter(114, 0, 6);
Map.setOptions('HYBRID');
Map.style().set('cursor', 'crosshair');

// Hide controls that take up mobile screen real estate.
// Keep zoomControl and layerList; layerList provides native opacity control.
Map.setControlVisibility({
  mapTypeControl    : false,
  fullscreenControl : false,
  scaleControl      : false,
  zoomControl       : true,
  layerList         : true
});

// ============================================================================
// 3. PROVINCE BOUNDARIES
// ============================================================================
function loadProvinces() {
  return ee.FeatureCollection('projects/sat-io/open-datasets/FAO/GAUL/GAUL_2024_L1')
    .filter(ee.Filter.eq('iso3_code', 'IDN'))
    .filter(ee.Filter.inList('gaul1_name', KALIMANTAN_ADM1));
}

// ============================================================================
// 4. VIIRS NRT -- Active fire detections (375m, Suomi-NPP VIIRS C2)
// ============================================================================
function loadVIIRS(geometry, startDate, endDate) {
  var col = ee.ImageCollection('NASA/LANCE/SNPP_VIIRS/C2')
    .filterDate(startDate, endDate)
    .filterBounds(geometry)
    .select('confidence');

  var maxConf = col.max().clip(geometry);

  var nominal = maxConf.updateMask(maxConf.eq(1));
  var high    = maxConf.updateMask(maxConf.eq(2));
  var allFire = maxConf.updateMask(maxConf.gte(1));

  function pixelCount(img) {
    return img.reduceRegion({
      reducer  : ee.Reducer.count(),
      geometry : geometry,
      scale    : 375,
      maxPixels: 1e10,
      tileScale: 4
    }).getNumber('confidence');
  }

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
    totalCount   : pixelCount(allFire),
    provinceCount: provinceCount
  };
}

// ============================================================================
// 5. LANDSAT 8/9 -- Cloud masking + median composite
// ============================================================================
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

function getLandsatComposite(geometry, startDate, endDate) {
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(geometry).filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', L_CLOUD_MAX))
    .map(maskLandsatClouds);
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterBounds(geometry).filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', L_CLOUD_MAX))
    .map(maskLandsatClouds);
  return l8.merge(l9).median().clip(geometry);
}

// ============================================================================
// 6. ATBI AND dATBI COMPUTATION
// ============================================================================
function computeATBI(lsImage) {
  var nir   = lsImage.select('SR_B5');
  var swir2 = lsImage.select('SR_B7');
  var validMask = nir.gt(0.05).and(swir2.gt(0)).and(nir.lt(1.0)).and(swir2.lt(1.0));
  var nbr = nir.subtract(swir2).divide(nir.add(swir2));
  return nbr.multiply(swir2.divide(nir)).updateMask(validMask).rename('ATBI');
}

function computeDATBI(atbiPre, atbiPost) {
  return atbiPre.subtract(atbiPost).rename('dATBI');
}

// ============================================================================
// 7. ADAPTIVE OTSU THRESHOLDING
// ============================================================================
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

  var maxB    = between.reduce(ee.Reducer.max(), [0]);
  var maxBRep = maxB.repeat(0, ee.Number(n));
  var isMax   = between.eq(maxBRep);
  var rawT    = vals.multiply(isMax).reduce(ee.Reducer.sum(), [0]).get([0]);

  return rawT.add(OTSU_BIAS).max(OTSU_MIN);
}

// ============================================================================
// 8. BURN SEVERITY CLASSIFICATION
// ============================================================================
function classifyBurnSeverity(datbi, burnThreshold) {
  var T   = burnThreshold;
  var mod = T.add(SEV_MOD_OFFSET);
  var hi  = T.add(SEV_HIGH_OFFSET);
  return datbi.multiply(0)
    .where(datbi.gte(T).and(datbi.lt(ee.Image(mod))),  1)
    .where(datbi.gte(ee.Image(mod)).and(datbi.lt(ee.Image(hi))), 2)
    .where(datbi.gte(ee.Image(hi)), 3)
    .updateMask(datbi.gte(ee.Image(T)))
    .rename('burn_severity')
    .toInt();
}

// ============================================================================
// 9. SENTINEL-1 SAR
// ============================================================================
function loadSAR(geometry, preStart, preEnd, postStart, postEnd) {
  var s1Col = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(geometry)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
    .select(['VV', 'VH']);

  var s1Pre  = s1Col.filterDate(preStart, preEnd).mean().clip(geometry);
  var s1Post = s1Col.filterDate(postStart, postEnd).mean().clip(geometry);
  return {pre: s1Pre, post: s1Post, change: s1Post.subtract(s1Pre).rename(['dVV', 'dVH'])};
}

// ============================================================================
// 10. BURN AREA CALCULATION
// ============================================================================
function computeBurnAreas(burnSeverity, geometry) {
  var pixelArea = ee.Image.pixelArea().divide(10000);
  function areaForClass(cls) {
    return pixelArea.updateMask(burnSeverity.eq(cls)).rename('area_ha')
      .reduceRegion({
        reducer: ee.Reducer.sum(), geometry: geometry,
        scale: 30, maxPixels: 1e10, tileScale: 4
      }).getNumber('area_ha');
  }
  return {lowHa: areaForClass(1), modHa: areaForClass(2), highHa: areaForClass(3)};
}

// ============================================================================
// 11. VECTORIZE BURN SCARS (optional -- may timeout on large AOIs)
// ============================================================================
function vectorizeBurnScars(burnSeverity, geometry) {
  var vectors = burnSeverity.reduceToVectors({
    geometry: geometry, scale: VECTORIZE_SCALE, geometryType: 'polygon',
    eightConnected: true, labelProperty: 'burn_severity',
    maxPixels: 1e10, tileScale: 4
  });
  return vectors.filter(ee.Filter.gte('area', BURN_MIN_HA * 10000));
}

// ============================================================================
// 12. HELPER WIDGETS
// ============================================================================
function divider() {
  return ui.Label('', {
    height: '1px', backgroundColor: '#ddd', margin: '6px 0', stretch: 'horizontal'
  });
}

// Compact colored swatch + label row
function swatch(color, text) {
  return ui.Panel([
    ui.Label('', {
      backgroundColor: color, width: '12px', height: '12px',
      margin: '2px 5px 2px 0', padding: '0', border: '1px solid #ccc'
    }),
    ui.Label(text, {fontSize: '10px', margin: '1px 0'})
  ], ui.Panel.Layout.flow('horizontal'), {margin: '1px 0'});
}

// Gradient thumbnail legend bar
function legendBar(palette) {
  return ui.Thumbnail({
    image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
               .visualize({min: 0, max: 1, palette: palette}),
    params: {bbox: [-180, -1, 180, 1], dimensions: '220x12'},
    style : {stretch: 'horizontal', height: '12px', margin: '2px 0 0 0', padding: '0'}
  });
}

// ============================================================================
// 13. LAYERS TAB PANEL
// ============================================================================
/**
 * Builds the Layers tab panel with per-layer descriptions suitable for a
 * general (non-geospatial) audience, compact legends, Otsu threshold, and
 * a peatland mask toggle.
 * No opacity sliders -- use the native Layers button (top-right of map).
 * Must be called AFTER all Map.addLayer() calls.
 *
 * @param {ee.Number} otsuVal  - Adaptive Otsu threshold for display
 * @param {ee.Image}  peatMask - Binary peatland mask (1 = peat present)
 * @param {Object}    viirsObj - {nominal, high, ...} from loadVIIRS()
 * @returns {ui.Panel} Floating panel (initially hidden)
 */
function buildLayersTab(otsuVal, peatMask, viirsObj) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position       : 'bottom-left',
      width          : '340px',
      maxHeight      : '370px',
      backgroundColor: 'rgba(255,255,255,0.96)',
      padding        : '10px 12px',
      shown          : false,
      // Offset upward to clear the ~44px tab bar
      margin         : '0 0 50px 0'
    }
  });

  panel.add(ui.Label('Map Layers',
    {fontWeight: 'bold', fontSize: '14px', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'Toggle layers on or off. Tap the Layers button (top-right of map) to adjust opacity.',
    {fontSize: '9px', color: '#888', margin: '0 0 8px 0'}
  ));

  // Helper: one checkbox row + description line
  function layerRow(layerIdx, description) {
    var layer = Map.layers().get(layerIdx);
    var cb = ui.Checkbox({
      label: layer.getName(),
      value: layer.getShown(),
      style: {fontSize: '11px', fontWeight: 'bold', margin: '2px 0 0 0'}
    });
    cb.onChange(function(checked) { layer.setShown(checked); });
    panel.add(cb);
    panel.add(ui.Label(description,
      {fontSize: '9px', color: '#666', margin: '0 0 4px 8px'}));
  }

  // --- Province Boundaries (index 0) ---
  layerRow(0,
    'Outlines the five Kalimantan provinces. Always visible as a geographic reference.'
  );

  panel.add(divider());
  panel.add(ui.Label('Landsat Optical Imagery',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));

  // --- Landsat Pre-fire (index 1) ---
  layerRow(1,
    'True-color satellite image from July 2026, before the fires. ' +
    'Appears green because Kalimantan is mostly intact tropical forest. ' +
    'Use as a baseline to compare against the post-fire image.'
  );

  // --- Landsat Post-fire True Color (index 2) ---
  layerRow(2,
    'True-color satellite image from August 2026 showing conditions after the fires. ' +
    'Green areas are still-forested; darker brownish patches indicate burn scars. ' +
    'Clouds appear as white patches where the satellite could not see the ground.'
  );

  // --- Landsat Post-fire False Color SWIR (index 3) ---
  layerRow(3,
    'False-color image using shortwave infrared (SWIR) bands that are more ' +
    'sensitive to fire and bare soil. Active burn scars appear bright red-orange, ' +
    'making them easier to spot than in the true-color view.'
  );

  panel.add(divider());
  panel.add(ui.Label('Peatland Data',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));

  // --- Peatland Extent (index 4) ---
  layerRow(4,
    'Global Peatland Map 2.0 (1 km). Light green = peat in soil mosaic; ' +
    'dark green = peat-dominated. Source: Greifswald Mire Centre / Global Peatlands Initiative.'
  );
  panel.add(ui.Panel([
    swatch(PALETTE_PEAT[0], 'Peat in soil mosaic'),
    swatch(PALETTE_PEAT[1], 'Peat-dominated')
  ], ui.Panel.Layout.flow('vertical'), {margin: '0 0 4px 0'}));

  panel.add(divider());
  panel.add(ui.Label('Burn Scar Analysis',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));

  // --- dATBI burn signal only (index 5) ---
  layerRow(5,
    'Diagnostic: shows only where the satellite burn signal is positive, before ' +
    'the severity threshold is applied. Yellow = weak signal; red = strong burn signal. ' +
    'Areas with no signal (water, intact forest, cloud shadow) are left transparent. ' +
    'Compare with the Burn Severity layer to see how the threshold was applied.'
  );
  panel.add(legendBar(PALETTE_SEVERITY));
  panel.add(ui.Panel([
    ui.Label('Weak signal', {fontSize: '8px'}),
    ui.Label('',            {fontSize: '8px', stretch: 'horizontal'}),
    ui.Label('Strong burn', {fontSize: '8px'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal', margin: '0 0 4px 0'}));

  // --- Burn Severity (index 6) ---
  layerRow(6,
    'Fire damage intensity classified into three levels using the satellite-derived ' +
    'burn index. Low severity = partial scorching; high severity = intense burn with ' +
    'major vegetation loss. Clouds and unburned areas are left transparent.'
  );
  panel.add(legendBar(PALETTE_SEVERITY));
  panel.add(ui.Panel([
    ui.Label('Low',      {fontSize: '9px'}),
    ui.Label('Moderate', {fontSize: '9px', stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('High',     {fontSize: '9px'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal', margin: '0 0 2px 0'}));

  var otsuLbl = ui.Label('Scene threshold T: computing...',
    {fontSize: '9px', color: '#888', margin: '0 0 4px 0'});
  panel.add(otsuLbl);
  otsuVal.evaluate(function(v) {
    otsuLbl.setValue('Scene-adaptive threshold T: ' + (v != null ? v.toFixed(4) : 'N/A'));
  });

  panel.add(divider());
  panel.add(ui.Label('Active Fire Detections',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));

  // --- VIIRS Nominal (index 7) ---
  layerRow(7,
    'NASA VIIRS satellite detected an active fire signal at this 375m pixel ' +
    'with nominal (likely) confidence. Captured during the August 2026 dry season.'
  );
  panel.add(swatch(PALETTE_HOTSPOT.nominal, 'Nominal confidence -- likely fire'));

  // --- VIIRS High Confidence (index 8) ---
  layerRow(8,
    'Active fire detection with a strong, high-confidence thermal signal. ' +
    'These represent the most certain fire locations.'
  );
  panel.add(swatch(PALETTE_HOTSPOT.high, 'High confidence -- strong fire signal'));

  panel.add(divider());
  panel.add(ui.Label('Radar (Cloud-Penetrating)',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));

  // --- SAR dVV (index 9) ---
  layerRow(9,
    'Unlike optical cameras, radar from the Sentinel-1 satellite passes through ' +
    'cloud cover. A drop in radar signal (red) between July and August suggests ' +
    'the forest canopy was lost -- a known fire indicator. Use this to look for ' +
    'fire damage in areas the optical layers could not see through clouds.'
  );
  panel.add(legendBar(PALETTE_SAR));
  panel.add(ui.Panel([
    ui.Label('Signal loss (fire?)', {fontSize: '8px'}),
    ui.Label('No change',           {fontSize: '8px', stretch: 'horizontal', textAlign: 'center'}),
    ui.Label('Signal gain',         {fontSize: '8px'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

  panel.add(divider());

  // --- Peatland mask toggle ---
  panel.add(ui.Label('Peatland Mask',
    {fontWeight: 'bold', fontSize: '11px', color: '#333', margin: '0 0 2px 0'}));
  panel.add(ui.Label(
    'Restricts burn severity, dATBI, and VIIRS hotspot layers to peatland ' +
    'areas only. SAR and Landsat layers are unaffected.',
    {fontSize: '9px', color: '#555', margin: '0 0 2px 0'}
  ));
  panel.add(ui.Label(
    'Note: peatland data is at 1 km -- edges may appear blocky.',
    {fontSize: '9px', color: '#aaa', margin: '0 0 4px 0'}
  ));

  var peatToggle = ui.Checkbox({
    label: 'Restrict fire layers to peatland areas',
    value: false,
    style: {fontSize: '11px', fontWeight: 'bold', margin: '0 0 4px 0'}
  });
  panel.add(peatToggle);

  // Layer indices: 5 = dATBI, 6 = burn severity, 7 = VIIRS nominal, 8 = VIIRS high
  peatToggle.onChange(function(checked) {
    if (checked) {
      Map.layers().get(5).setEeObject(datbi.updateMask(datbi.gt(0)).updateMask(peatMask));
      Map.layers().get(6).setEeObject(burnSeverity.updateMask(peatMask));
      Map.layers().get(7).setEeObject(viirsObj.nominal.updateMask(peatMask));
      Map.layers().get(8).setEeObject(viirsObj.high.updateMask(peatMask));
    } else {
      Map.layers().get(5).setEeObject(datbi.updateMask(datbi.gt(0)));
      Map.layers().get(6).setEeObject(burnSeverity);
      Map.layers().get(7).setEeObject(viirsObj.nominal);
      Map.layers().get(8).setEeObject(viirsObj.high);
    }
  });

  return panel;
}

// ============================================================================
// 14. STATS TAB PANEL
// ============================================================================
/**
 * Builds the Stats tab panel: VIIRS fire pixel counts (total + per province),
 * estimated burned area per severity class, and key limitations note.
 *
 * @param {Object}               viirs      - From loadVIIRS()
 * @param {ee.FeatureCollection} provinces  - Kalimantan province polygons
 * @param {Object}               burnAreas  - From computeBurnAreas()
 * @param {ee.Number}            otsuVal    - Adaptive Otsu threshold
 * @returns {ui.Panel} Floating panel (initially hidden)
 */
function buildStatsTab(viirs, provinces, burnAreas, otsuVal) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position       : 'bottom-left',
      width          : '340px',
      maxHeight      : '350px',
      backgroundColor: 'rgba(255,255,255,0.96)',
      padding        : '10px 12px',
      shown          : false,
      margin         : '0 0 50px 0'
    }
  });

  panel.add(ui.Label('Fire Statistics',
    {fontWeight: 'bold', fontSize: '14px', margin: '0 0 2px 0'}));
  panel.add(ui.Label('Kalimantan -- August 2026',
    {fontSize: '10px', color: '#666', margin: '0 0 4px 0'}));

  // Event narrative for general audience
  panel.add(ui.Label(
    'In August 2026, widespread fires were detected across Kalimantan (Indonesian Borneo) ' +
    'during the annual dry season. Fires affected peatland and forest areas across ' +
    'multiple provinces. The causes are under investigation by relevant authorities. ' +
    'This map combines three satellite systems to track where fires occurred, ' +
    'how severely vegetation was burned, and where radar data reveals fire ' +
    'signals hidden under cloud cover.',
    {fontSize: '10px', color: '#333', margin: '0 0 6px 0'}
  ));

  panel.add(ui.Label(
    'Post-fire: ' + START_DATE + ' to ' + END_DATE + ' | ' +
    'Pre-fire ref: ' + PRE_FIRE_START + ' to ' + PRE_FIRE_END,
    {fontSize: '9px', color: '#888', margin: '0 0 6px 0'}
  ));

  panel.add(divider());

  // VIIRS fire detections
  panel.add(ui.Label('Active Fire Detections (VIIRS 375m)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 3px 0'}));
  panel.add(ui.Label(
    'Each pixel = 375m x 375m area with confirmed fire thermal signal. ' +
    'Nominal and high confidence only.',
    {fontSize: '9px', color: '#555', margin: '0 0 4px 0'}
  ));

  var totalLabel = ui.Label('Total: computing...',
    {fontSize: '13px', fontWeight: 'bold', color: '#e74c3c', margin: '0 0 4px 0'});
  panel.add(totalLabel);
  viirs.totalCount.evaluate(function(n) {
    totalLabel.setValue('Total: ' + (n != null ? n.toLocaleString() : '0') + ' fire pixels');
  });

  // Per-province breakdown
  panel.add(ui.Label('By province:', {fontWeight: 'bold', fontSize: '10px', margin: '0 0 2px 0'}));
  KALIMANTAN_ADM1.forEach(function(name) {
    var short = name.replace('Kalimantan ', 'Kal. ');
    var lbl = ui.Label(short + ': computing...', {fontSize: '10px', margin: '1px 0'});
    panel.add(lbl);
    var filtered = provinces.filter(ee.Filter.eq('gaul1_name', name));
    var provGeom = ee.Geometry(ee.Algorithms.If(
      filtered.size().gt(0),
      filtered.first().geometry(),
      ee.Geometry.Point([0, 0])
    ));
    viirs.provinceCount(provGeom).evaluate(function(n) {
      lbl.setValue(short + ': ' + (n != null ? n.toLocaleString() : '0') + ' px');
    });
  });

  panel.add(divider());

  // Burn area by severity
  panel.add(ui.Label('Estimated Burned Area (Landsat)',
    {fontWeight: 'bold', fontSize: '12px', margin: '0 0 3px 0'}));
  panel.add(ui.Label(
    'Mapped using dATBI index (Waleed & Bilal 2026). Cloud-covered areas excluded.',
    {fontSize: '9px', color: '#555', margin: '0 0 2px 0'}
  ));

  var otsuStatLbl = ui.Label('Scene-adaptive T: computing...',
    {fontSize: '9px', color: '#888', margin: '0 0 4px 0'});
  panel.add(otsuStatLbl);
  otsuVal.evaluate(function(v) {
    otsuStatLbl.setValue('Scene-adaptive T = ' + (v != null ? v.toFixed(4) : 'N/A'));
  });

  var lowLbl  = ui.Label('Low: computing...',
    {fontSize: '11px', color: '#e67e22', margin: '1px 0'});
  var modLbl  = ui.Label('Moderate: computing...',
    {fontSize: '11px', color: '#e74c3c', margin: '1px 0'});
  var highLbl = ui.Label('High: computing...',
    {fontSize: '11px', color: '#922b21', margin: '1px 0'});
  var totLbl  = ui.Label('Total: computing...',
    {fontSize: '13px', fontWeight: 'bold', color: '#c0392b', margin: '4px 0 2px 0'});
  panel.add(lowLbl); panel.add(modLbl); panel.add(highLbl); panel.add(totLbl);

  burnAreas.lowHa.evaluate(function(v) {
    lowLbl.setValue('Low: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.modHa.evaluate(function(v) {
    modLbl.setValue('Moderate: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.highHa.evaluate(function(v) {
    highLbl.setValue('High: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });
  burnAreas.lowHa.add(burnAreas.modHa).add(burnAreas.highHa).evaluate(function(v) {
    totLbl.setValue('Total: ' + (v != null ? Math.round(v).toLocaleString() : '0') + ' ha');
  });

  panel.add(ui.Label(
    'Cloud gaps excluded -- actual burned extent may be larger.',
    {fontSize: '9px', color: '#c0392b', margin: '2px 0 0 0'}
  ));

  panel.add(divider());

  // Compact limitations
  panel.add(ui.Label('Key Limitations',
    {fontWeight: 'bold', fontSize: '11px', color: '#922b21', margin: '0 0 3px 0'}));
  [
    '1. Cloud gaps: burned areas under clouds are missed entirely.',
    '2. Timing lag: recent fires may not yet appear (1-2 day ingestion delay).',
    '3. False positives: bare agricultural land can mimic burn scars.',
    '4. Peatland sub-surface fires may be underestimated by optical sensors.'
  ].forEach(function(txt) {
    panel.add(ui.Label(txt, {fontSize: '9px', color: '#555', margin: '1px 0 1px 4px'}));
  });

  panel.add(divider());

  panel.add(ui.Label('Potential Improvements',
    {fontWeight: 'bold', fontSize: '11px', color: '#1a5276', margin: '0 0 3px 0'}));
  [
    '1. SAR-optical fusion: combining radar and optical burn signals would ' +
    'fill cloud gaps and reduce false positives (partially addressed by the SAR layer).',
    '2. Carbon emission estimate: combining mapped burn area with the ' +
    'companion mangrove biomass analysis would allow a first-order estimate ' +
    'of CO₂ released.',
    '3. Time-series monitoring: tracking VIIRS fire counts day by day through ' +
    'the dry season would show how the event evolved over time.'
  ].forEach(function(txt) {
    panel.add(ui.Label(txt, {fontSize: '9px', color: '#555', margin: '1px 0 3px 4px'}));
  });

  return panel;
}

// ============================================================================
// 15. REFERENCES TAB PANEL
// ============================================================================
/**
 * Builds the References tab panel: data sources and full citations.
 *
 * @returns {ui.Panel} Floating panel (initially hidden)
 */
function buildRefsTab() {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position       : 'bottom-left',
      width          : '340px',
      maxHeight      : '350px',
      backgroundColor: 'rgba(255,255,255,0.96)',
      padding        : '10px 12px',
      shown          : false,
      margin         : '0 0 50px 0'
    }
  });

  panel.add(ui.Label('Data & References',
    {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}));

  // Data sources
  panel.add(ui.Label('Data Sources', {fontWeight: 'bold', fontSize: '11px', margin: '0 0 2px 0'}));
  [
    'VIIRS NRT: NASA/LANCE/SNPP_VIIRS/C2 (375m)',
    'Landsat 8/9: USGS Collection 2 Level-2 SR (30m)',
    'Sentinel-1: ESA Copernicus GRD IW (10m)',
    'Boundaries: FAO GAUL 2024 Level 1 (Franceschini et al. 2025)',
    'Peatlands: Global Peatland Map 2.0 -- Global Peatlands Initiative / COP26 (1 km)'
  ].forEach(function(s) {
    panel.add(ui.Label(s, {fontSize: '9px', color: '#444', margin: '1px 0 1px 4px'}));
  });

  panel.add(divider());

  // Key references
  panel.add(ui.Label('References', {fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0'}));
  [
    {text: 'Afira, N., & Wijayanto, A. W. (2022). Mono-temporal and multi-temporal approaches for burnt area detection using Sentinel-2 satellite imagery (a case study of Rokan Hilir Regency, Indonesia). Ecological Informatics, 69, 101677.',
     doi:  'https://doi.org/10.1016/j.ecoinf.2022.101677'},
    {text: 'Franceschini, G., Khan, A., Moretti, L., Nyabuti, K., Asif, M., Bezuidenhoudt, E., & Morteo, K. (2025). The Global Administrative Unit Layers (GAUL) 2024. Technical guidelines. Rome, FAO.',
     doi:  'https://doi.org/10.4060/cd4262en'},
    {text: 'Giglio, L., Boschetti, L., Roy, D. P., Hall, J. V., Zubkova, M., Humber, M., Huang, H., & Oles, V. (2025). The NASA VIIRS burned area product, global validation, and intercomparison with the NASA MODIS burned area product. Remote Sensing of Environment, 331, 115006.',
     doi:  'https://doi.org/10.1016/j.rse.2025.115006'},
    {text: 'Greifswald Mire Centre (2022). Global Peatland Map 2.0. Underlying dataset of the UNEP Global Peatland Assessment -- The State of the World\'s Peatlands: Evidence for action toward the conservation, restoration, and sustainable management of peatlands, Global Peatlands Initiative, United Nations Environment Programme, Nairobi.',
     doi:  'https://www.greifswald-moor-centrum.de/en/services/gis-data/global-peatland-map-2-0/'},
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
     doi:  'https://doi.org/10.1016/j.jag.2026.105517'}
  ].forEach(function(ref) {
    panel.add(ui.Label(ref.text,
      {fontSize: '9px', color: '#333', margin: '4px 0 0 0', fontWeight: 'bold'}));
    panel.add(ui.Label(ref.doi,
      {fontSize: '9px', color: '#1a73e8', margin: '0 0 2px 4px', fontFamily: 'monospace'}));
  });

  panel.add(divider());
  panel.add(ui.Label('Author: Muhammad Wahyu Ramadhan',
    {fontSize: '9px', color: '#555', margin: '0'}));
  panel.add(ui.Label('github.com/mwahyur46',
    {fontSize: '9px', color: '#1a73e8', margin: '0'}));

  return panel;
}

// ============================================================================
// 16. PIXEL INSPECTOR PANEL (map-click overlay, top-right)
// ============================================================================
/**
 * Builds the pixel inspector overlay. Wires a Map.onClick handler to
 * sample dATBI, burn severity, and SAR dVV at the tapped location and
 * display results in a compact floating card (top-right).
 *
 * @param {ee.Image} datbi        - dATBI image
 * @param {ee.Image} burnSeverity - Classified severity image
 * @param {ee.Image} sarChange    - SAR dVV/dVH change image
 * @returns {ui.Panel} Inspector panel added to Map
 */
function buildInspectorPanel(datbi, burnSeverity, sarChange) {
  var panel = ui.Panel({
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position       : 'top-right',
      width          : '155px',
      backgroundColor: 'rgba(255,255,255,0.95)',
      padding        : '6px 8px',
      shown          : false
    }
  });

  var SEVERITY_LABELS = ['', 'Low', 'Moderate', 'High'];
  var SEVERITY_COLORS = ['', '#e67e22', '#e74c3c', '#922b21'];

  Map.onClick(function(coords) {
    panel.clear();
    panel.style().set('shown', true);

    panel.add(ui.Label(
      coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4),
      {fontSize: '9px', color: '#888', margin: '0 0 3px 0'}
    ));
    panel.add(ui.Label('Sampling...', {fontSize: '10px', color: '#aaa'}));

    var point     = ee.Geometry.Point([coords.lon, coords.lat]);
    var sampleImg = datbi.addBands(burnSeverity).addBands(sarChange.select('dVV'));

    sampleImg.reduceRegion({
      reducer: ee.Reducer.first(), geometry: point, scale: 30
    }).evaluate(function(vals) {
      panel.clear();
      panel.add(ui.Label(
        coords.lon.toFixed(4) + ', ' + coords.lat.toFixed(4),
        {fontSize: '9px', color: '#888', margin: '0 0 3px 0'}
      ));

      if (!vals || vals['dATBI'] === null || vals['dATBI'] === undefined) {
        panel.add(ui.Label('No optical data here.', {fontSize: '10px', color: '#c00'}));
        return;
      }

      var sevCode  = vals['burn_severity'];
      var sevLabel = (sevCode == null) ? 'Unburned' : SEVERITY_LABELS[sevCode];
      var sevColor = (sevCode == null) ? '#27ae60'  : SEVERITY_COLORS[sevCode];

      panel.add(ui.Label('dATBI: ' + vals['dATBI'].toFixed(4),
        {fontSize: '11px', fontWeight: 'bold', color: '#c0392b', margin: '1px 0'}));
      panel.add(ui.Label('Severity: ' + sevLabel,
        {fontSize: '11px', fontWeight: 'bold', color: sevColor, margin: '1px 0'}));
      if (vals['dVV'] != null) {
        panel.add(ui.Label('SAR dVV: ' + vals['dVV'].toFixed(2) + ' dB',
          {fontSize: '11px', fontWeight: 'bold', color: '#1a73e8', margin: '1px 0'}));
      }

      // Dismiss button
      var closeBtn = ui.Button({
        label: 'x',
        style: {fontSize: '9px', padding: '1px 4px', margin: '4px 0 0 0',
                backgroundColor: '#eee', color: '#555'}
      });
      closeBtn.onClick(function() { panel.style().set('shown', false); });
      panel.add(closeBtn);
    });
  });

  return panel;
}

// ============================================================================
// 17. MAIN -- load data, compute products, add layers, build mobile UI
// ============================================================================

// --- Compute analysis products ---
var provinces    = loadProvinces();

// Global Peatland Map 2.0 (1 km resolution)
// Pixel values: 1 = peat-dominated, 2 = peat in soil mosaic
var peatRaw  = ee.Image('projects/sat-io/open-datasets/ML-GLOBAL-PEATLAND-EXTENT')
                 .clip(aoi)
                 .unmask(0);
var peatMask = peatRaw.gte(1);

var viirs        = loadVIIRS(aoi, START_DATE, END_DATE);
var lsPre        = getLandsatComposite(aoi, PRE_FIRE_START, PRE_FIRE_END);
var lsPost       = getLandsatComposite(aoi, START_DATE, END_DATE);
var atbiPre      = computeATBI(lsPre);
var atbiPost     = computeATBI(lsPost);
var datbi        = computeDATBI(atbiPre, atbiPost);
var otsuVal      = otsuThreshold(datbi, aoi);
var burnSeverity = classifyBurnSeverity(datbi, otsuVal);
var burnAreas    = computeBurnAreas(burnSeverity, aoi);
var sar          = loadSAR(aoi, PRE_FIRE_START, PRE_FIRE_END, START_DATE, END_DATE);

// SAR speckle reduction: 60m focal median suppresses pixel-level coherent noise
// (speckle) while preserving real change signals that span fire-affected patches
// (typically >100m). Applied to the change image only, not the raw composites.
var sarChangeSm  = sar.change.focal_median(60, 'circle', 'meters');

// ============================================================================
// ADD MAP LAYERS
// ============================================================================

Map.addLayer(
  provinces.style({color: '#ffffff', fillColor: '00000000', width: 1.5}),
  {}, 'Province Boundaries', true
);

Map.addLayer(
  lsPre,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0, max: 0.25, gamma: 0.9},
  'Landsat Pre-fire True Color (Jul 2026)', false
);

Map.addLayer(
  lsPost,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0, max: 0.25, gamma: 0.9},
  'Landsat Post-fire True Color (Aug 2026)', true
);

Map.addLayer(
  lsPost,
  {bands: ['SR_B7','SR_B5','SR_B4'], min: 0, max: 0.4, gamma: 0.9},
  'Landsat Post-fire False Color SWIR', false
);

// Peatland extent (Global Peatland Map 2.0, 1 km) -- index 4
// Placed above Landsat base imagery but below analysis layers.
Map.addLayer(
  peatRaw.updateMask(peatMask)
         .visualize({min: 1, max: 2, palette: PALETTE_PEAT, opacity: 0.35}),
  {}, 'Peatland Extent (Global Peatland Map 2.0)', false
);

// dATBI diagnostic -- index 5
Map.addLayer(
  datbi.updateMask(datbi.gt(0)),
  {min: 0, max: 0.4, palette: PALETTE_SEVERITY},
  'dATBI (burn signal only)', false
);

// Burn Severity -- index 6
Map.addLayer(
  burnSeverity,
  {min: 1, max: 3, palette: PALETTE_SEVERITY},
  'Burn Severity (classified)', true, 0.85
);

// VIIRS Nominal -- index 7
Map.addLayer(
  viirs.nominal,
  {min: 1, max: 1, palette: [PALETTE_HOTSPOT.nominal]},
  'VIIRS Hotspots -- Nominal (conf=1)', true
);

// VIIRS High Confidence -- index 8
Map.addLayer(
  viirs.high,
  {min: 2, max: 2, palette: [PALETTE_HOTSPOT.high]},
  'VIIRS Hotspots -- High Confidence (conf=2)', true
);

// SAR dVV: speckle-filtered change image, ±3 dB range so fire-related canopy
// loss signals (-1 to -3 dB) render as clearly visible red. -- index 9
Map.addLayer(
  sarChangeSm.select('dVV'),
  {min: -3, max: 3, palette: PALETTE_SAR},
  'SAR Backscatter Change dVV', false, 0.7
);

// ============================================================================
// BUILD AND MOUNT MOBILE UI OVERLAYS
// ============================================================================

// --- Tab content panels (hidden by default, mounted to map) ---
// Must be built AFTER addLayer calls so layer checkboxes reflect Map.layers()
var layersTab = buildLayersTab(otsuVal, peatMask, viirs);
var statsTab  = buildStatsTab(viirs, provinces, burnAreas, otsuVal);
var refsTab   = buildRefsTab();

Map.add(layersTab);
Map.add(statsTab);
Map.add(refsTab);

// --- Pixel inspector (shown on map tap) ---
var inspectorPanel = buildInspectorPanel(datbi, burnSeverity, sarChangeSm);
Map.add(inspectorPanel);

// --- Compact title bar (top-center) ---
var titleBar = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    position       : 'top-center',
    padding        : '5px 12px',
    backgroundColor: 'rgba(0,0,0,0.68)',
    margin         : '0'
  }
});
titleBar.add(ui.Label('Kalimantan Wildfire Monitor',
  {fontWeight: 'bold', fontSize: '13px', color: 'white', margin: '0'}));
titleBar.add(ui.Label('Aug 2026 -- VIIRS | Landsat | SAR',
  {fontSize: '10px', color: '#ccc', margin: '0'}));
Map.add(titleBar);

// --- Tab bar (bottom-center) ---
// Each button toggles its target panel; tapping an open tab closes it.
var openPanel = null;

function makeTabButton(label, targetPanel, allPanels) {
  var btn = ui.Button({
    label: label,
    style: {
      fontSize       : '11px',
      padding        : '6px 12px',
      margin         : '0 3px',
      backgroundColor: '#f8f9fa',
      color          : '#333'
    }
  });
  btn.onClick(function() {
    if (openPanel === targetPanel) {
      targetPanel.style().set('shown', false);
      openPanel = null;
    } else {
      allPanels.forEach(function(p) { p.style().set('shown', false); });
      targetPanel.style().set('shown', true);
      openPanel = targetPanel;
    }
  });
  return btn;
}

var allTabPanels = [layersTab, statsTab, refsTab];

var tabBar = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {
    position       : 'bottom-center',
    padding        : '5px 8px',
    backgroundColor: 'rgba(255,255,255,0.95)',
    margin         : '0'
  }
});
tabBar.add(makeTabButton('Layers',     layersTab, allTabPanels));
tabBar.add(makeTabButton('Stats',      statsTab,  allTabPanels));
tabBar.add(makeTabButton('References', refsTab,   allTabPanels));
Map.add(tabBar);

// ============================================================================
// END OF SCRIPT
// ============================================================================
