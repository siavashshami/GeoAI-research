// ==========================
// 0. AOI & TEMPORAL SETTINGS
// ==========================
var zone = ee.FeatureCollection('projects/ee-siavashshami/assets/Cs');
var zoneGeom = zone.geometry();

var startYear = 2016;
var endYear   = 2016;   
var years = ee.List.sequence(startYear, endYear);

var climateZoneName = 'Cs';
var analysisScale = 2000;  

Map.centerObject(zone, 5);
Map.addLayer(zone, {color: 'red'}, 'AOI: ' + climateZoneName);

// ==========================
// 1. LAND COVER FUNCTIONS (MODIS MCD12Q1 → 12 merged classes)
// ==========================
var mergedClassNames = [
  'EvergreenForest', 'DeciduousForest', 'MixedForest', 'Shrublands',
  'Savanna', 'Grasslands', 'Wetlands', 'Agriculture', 'Urban',
  'SnowIce', 'Barren', 'Water'
];

function getLandCoverImage(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = startDate.advance(1, 'year');
  return ee.ImageCollection('MODIS/061/MCD12Q1')
    .filterDate(startDate, endDate)
    .filterBounds(zoneGeom)
    .first()
    .select('LC_Type1')
    .clip(zoneGeom);
}

function reclassifyToMerged(image) {
  return image
    .updateMask(image.neq(0).and(image.neq(254)).and(image.neq(255)))
    .remap(
      [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
      [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 8, 9, 8, 10, 11, 12]
    )
    .rename('merged_LC');
}

// ==========================
// 1A. LAND COVER VALIDATION (FIRST YEAR)
// ==========================
var lcOriginal = getLandCoverImage(startYear);
var lcMerged   = reclassifyToMerged(lcOriginal);

print('LC ' + startYear + ' original:', lcOriginal);
print('LC ' + startYear + ' merged:', lcMerged);

print('ORIGINAL LC HISTOGRAM (17 classes):',
  lcOriginal.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: zoneGeom,
    scale: analysisScale,
    maxPixels: 1e8
  })
);

print('MERGED LC HISTOGRAM (12 classes):',
  lcMerged.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: zoneGeom,
    scale: analysisScale,
    maxPixels: 1e8
  })
);

Map.addLayer(
  lcOriginal,
  {min: 0, max: 17, palette: [
    '000000','006400','008000','00FF00','8B4513','FFA500','FFFF00',
    'FFC0CB','808080','00FFFF','FF00FF','00FF7F','FF0000','000080',
    'FFFFFF','0000FF','800080'
  ]},
  'LC Original (17 classes)'
);

Map.addLayer(
  lcMerged,
  {min: 1, max: 12, palette: [
    '004000','008000','00C000','8B4513','FFA500','FFFF00',
    '87CEEB','FFC0CB','A9A9A9','FFFFFF','FF0000','0000FF'
  ]},
  'LC Merged (12 classes)'
);

// ==========================
// K-FACTOR: USDA NRCS VALIDATED (Renard et al. 1997)
// ==========================
function computeKFactorScientific(geometry, scale) {
  var texture_raw = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").select('b0');
  
  var texture = texture_raw
    .clip(geometry.buffer(10000))
    .resample('bilinear')
    .reproject({crs: 'EPSG:4326', scale: scale});
  
  var K = texture.remap(
    [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12],  
    [0.41,0.45,0.31,0.34,0.36,0.29,0.28,0.26,0.20,0.25,0.17,0.15], 
    0
  ).rename('K');
  
  var K_valid = K.updateMask(K.gt(0).and(K.lte(0.69)));
  
  return K_valid.addBands([
    texture.rename('texture_class'),
    K_valid.gt(0).rename('K_valid_mask')
  ]).set('source', 'OpenLandMap_v02 + USDA_NRCS_1997');
}

// ==========================
// R-FACTOR: CHIRPS DAILY (Renard et al. 1997)
// ==========================
function temporalMonthlySum(collection, startDate) {
  var months = ee.List.sequence(0, 11);
  var start  = ee.Date(startDate);
  return ee.ImageCollection(
    months.map(function(m) {
      m = ee.Number(m);
      var s = start.advance(m, 'month');
      var e = s.advance(1, 'month');
      return collection.filterDate(s, e).sum().set('system:time_start', s.millis());
    })
  );
}

function computeRFactor(year, geometry) {
  var startYearDate = ee.Date.fromYMD(year, 1, 1);
  var endYearDate   = startYearDate.advance(1, 'year');
  var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    .filterBounds(geometry)
    .filterDate(startYearDate, endYearDate)
    .select('precipitation');
  var monthlyRain = temporalMonthlySum(chirps, startYearDate);
  var annualRain = monthlyRain.reduce(ee.Reducer.sum()).rename('P').clip(geometry);
  var P_safe = annualRain.max(1);
  var monthlySquares = monthlyRain.map(function(img) { return img.pow(2); });
  var sumMonthlySquares = monthlySquares.reduce(ee.Reducer.sum());
  var MFI = sumMonthlySquares.divide(P_safe).rename('MFI');
  var R_raw = MFI.expression(
    '1.735 * pow(10, 1.5 * log10(MFI) - 0.8188)',
    {MFI: MFI}
  ).rename('R');
  var R = R_raw.where(annualRain.eq(0), 0).rename('R');
  return R.clip(geometry);
}

// ==========================
// C-FACTOR: HEMISPHERE-AWARE (Van der Knijff 2000)
// ==========================
function computeCfactorScientific(year, geometry) {
  var startYear = ee.Date.fromYMD(year, 1, 1);
  var endYear   = startYear.advance(1, 'year');
  var endYear2  = startYear.advance(2, 'year');
  
  var ndviAnnualRaw = ee.ImageCollection('MODIS/061/MOD13Q1')
    .filterBounds(geometry)
    .filterDate(startYear, endYear)
    .select('NDVI')
    .median()
    .clip(geometry);
  var ndviAnnual = ndviAnnualRaw.multiply(0.0001).clamp(0, 1);
  
  var centroid = geometry.bounds().centroid({maxError: 100});
  var coords   = ee.List(centroid.coordinates());
  var latitude = ee.Number(coords.get(1));
  var isSouthern = latitude.lt(0);
  
  var growingSeasonNDVI = ee.Image(ee.Algorithms.If(
    isSouthern,
    ee.ImageCollection('MODIS/061/MOD13Q1')
      .filterBounds(geometry)
      .filter(ee.Filter.and(
        ee.Filter.date(startYear, endYear),
        ee.Filter.calendarRange(11, 12, 'month')
      ))
      .merge(
        ee.ImageCollection('MODIS/061/MOD13Q1')
          .filterBounds(geometry)
          .filter(ee.Filter.and(
            ee.Filter.date(endYear, endYear2),
            ee.Filter.calendarRange(1, 5, 'month')
          ))
      )
      .select('NDVI')
      .median()
      .clip(geometry),
    ee.ImageCollection('MODIS/061/MOD13Q1')
      .filterBounds(geometry)
      .filterDate(startYear, endYear)
      .filter(ee.Filter.calendarRange(4, 10, 'month'))
      .select('NDVI')
      .median()
      .clip(geometry)
  ));
  
  var growingNdviRaw = growingSeasonNDVI;
  var growingNdvi    = growingNdviRaw.multiply(0.0001).clamp(0, 1);
  
  var C_annual = ndviAnnual.expression('exp(-2 * ndvi / (1 - ndvi))', {ndvi: ndviAnnual});
  var C_growing = growingNdvi.expression('exp(-2 * ndvi / (1 - ndvi))', {ndvi: growingNdvi});
  var C_scientific = C_annual.multiply(0.3).add(C_growing.multiply(0.7)).clamp(0.001, 0.95).rename('C_scientific');
  
  return C_scientific.addBands([
    growingNdviRaw.rename('growing_NDVI_raw'),
    ndviAnnualRaw.rename('annual_NDVI_raw')
  ]);
}

// ==========================
// LS-FACTOR: MOORE & BURCH (1986) - CORRECTED VERSION
// Uses HydroSHEDS flow accumulation in number of cells,
// cell size in meters, and As = N * cellSize.
// ==========================
function computeLSFactorScientific(geometry, analysisScale) {
  var dem = ee.Image('USGS/SRTMGL1_003')
    .clip(geometry)
    .resample('bilinear')
    .reproject({crs: 'EPSG:4326', scale: analysisScale});
  
  var slope = ee.Terrain.slope(dem);
  var slope_rad = slope.multiply(Math.PI/180);
  
  // HydroSHEDS flow accumulation (units: 100 * number of cells)
  var flowAcc_raw = ee.Image('WWF/HydroSHEDS/15ACC').clip(geometry);
  var cellSize_m = ee.Number(analysisScale); // e.g., 1000 m
  
  // Convert to number of upstream cells + current cell
  var N = flowAcc_raw.divide(100).add(1).rename('N_cells');
  
  // Specific catchment area As = N * cellSize (in meters)
  var As = N.multiply(cellSize_m).rename('As_m');
  
  // Moore & Burch (1986) L-factor: L = (As / 22.13)^m
  var m = 0.4;
  var L_factor = As.expression(
    'As > 1 ? pow(As / 22.13, m) : 1',
    {
      'As': As,
      'm': m
    }
  ).rename('L_moore_burch');
  
  // S-factor: McCool et al. (1987) as used in RUSLE
  var S_factor = slope_rad.expression(
    's > 0.0896 ? (16.8 * sin(s) - 0.50) : (10.8 * sin(s) + 0.03)', 
    {s: slope_rad}
  ).rename('S_rusle');
  
  var LS_scientific = L_factor.multiply(S_factor).rename('LS_moore_burch');
  
  return LS_scientific.addBands([L_factor, S_factor, As, N])
    .set('method', 'Moore_Burch_1986_HydroSHEDS_Corrected');
}

// ==========================
// P-FACTOR: UNIFORM FOR ALL LAND COVER CLASSES (SLOPE-BASED ONLY)
// ==========================
function computePfactorScientific(lcMerged, slopePercent, cFactorImage, climateZoneName) {
  // P-factor based solely on slope (Wischmeier & Smith 1978)
  // This applies uniformly to all land cover classes.
  var P_scientific = slopePercent.expression(
    '(s < 2) ? 0.10 : ' +
    '(s < 5) ? 0.12 : ' +
    '(s < 8) ? 0.14 : ' +
    '(s < 12) ? 0.25 : ' +
    '(s < 16) ? 0.30 : ' +
    '(s < 20) ? 0.50 : 1.00',
    {s: slopePercent.divide(100)}  // Convert percent to fraction
  ).rename('P_scientific').clip(lcMerged.geometry());  // Ensure clipping to AOI
  
  return P_scientific.addBands([
    P_scientific.rename('P_slope_only')
  ]).set('method', 'Wischmeier_1978_slope_only_uniform');
}

// ==========================
// COMPUTE ALL FACTORS (FIRST YEAR DISPLAY)
// ==========================
print('Computing R-factor for ' + startYear + ' ...');
var R_first = computeRFactor(startYear, zoneGeom);

print('Computing hemisphere-aware C-factor for ' + startYear + ' ...');
var cFactor_first = computeCfactorScientific(startYear, zoneGeom);

print('Computing SCIENTIFIC K-factor (USDA NRCS validated) ...');
var K_first = computeKFactorScientific(zoneGeom, analysisScale);

var lcMergedFirst = reclassifyToMerged(getLandCoverImage(startYear));
var ls_first = computeLSFactorScientific(zoneGeom, analysisScale);
var slopePercent = ls_first.select('S_rusle').multiply(100).rename('slope_percent');
var P_first = computePfactorScientific(lcMergedFirst, slopePercent, cFactor_first, climateZoneName);

// ==========================
// MAP VISUALIZATION - ALL FACTORS
// ==========================
Map.addLayer(R_first, {min: 0, max: 800, palette: ['lightblue','blue','yellow','orange','red','darkred']}, 'R-factor ' + startYear);

Map.addLayer(K_first.select('K'), {min: 0.1, max: 0.5, palette: ['green','yellow','orange','red']}, 'K-factor SCIENTIFIC (NRCS) ' + startYear);

Map.addLayer(K_first.select('texture_class'), {min: 1, max: 12, palette: [
  'darkblue','blue','lightblue','green','yellow','orange','red','magenta','cyan','purple','brown','black'
]}, 'USDA Texture Classes');

Map.addLayer(cFactor_first.select('C_scientific'), {min: 0.001, max: 0.5, palette: ['red','orange','yellow','lightgreen','darkgreen']}, 'C-factor ' + startYear + ' (Van der Knijff)');

Map.addLayer(ls_first.select('LS_moore_burch'), {min: 0, max: 20, palette: ['green','yellow','orange','red','purple']}, 'LS-factor Moore&Burch CORRECTED ' + startYear);

Map.addLayer(P_first.select('P_scientific'), {min: 0.0, max: 1.0, palette: ['darkgreen','green','yellow','orange','red']}, 'P-factor UNIFORM (slope-based) ' + startYear);

// ==========================
// FINAL RUSLE CALCULATION
// ==========================
var rusle_first = R_first.multiply(K_first.select('K'))
  .multiply(ls_first.select('LS_moore_burch'))
  .multiply(cFactor_first.select('C_scientific'))
  .multiply(P_first.select('P_scientific'))
  .rename('RUSLE_ton_ha');

Map.addLayer(rusle_first, {min: 0, max: 50, palette: ['green','yellow','orange','red','purple','black']}, 'COMPLETE RUSLE Soil Loss ' + startYear + ' (t/ha/yr)');

// ==========================
// SANITY CHECKS: SLOPE vs RUSLE & NDVI vs RUSLE
// ==========================
print('Starting model sanity checks...');

// --- Prepare necessary layers ---
var soilLoss = rusle_first;  // RUSLE in t/ha/yr
var slopePct = ls_first.select('S_rusle').multiply(100).rename('slope_percent');

// Extract NDVI (annual) from C-factor image
var ndviRaw = cFactor_first.select('annual_NDVI_raw');
var ndvi = ndviRaw.multiply(0.0001).rename('NDVI'); // MODIS NDVI scaling factor

// ------------------------------------------------------------------
// 1. Slope classes and mean RUSLE
// ------------------------------------------------------------------
// Define slope breakpoints (percent) as ee.List
var slopeBreaks = ee.List([0, 2, 5, 8, 12, 16, 20, 100]);
var slopeClassNames = ee.List(['0-2%', '2-5%', '5-8%', '8-12%', '12-16%', '16-20%', '>20%']);

// Create a slope class image using iterate (server-side safe)
var slopeClass = ee.Image(
  ee.List.sequence(0, slopeBreaks.length().subtract(2))
    .iterate(function(i, img) {
      i = ee.Number(i);
      var lower = slopeBreaks.get(i);
      var upper = slopeBreaks.get(i.add(1));
      var mask = slopePct.gte(ee.Number(lower)).and(slopePct.lt(ee.Number(upper)));
      return ee.Image(img).where(mask, i.add(1));
    }, ee.Image(0))
).rename('slope_class').toInt();

// Mask out unclassified pixels (class = 0)
slopeClass = slopeClass.updateMask(slopeClass.gt(0));

// Compute mean RUSLE per slope class
var slopeStats = soilLoss.addBands(slopeClass).reduceRegion({
  reducer: ee.Reducer.mean().group({
    groupField: 1,        // group by band index 1 (slope_class)
    groupName: 'slope_class'
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9,
  bestEffort: true
});

// Extract and print results as text
var slopeGroups = ee.List(slopeStats.get('groups'));
print('Sanity Check #1: Mean Soil Loss by Slope Class (should increase with slope)', 
  slopeGroups.map(function(g) {
    var dict = ee.Dictionary(g);
    var classNum = ee.Number(dict.get('slope_class')).int();
    var mean = dict.get('mean', 0);  // default 0 if key missing
    var className = ee.String(slopeClassNames.get(classNum.subtract(1)));
    return ee.Dictionary({
      'Slope Class': className,
      'Mean RUSLE (t/ha/yr)': ee.Number(mean).format('%.2f')
    });
  })
);

// Create a feature collection for charting
var slopeChartFeatures = ee.FeatureCollection(
  slopeGroups.map(function(g) {
    var dict = ee.Dictionary(g);
    var classNum = ee.Number(dict.get('slope_class')).int();
    var mean = dict.get('mean', 0);
    var className = ee.String(slopeClassNames.get(classNum.subtract(1)));
    return ee.Feature(null, {
      'slope_class': className,
      'mean_erosion': mean
    });
  })
);

// Display chart
var slopeChart = ui.Chart.feature.byFeature(slopeChartFeatures, 'slope_class', 'mean_erosion')
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Mean Soil Loss by Slope Class',
    hAxis: {title: 'Slope Class'},
    vAxis: {title: 'Mean RUSLE (t/ha/yr)'},
    legend: {position: 'none'},
    colors: ['blue']
  });
print('Slope Class vs. Mean Soil Loss', slopeChart);

// ------------------------------------------------------------------
// 2. NDVI classes and mean RUSLE
// ------------------------------------------------------------------
// Define NDVI breakpoints (from -1 to 1)
var ndviBreaks = ee.List([-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0]);
var ndviClassNames = ee.List(['<0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6', '0.6-0.7', '>0.7']);

// Create NDVI class image using iterate
var ndviClass = ee.Image(
  ee.List.sequence(0, ndviBreaks.length().subtract(2))
    .iterate(function(i, img) {
      i = ee.Number(i);
      var low = ndviBreaks.get(i);
      var high = ndviBreaks.get(i.add(1));
      var mask = ndvi.gte(ee.Number(low)).and(ndvi.lt(ee.Number(high)));
      return ee.Image(img).where(mask, i.add(1));
    }, ee.Image(0))
).rename('ndvi_class').toInt();

// Mask out unclassified pixels (class = 0)
ndviClass = ndviClass.updateMask(ndviClass.gt(0));

// Compute mean RUSLE per NDVI class
var ndviStats = soilLoss.addBands(ndviClass).reduceRegion({
  reducer: ee.Reducer.mean().group({
    groupField: 1,
    groupName: 'ndvi_class'
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9,
  bestEffort: true
});

var ndviGroups = ee.List(ndviStats.get('groups'));
print('Sanity Check #2: Mean Soil Loss by NDVI Class (should decrease with increasing NDVI)', 
  ndviGroups.map(function(g) {
    var dict = ee.Dictionary(g);
    var classNum = ee.Number(dict.get('ndvi_class')).int();
    var mean = dict.get('mean', 0);
    var className = ee.String(ndviClassNames.get(classNum.subtract(1)));
    return ee.Dictionary({
      'NDVI Class': className,
      'Mean RUSLE (t/ha/yr)': ee.Number(mean).format('%.2f')
    });
  })
);

// Feature collection for chart
var ndviChartFeatures = ee.FeatureCollection(
  ndviGroups.map(function(g) {
    var dict = ee.Dictionary(g);
    var classNum = ee.Number(dict.get('ndvi_class')).int();
    var mean = dict.get('mean', 0);
    var className = ee.String(ndviClassNames.get(classNum.subtract(1)));
    return ee.Feature(null, {
      'ndvi_class': className,
      'mean_erosion': mean
    });
  })
);

var ndviChart = ui.Chart.feature.byFeature(ndviChartFeatures, 'ndvi_class', 'mean_erosion')
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Mean Soil Loss by NDVI Class',
    hAxis: {title: 'NDVI Class'},
    vAxis: {title: 'Mean RUSLE (t/ha/yr)'},
    legend: {position: 'none'},
    colors: ['green']
  });
print('NDVI Class vs. Mean Soil Loss', ndviChart);

print('Sanity checks completed. If both trends are as expected, the model is behaving logically.');

// ==========================
// STATISTICS COLLECTION (with standard deviation for all factors)
// ==========================
var Rstats = R_first.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

var Kstats = K_first.select('K').reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

var cStats = cFactor_first.select('C_scientific').reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

var LSstats = ls_first.select('LS_moore_burch').reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

var Pstats = P_first.select('P_scientific').reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

var rusleStats = rusle_first.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.median(),
      sharedInputs: true
    }).combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: zoneGeom,
  scale: analysisScale,
  maxPixels: 1e9
});

// Hemisphere detection
var centroidGeom = zoneGeom.bounds().centroid({maxError: 100});
centroidGeom.evaluate(function(c) {
  var lat = c.coordinates[1];
  print('Zone centroid latitude (deg):', lat);
  print('Hemisphere:', lat < 0 ? 'Southern' : 'Northern');
});

// ==========================
// LAND COVER EROSION STATISTICS (MULTI-YEAR)
// ==========================
function safeNumber(value, defaultVal) {
  value = ee.Algorithms.If(value, value, defaultVal);
  return ee.Number(value);
}

function getLandCoverMask(lcImage, className) {
  var classIndex = ee.List(mergedClassNames).indexOf(className).add(1);
  return lcImage.eq(ee.Image.constant(classIndex));
}

function getLandCoverErosionStats(soilLossImage, lcImage, year, landcoverClass, climateZoneName) {
  var pixelArea = ee.Image.pixelArea();
  var classMask = getLandCoverMask(lcImage, landcoverClass);
  
  var pixelCount = classMask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: zoneGeom,
    scale: analysisScale,
    maxPixels: 1e10
  }).values().get(0);
  var pixelCountNum = ee.Number(pixelCount);
  
  var total_class_area_m2 = ee.Number(0);
  var mean_erosion = ee.Number(0);
  var hasPixels = pixelCountNum.gt(0);
  
  var areaDict = ee.Algorithms.If(hasPixels, pixelArea.updateMask(classMask).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: zoneGeom,
    scale: analysisScale,
    maxPixels: 1e10
  }), ee.Dictionary({area: 0}));
  total_class_area_m2 = safeNumber(ee.Dictionary(areaDict).get('area'), 0);
  
  var meanDict = ee.Algorithms.If(hasPixels, soilLossImage.updateMask(classMask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: zoneGeom,
    scale: analysisScale,
    maxPixels: 1e10
  }), ee.Dictionary({RUSLE_ton_ha: 0}));
  var mean_erosion_raw = safeNumber(ee.Dictionary(meanDict).get('RUSLE_ton_ha'), 0);
  
  mean_erosion = ee.Number(ee.Algorithms.If(
    ee.List(['Water', 'SnowIce', 'Wetlands', 'Urban']).contains(landcoverClass), 0, mean_erosion_raw
  ));
  
  var total_class_area_km2 = total_class_area_m2.divide(1e6);
  var total_class_area_ha = total_class_area_m2.divide(1e4);
  var total_erosion = mean_erosion.multiply(total_class_area_ha);
  
  return ee.Feature(null, {
    year: year,
    climate_zone: climateZoneName,
    landcoverClass: landcoverClass,
    pixel_count: pixelCountNum,
    total_class_area_km2: total_class_area_km2,
    total_class_area_ha: total_class_area_ha,
    mean_erosion_ton_ha_yr: mean_erosion,
    total_erosion_ton_yr: total_erosion
  });
}

// ==========================
// MULTI-YEAR RUSLE COMPUTATION
// ==========================
print('Computing multi-year LAND COVER erosion statistics for years: ' + startYear + ' to ' + endYear + ' ...');

function computeRUSLEForYear(year) {
  year = ee.Number(year).toInt();
  
  var lcImage = getLandCoverImage(year);
  var lcMergedYear = reclassifyToMerged(lcImage);
  
  var R = computeRFactor(year, zoneGeom);
  var K = computeKFactorScientific(zoneGeom, analysisScale).select('K');
  var C = computeCfactorScientific(year, zoneGeom);
  var LS = computeLSFactorScientific(zoneGeom, analysisScale);
  var slopePercent = LS.select('S_rusle').multiply(100).rename('slope_percent');
  var P = computePfactorScientific(lcMergedYear, slopePercent, C, climateZoneName);
  
  var soilLoss = R.multiply(K)
                  .multiply(LS.select('LS_moore_burch'))
                  .multiply(C.select('C_scientific'))
                  .multiply(P.select('P_scientific'))
                  .rename('RUSLE_ton_ha');
  
  return {soilLoss: soilLoss, lcMerged: lcMergedYear};
}

var allYearStats = ee.FeatureCollection(
  years.map(function(year) {
    var result = computeRUSLEForYear(year);
    var soilLoss = result.soilLoss;
    var lcMergedYear = result.lcMerged;
    
    var statsForYear = ee.List(mergedClassNames).map(function(landcoverClass) {
      landcoverClass = ee.String(landcoverClass);
      return getLandCoverErosionStats(soilLoss, lcMergedYear, year, landcoverClass, climateZoneName);
    });
    return statsForYear;
  }).flatten()
);

print('MULTI-YEAR LAND COVER EROSION RESULTS:', allYearStats.limit(10));

// ==========================
// VISUALIZATION CHARTS
// ==========================
var areaChartGrouped = ui.Chart.feature.groups({
  features: allYearStats,
  xProperty: 'year',
  yProperty: 'total_class_area_km2',
  seriesProperty: 'landcoverClass'
}).setChartType('ColumnChart').setOptions({
  title: 'Land Cover Area by Class (km²) - ' + climateZoneName + ' (' + startYear + '–' + endYear + ')',
  hAxis: {title: 'Year'},
  vAxis: {title: 'Area (km²)'},
  legend: {position: 'right'}
});
print('AREA CHART (grouped):', areaChartGrouped);

var meanErosionGrouped = ui.Chart.feature.groups({
  features: allYearStats,
  xProperty: 'year',
  yProperty: 'mean_erosion_ton_ha_yr',
  seriesProperty: 'landcoverClass'
}).setChartType('LineChart').setOptions({
  title: 'Mean Soil Erosion Rate by Land Cover (t/ha/yr) - ' + climateZoneName + ' (' + startYear + '–' + endYear + ')',
  hAxis: {title: 'Year'},
  vAxis: {title: 'Mean Erosion (t/ha/yr)'},
  lineWidth: 2,
  pointSize: 3,
  legend: {position: 'right'}
});
print('EROSION RATE CHART (grouped):', meanErosionGrouped);

// ==========================
// FINAL SUMMARY TABLE (with standard deviation for all)
// ==========================
print('================================================================================');
print('COMPLETE SCIENTIFIC RUSLE FACTOR SUMMARY (YEAR: ' + startYear + ') - ' + climateZoneName);
print('================================================================================');
print('R-FACTOR (CHIRPS, Renard 1997):', 'Min:', Rstats.get('R_min'), 'Max:', Rstats.get('R_max'), 'Mean:', Rstats.get('R_mean'), 'StdDev:', Rstats.get('R_stdDev'));
print('K-FACTOR (USDA NRCS, Renard 1997):', 'Min:', Kstats.get('K_min'), 'Max:', Kstats.get('K_max'), 'Mean:', Kstats.get('K_mean'), 'StdDev:', Kstats.get('K_stdDev'));
print('LS-FACTOR (Moore & Burch 1986 - CORRECTED):', 'Min:', LSstats.get('LS_moore_burch_min'), 'Max:', LSstats.get('LS_moore_burch_max'), 'Mean:', LSstats.get('LS_moore_burch_mean'), 'StdDev:', LSstats.get('LS_moore_burch_stdDev'));
print('C-FACTOR (Van der Knijff 2000):', 'Min:', cStats.get('C_scientific_min'), 'Max:', cStats.get('C_scientific_max'), 'Mean:', cStats.get('C_scientific_mean'), 'StdDev:', cStats.get('C_scientific_stdDev'));
print('P-FACTOR (Wischmeier 1978, slope-based, uniform):', 'Min:', Pstats.get('P_scientific_min'), 'Max:', Pstats.get('P_scientific_max'), 'Mean:', Pstats.get('P_scientific_mean'), 'StdDev:', Pstats.get('P_scientific_stdDev'));
print('RUSLE TOTAL SOIL LOSS:', 'Min:', rusleStats.get('RUSLE_ton_ha_min'), '| Max:', rusleStats.get('RUSLE_ton_ha_max'), '| Mean:', rusleStats.get('RUSLE_ton_ha_mean'), '| Median:', rusleStats.get('RUSLE_ton_ha_median'), '| StdDev:', rusleStats.get('RUSLE_ton_ha_stdDev'));
print('================================================================================');
print('SCIENTIFIC METHOD REFERENCES:');
print('K:  Renard et al. (1997) USDA-ARS Handbook 703');
print('C:  Van der Knijff et al. (2000) NDVI-based');
print('LS: Moore & Burch (1986) Soil Sci. Soc. Am. J. 50:1294-1298 + HydroSHEDS corrected v2.1');
print('P:  Wischmeier & Smith (1978) USDA Handbook 537 (slope-based, uniform across all land cover)');
print('R:  Renard et al. (1997) CHIRPS implementation');
print('================================================================================');