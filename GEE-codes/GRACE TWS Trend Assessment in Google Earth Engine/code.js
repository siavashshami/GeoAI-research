// ==========================
var START_DATE = ee.Date('2009-01-01');
var END_DATE   = ee.Date('2010-01-01');        // exclusive
var BAND_NAME  = 'lwe_thickness_csr';
var EXPORT_CRS = 'EPSG:4326';                  // change to e.g. EPSG:32639 (Iran)
var EXPORT_SCALE = 1000;                       // metres — realistic for GRACE
var ASSET_ID   = 'users/YOUR_USERNAME/GRACE_mean_TWS_2009'; // <-- edit this
var MAP_PALETTE = ['red', 'yellow', 'green', 'cyan', 'purple', 'blue', 'black'];


// ---------- 1. REGION OF INTEREST ------------------------------------------
var roi = table.geometry();

Map.centerObject(roi, 7);
Map.addLayer(
  table.style({color: 'red', fillColor: '00000000', width: 2}),
  {}, 'ROI'
);


// ---------- 2. LOAD GRACE COLLECTION ---------------------------------------
var grace = ee.ImageCollection('NASA/GRACE/MASS_GRIDS/LAND')
  .filterDate(START_DATE, END_DATE)
  .select(BAND_NAME);

print('GRACE collection (filtered):', grace);
print('Number of monthly images :', grace.size());
print('First image metadata     :', ee.Image(grace.first()));


// ---------- 3. MONTHLY TIME SERIES (native resolution) ---------------------
var monthlyChart = ui.Chart.image.series({
  imageCollection: grace,
  region          : roi,
  reducer         : ee.Reducer.mean(),
  scale           : 1000,
  xProperty       : 'system:time_start'
})
.setOptions({
  title             : 'GRACE Terrestrial Water Storage Anomaly (lwe_thickness_csr)',
  hAxis             : {title: 'Date'},
  vAxis             : {title: 'Equivalent water thickness (cm)'},
  lineWidth         : 2,
  pointSize         : 3,
  interpolateNulls  : true
});
print('Monthly TWS time series:', monthlyChart);


// ---------- 4. PIXEL-WISE LINEAR TREND (cm / year) -------------------------
// Convert time to "years since START_DATE" and fit y = a + b*t per pixel.
var withTime = grace.map(function(img) {
  var t = ee.Date(img.get('system:time_start'))
            .difference(START_DATE, 'year');
  return img
    .addBands(ee.Image.constant(t).rename('t').toFloat())
    .copyProperties(img, ['system:time_start']);
});

// reduce(linearFit) -> bands: 'scale' (slope) and 'offset' (intercept)
var trend  = withTime.select(['t', BAND_NAME])
                     .reduce(ee.Reducer.linearFit());
var slope  = trend.select('scale').rename('tws_trend_cm_per_yr');

print('TWS trend image (cm/yr):', slope);


// ---------- 5. MAPS ---------------------------------------------------------
// 5a. Mean TWS over the period, clipped to the ROI
var meanTWS = grace.mean().clip(roi);

Map.addLayer(
  meanTWS,
  {min: -8.554, max: -1.054, palette: MAP_PALETTE},
  'Mean TWS (cm)'
);

// 5b. Linear trend, clipped to the ROI
Map.addLayer(
  slope.clip(roi),
  {min: -5, max: 5, palette: ['red', 'white', 'blue']},
  'TWS trend (cm/yr)'
);


// ---------- 6. EXPORT MEAN IMAGE TO AN EARTH ENGINE ASSET ------------------
Export.image.toAsset({
  image       : meanTWS,
  description : 'GRACE_mean_TWS_2009',
  assetId     : ASSET_ID,
  region      : roi,
  scale       : EXPORT_SCALE,
  crs         : EXPORT_CRS,
  maxPixels   : 1e13
});